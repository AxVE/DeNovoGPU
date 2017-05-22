#include "workerCL.hpp"


#define __CL_DEFINE_EXCEPTIONS
#include <CL/cl.hpp>

#include <vector>
#include <string>
#include <iostream>

#include "log.hpp"
#include "reads.hpp"

using namespace std;

WorkerCL::WorkerCL(size_t platform_id, size_t device_id, Log& log){
	//Get the log output manager
	m_log = &log;	

	/* GPU environment preparation */
	
	//Get platforms (drivers) list and keep the one to be use
	vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if(!all_platforms.size()){
		//No platform ? Send error
		string error = "OPENCL: no platforms found.";
		throw(error);
	}
	m_platform = all_platforms[platform_id];

	//Get devices list and keep the one to use
	vector<cl::Device> devices;
	m_platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
	if(device_id >= devices.size()){
		string error = "OPENCL: no device of id "+to_string(device_id)+
			" on platform "+to_string(platform_id)+" ("+to_string(devices.size())+" devices)";
		throw(error);
	}
	m_device = devices[device_id];

	//Create context
	m_context = cl::Context({m_device});
	
	//initialize commands queue
	m_commandqueue = cl::CommandQueue(m_context, m_device);

	//Build kernel sources
	m_program = cl::Program(m_context, kernel_cmp_2_contigs);
	if(m_program.build({m_device}) != CL_SUCCESS){
		string error = "Error building: ";
		error += m_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_device);
		throw(error);
	}
	

	//Make the kernel
	m_kernel = cl::Kernel(m_program, "cmp_2_contigs");

}

WorkerCL::~WorkerCL(){
}

void WorkerCL::run(const Contigs& contigs, size_t work_group_size){
	/*
	Create the string containing all contigs sequences
	and store the list of contigs size (in same orders)
	to be able to retrieve contigs. Size is store in
	an uint64_t (ulong) to be sure it is 64bits as cl_ulong
	*/

	//Get list of contigs size, the total length of the
	//contigs concatenation and the longuest contig size
	uint64_t nbContigs = contigs.get_nbContigs();
	uint64_t ultraSequenceSize = 0;
	uint64_t longuest_contig_size = 0;
	vector<uint64_t> contigs_size (nbContigs, 0);
	for(uint64_t i=0; i < nbContigs; i++){
		contigs_size[i] = contigs.get_sizeContig(i);
		ultraSequenceSize += contigs_size[i];
		if(contigs_size[i] > longuest_contig_size){longuest_contig_size = contigs_size[i];}
	}


	//Prepare GPU for the run
	cl::Event ev;
		//infos buffer (64bits): number of contigs, size of the ultrasequence, size of longuest contig
		//buffer only accepts non-dynamics arrays (even of size 1)
	m_log->write("Prepare infos buffer");
	uint64_t infos[3] = {nbContigs, ultraSequenceSize, longuest_contig_size};
	cl::Buffer buf_infos (m_context, CL_MEM_READ_ONLY, sizeof(uint64_t));
	m_commandqueue.enqueueWriteBuffer(buf_infos, CL_TRUE, 0, sizeof(uint64_t), &infos);
	m_kernel.setArg(0, buf_infos); //update kernel

		//Prepare the buffer for the results matrix (it will be 1D so an id of {x,y} is id=x+y*x_size)
		//The size of the buffer = char * x_size * y_size. Note: x_size == y_size == nb_contigs
	m_log->write("Prepare scores matrix buffer");
	unsigned int scores_size = sizeof(char)*nbContigs*nbContigs;
	cl::Buffer buf_scores (m_context, CL_MEM_WRITE_ONLY, scores_size);
	m_kernel.setArg(1, buf_scores);

		//sequences sizes (array of 64bits) buffer
	m_log->write("Prepare contigs sizes buffer");
	cl::Buffer buf_sizes (m_context, CL_MEM_READ_ONLY, sizeof(uint64_t)*nbContigs);
	m_commandqueue.enqueueWriteBuffer(buf_sizes, CL_TRUE, 0, sizeof(uint64_t)*nbContigs, &contigs_size[0]);
	m_kernel.setArg(2, buf_sizes);

		//ultrasequence, get each contigs sequence and add it in ultrasequence
	m_log->write("Prepare ultraSequence buffer");
	string txt = "ultraSeqSize = "+to_string(ultraSequenceSize*sizeof(char))+"B";
	m_log->write(txt);
	char* ultraSequence = new char[ultraSequenceSize];
	uint64_t i = 0;
	for(uint64_t c=0; c < nbContigs; c++){
		string seq = contigs.get_seqContig(c);
		for(size_t j=0; j < seq.size(); j++){
				ultraSequence[i] = seq[j];
				i++;
		}
	}
	cl::Buffer buf_ultraseq (m_context, CL_MEM_READ_ONLY, sizeof(char)*ultraSequenceSize);
	m_commandqueue.enqueueWriteBuffer(buf_ultraseq, CL_TRUE, 0, sizeof(char)*ultraSequenceSize, ultraSequence);
	m_kernel.setArg(3, buf_ultraseq);

	delete ultraSequence; //Clean  the memory
	ultraSequence = nullptr;

		//buffer for work items : they need 2 arrays for need2a and 2 array for each sequences. These arrays are of size of longuest contig.
		/*
			Each work item need its own array (for each arrays) allocated on local.
			So a local array (of a work group) contains all concatenated array of the work items of the same group.
			This local array have a size of longuest_contig_size*work_group_size number of elements.
		*/

	//Run the kernel and wait the end
	m_log->write("Run kernel");
	m_commandqueue.enqueueNDRangeKernel(m_kernel,cl::NullRange, cl::NDRange(nbContigs, nbContigs), cl::NullRange, NULL, &ev);
	ev.wait();

	//Get the score matrix: get the buffer into a 1D array then convert de 2D vectors array
	m_log->write("Get scores matrix");
	char* scores_1D = new char[nbContigs*nbContigs];
	m_commandqueue.enqueueReadBuffer(buf_scores, CL_TRUE, 0, scores_size, scores_1D);
	vector< vector<char> > scores = vector< vector<char> >(nbContigs, vector<char>(nbContigs, 0));
	for(size_t j=0; j<nbContigs; j++){
		for(size_t i=0; i<nbContigs; i++){
			scores[i][j] = scores_1D[i+nbContigs*j];
		}
	}

	//TEST: display scores
	/*
	txt= "Seqs";
	for(size_t i=0; i < nbContigs; i++){txt += "\t" + to_string(i);}
	m_log->write(txt);
	for(size_t j=0; j < nbContigs; j++){
		txt = to_string(j);
		for(size_t i=0; i < nbContigs; i++){
			txt += "\t"+to_string(scores[i][j]);
		}
		m_log->write(txt);
	}
	*/

	//Clean the memory
	delete scores_1D;
	scores_1D = nullptr;
}

void WorkerCL::list_infos(Log& output){
	string txt;

	//Get platforms
	txt = "\n\t[Platforms list]";
	output.write(txt);
	vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	if(!platforms.size()){
		txt = "No platform detected";
		output.write(txt);
		return;
	}
	txt = "platform_id\tplatform_name\tplatform_profile\tplatform_vendor\tplatform_version";
	output.write(txt);
	for(size_t i=0; i < platforms.size(); i++){
		cl::Platform& p = platforms[i];
		txt = to_string(i);
		txt += "\t" + p.getInfo<CL_PLATFORM_NAME>();
		txt += "\t" + p.getInfo<CL_PLATFORM_PROFILE>();
		txt += "\t" + p.getInfo<CL_PLATFORM_VENDOR>();
		txt += "\t" + p.getInfo<CL_PLATFORM_VERSION>();
		//txt += "\t" + p.getInfo<CL_PLATFORM_EXTENSIONS>();
		output.write(txt);
	}

	//Get devices
	txt = "\n\t[Devices list]";
	output.write(txt);
	txt = "platform_id\tdevice_id\tdevice_name\tdevice_vendor\tdevice_profile\tdevice_version\tdevice_globalmem\tdevice_localmem\tdevice_maxgroupsize\tdriver_version\topencl_c_version";
	output.write(txt);
	for(size_t p = 0; p < platforms.size(); p++){
		vector<cl::Device> devices;
		platforms[p].getDevices(CL_DEVICE_TYPE_ALL, &devices);
		for(size_t d = 0; d < devices.size(); d++){
			cl::Device& device = devices[d];
			txt = to_string(p)+"\t"+to_string(d);
			txt += "\t" + device.getInfo<CL_DEVICE_NAME>();
			txt += "\t" + device.getInfo<CL_DEVICE_VENDOR>();
			txt += "\t" + device.getInfo<CL_DEVICE_PROFILE>();
			txt += "\t" + device.getInfo<CL_DEVICE_VERSION>();
			txt += "\t" + to_string(device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>())+"B";
			txt += "\t" + to_string(device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>())+"B";
			txt += "\t" + to_string(device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
			txt += "\t" + device.getInfo<CL_DRIVER_VERSION>();
			txt += "\t" + device.getInfo<CL_DEVICE_OPENCL_C_VERSION>();
			output.write(txt);
		}
	}
}	

/*
Les balises 'R"CLCODE(' et ')CLCODE'   (du type R"NAME( ... )NAME") servent à définir un
string litéral brut. C'est utile pour avoir un string sur plusieurs ligne, comme un code,
et cela évite d'avoir à ouvrir puis fermer les guillemets à chaque ligne.
*/

string WorkerCL::kernel_cmp_2_contigs = R"CLCODE(
	kernel void cmp_2_contigs(global unsigned long *infos, global char *scores, global unsigned long *seqs_sizes, global char *ultraseq){
		int seq1_id = get_global_id(0);
		int seq2_id = get_global_id(1);

		/* Old way: basic memory. Work
		//Get the position of the seqs in the ultraseq
		unsigned long seq1_begin = 0;
		unsigned long seq2_begin = 0;
		unsigned long i=0;
		while(i < seq1_id || i < seq2_id){
			if(i < seq1_id){seq1_begin += seqs_sizes[i];}
			if(i < seq2_id){seq2_begin += seqs_sizes[i];}
			i++;
		}

		//Get the sizes and the end pos of the seqs
		unsigned long seq1_size = seqs_sizes[seq1_id];
		unsigned long seq2_size = seqs_sizes[seq2_id];
		unsigned long seq1_end= seq1_begin + seq1_size - 1;
		unsigned long seq2_end= seq2_begin + seq1_size - 1;

		//Test if same seq : =1 else =0
		scores[seq1_id + infos[0]*seq2_id] = 0;
		if(seq1_size == seq2_size){
			bool same = true;
			for(unsigned int i=0; i < seq1_size; i++){
				if(ultraseq[seq1_begin+i] != ultraseq[seq2_begin+i]){
					same = false;
					i = seq1_size;
				}
			}
			if(same){scores[seq1_id+infos[0]*seq2_id]=1;}
		}
		*/

		//Get the position of seqs and sizes
		global char* seq1_ptr = NULL;
		global char* seq2_ptr = NULL;
			//Calculate the begin for each item (the for boucle is on all contigs (infos[0]) but it's stop before)
		unsigned long start = 0;
		for(unsigned long i = 0; i < seq1_id; i++){start += seqs_sizes[i];}
		seq1_ptr = &ultraseq[start];
		start = 0;
		for(unsigned long i = 0; i < seq2_id; i++){start += seqs_sizes[i];}
		seq2_ptr = &ultraseq[start];
			//get sizes
		unsigned long seq1_size = seqs_sizes[seq1_id];
		unsigned long seq2_size = seqs_sizes[seq2_id];
			
		//Test: if same seq then =1 else =0
		scores[seq1_id + infos[0]*seq2_id] = 0;
		if(seq1_size == seq2_size){
			bool same=true;
			for(unsigned int i=0; i < seq1_size; i++){
				if(seq1_ptr[i] != seq2_ptr[i]){
					same = false;
					i = seq1_size;
				}
			}
			if(same){scores[seq1_id+infos[0]*seq2_id]=1;}
		}

	}
)CLCODE";
