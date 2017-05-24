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

vector< vector<int8_t> > WorkerCL::run(const Contigs& contigs, size_t work_group_size){
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

		//buffers for work items
		/*
			They need 1 array for need1a and 1 array for the one contig as the another one is only read once per work-item,
			there is no necessity to copy it in local.
			These arrays are of size of longuest contig.
			Each work item need its own array (for each arrays) allocated on local.
			So a local array (of a work group) contains all concatenated array of the work items of the same group.
			This local array have a size of longuest_contig_size*work_group_size number of elements.
		*/
	m_log->write("Prepare work-items buffers");
	txt = "itemBufferSize = "+to_string(longuest_contig_size);
	m_log->write(txt);
	txt = "workGroupSize = "+to_string(work_group_size);
	m_log->write(txt);
	size_t locSize = longuest_contig_size*work_group_size*sizeof(int8_t);
	txt= "charBufferGroup = "+to_string(locSize)+"B";
	m_log->write(txt);
	m_kernel.setArg(4, locSize, NULL); //Declare only the space size so it can be on local
	locSize = longuest_contig_size*work_group_size*sizeof(int64_t);
	txt= "intBufferGroup = "+to_string(locSize)+"B";
	m_log->write(txt);
	m_kernel.setArg(5, locSize, NULL); //Declare only the space size so it can be on local

	//Run the kernel and wait the end (global ids (contig1_id, contig2_id) is 2D to 1D (nbContigs*nbContigs), local is 1D (work group item id)
	//So contig2_id = global_id/nbContigs and contig1_id=global_id-contig2_id*nbContigs
	m_log->write("Run kernel");
	size_t nbGlobalElem = nbContigs*nbContigs;
	txt="global_nb_elems="+to_string(nbGlobalElem); m_log->write(txt);
		//Important: the number of global ids MUST be a multiplicative of the number of local ids. So It is the greastest divisor of nbGlobalElem below or equal to work_group_size
	size_t nbLocalElem = work_group_size; while(nbGlobalElem%nbLocalElem){nbLocalElem--;}
	txt="local_nb_elems="+to_string(nbLocalElem); m_log->write(txt);
		//Launch kernel
	m_commandqueue.enqueueNDRangeKernel(m_kernel,cl::NullRange, cl::NDRange(nbGlobalElem), cl::NDRange(nbLocalElem), NULL, &ev);
	ev.wait();

	//Get the scores matrix: get the buffer into a 1D array then convert de 2D vectors array
	m_log->write("Get scores matrix");
	int8_t* scores_1D = new int8_t[nbContigs*nbContigs];
	m_commandqueue.enqueueReadBuffer(buf_scores, CL_TRUE, 0, scores_size, scores_1D);
	vector< vector<int8_t> > scores = vector< vector<int8_t> >(nbContigs, vector<int8_t>(nbContigs, 0));
	for(size_t j=0; j<nbContigs; j++){
		for(size_t i=0; i<nbContigs; i++){
			scores[i][j] = scores_1D[i+nbContigs*j];
		}
	}

	//Clean the memory
	delete scores_1D;
	scores_1D = nullptr;

	//Return the scores matrix
	return scores;
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

Reminder:
	infos[0]: number of contigs
	infos[1]: size of ultrasequence
	infos[2]: size of longuest contig (so the size of a lot of elements of buffers arrays)
*/

string WorkerCL::kernel_cmp_2_contigs = R"CLCODE(
	//This function return the match score of seq2 on seq1. The local array buffer intarray must be (at least) of the size of seq1 (so seq1_size).
	long score_2_seq(local char *seq1, unsigned long seq1_size, global char *seq2, unsigned long seq2_size, local long *intarray){
		//Test: if same seq then =1 else =0
		if(seq1_size == seq2_size){
			bool same=true;
			for(unsigned int i=0; i < seq1_size; i++){
				if(seq1[i] != seq2[i]){
					same = false;
					i = seq1_size;
				}
			}
			if(same){return 1;}
		}
		return 0;
	}

	kernel void cmp_2_contigs(__global unsigned long *infos, __global char *scores, __global unsigned long *seqs_sizes, __global char *ultraseq, __local char *charbufloc, __local long *intbufloc){
		size_t gid = get_global_id(0);
		size_t seq2_id = gid/infos[0];
		size_t seq1_id = gid - seq2_id*infos[0];
		size_t work_id = get_local_id(0);
		unsigned long nbContigs = infos[0];

		//Get the first contig sequence and its infos. As it will be read multiple times, copy it in local-item buffer
			//Get size
		unsigned long seq1_size = seqs_sizes[seq1_id];
			//Prepare the buffer to use
		local char* seq1 = &charbufloc[infos[2]*work_id];
			//Calculate the begin for each item (the for boucle is on all contigs (nbContigs) but it's stop before (seq1_id < nbContigs))
		unsigned long start = 0;
		for(unsigned long i=0; i < seq1_id; i++){start += seqs_sizes[i];}
			//Copy the contig sequence in local buffer (copy each char one by one from global to local)
		for(unsigned long c=0; c < seq1_size; c++){
			seq1[c] = ultraseq[start+c];
		}

		//Get the second contig sequence and its infos. It is read only once, so there is no need to copy it in local.
			//Get size
		unsigned long seq2_size = seqs_sizes[seq2_id];
			//Calculate the begin of this seq in ultraseq
		start = 0;
		for(unsigned long i = 0; i < seq2_id; i++){start += seqs_sizes[i];}
			//Get seq
		global char* seq2 = &ultraseq[start];

		//Get the local int array buffer
		local long *inta = &intbufloc[infos[2]*work_id];

		//Get match score of seq2 on seq1
		scores[seq1_id+nbContigs*seq2_id]=score_2_seq(seq1, seq1_size, seq2, seq2_size, inta);
	}
)CLCODE";
