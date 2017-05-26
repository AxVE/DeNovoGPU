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
	string txt;

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

	//Get device infos
	m_device_name = m_device.getInfo<CL_DEVICE_NAME>();
	txt = "gpu_name = "+m_device_name;
	m_log->write(txt);
	m_device_global_bytes = m_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
	txt = "gpu_global_mem_size = "+to_string(m_device_global_bytes)+"B";
	m_log->write(txt);
	m_device_local_bytes = m_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
	txt = "gpu_local_mem_size = "+to_string(m_device_local_bytes)+"B";
	m_log->write(txt);
	m_device_max_workgroupsize = m_device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
	txt = "gpu_max_workgroup_size = "+to_string(m_device_max_workgroupsize);
	m_log->write(txt);

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

	string txt = ""; //Used to create then output messages
	size_t buf_size=0; //Usr to calculate each buffer mem size
	size_t bufferGlobalUsage = 0; //Use to know the total of buffer usage
	size_t bufferLocalUsage = 0; //Use to know the total of buffer usage

	//Get list of contigs size, the total length of the
	//contigs concatenation and the longuest contig size
	m_log->write("Get contigs informations");
	uint64_t nbContigs = contigs.get_nbContigs();
	txt = "nbContigs = "+to_string(nbContigs); m_log->write(txt);
	uint64_t ultraSequenceSize = 0;
	uint64_t longuest_contig_size = 0;
	vector<uint64_t> contigs_size (nbContigs, 0);
	for(uint64_t i=0; i < nbContigs; i++){
		contigs_size[i] = contigs.get_sizeContig(i);
		ultraSequenceSize += contigs_size[i];
		if(contigs_size[i] > longuest_contig_size){longuest_contig_size = contigs_size[i];}
	}

	size_t nbGlobalElem = nbContigs*nbContigs;
	txt="global_nb_elems = "+to_string(nbGlobalElem); m_log->write(txt);
		//Important: the size of global range MUST be a multiplicative of local range.
		//So the size of local range is the greatest divisor of nbGlobalElem below or equal to work_group_size
	while(nbGlobalElem%work_group_size){work_group_size--;}
	txt="work_group_size = "+to_string(work_group_size); m_log->write(txt);
	
	txt = "buffer_1item_Size = "+to_string(longuest_contig_size);
	m_log->write(txt);

	//Prepare GPU for the run
	cl::Event ev;
		//infos buffer (64bits): number of contigs, size of the ultrasequence, size of longuest contig
		//buffer only accepts non-dynamics arrays (even of size 1)
	m_log->write("Prepare infos buffer");
	buf_size = sizeof(uint64_t)*3;
	txt = "infosBuf = "+to_string(buf_size)+"B";
	m_log->write(txt);
	uint64_t infos[3] = {nbContigs, ultraSequenceSize, longuest_contig_size};
	cl::Buffer buf_infos (m_context, CL_MEM_READ_ONLY, buf_size);
	m_commandqueue.enqueueWriteBuffer(buf_infos, CL_TRUE, 0, buf_size, &infos);
	m_kernel.setArg(0, buf_infos); //update kernel
	bufferGlobalUsage += buf_size;

		//Prepare the buffer for the results matrix (it will be 1D so an id of {x,y} is id=x+y*x_size)
		//The size of the buffer = char * x_size * y_size. Note: x_size == y_size == nb_contigs
	m_log->write("Prepare scores matrix buffer");
	buf_size = sizeof(char)*nbGlobalElem;
	txt = "matrixScoresBuf = "+to_string(buf_size)+"B";
	m_log->write(txt);
	size_t scores_size = buf_size;
	cl::Buffer buf_scores (m_context, CL_MEM_WRITE_ONLY, buf_size);
	m_kernel.setArg(1, buf_scores);
	bufferGlobalUsage += buf_size;

		//sequences sizes (array of 64bits) buffer
	m_log->write("Prepare contigs sizes buffer");
	buf_size = sizeof(uint64_t)*nbContigs;
	txt = "contigsSizesBuf = "+to_string(buf_size)+"B";
	m_log->write(txt);
	cl::Buffer buf_sizes (m_context, CL_MEM_READ_ONLY, buf_size);
	m_commandqueue.enqueueWriteBuffer(buf_sizes, CL_TRUE, 0, buf_size, &contigs_size[0]);
	m_kernel.setArg(2, buf_sizes);
	bufferGlobalUsage += buf_size;

		//ultrasequence, get each contigs sequence and add it in ultrasequence
	m_log->write("Prepare ultraSequence buffer");
	buf_size = sizeof(cl_char)*ultraSequenceSize;
	txt = "ultraSeqSize = "+to_string(buf_size)+"B";
	m_log->write(txt);
	cl_char* ultraSequence = new cl_char[ultraSequenceSize];
	uint64_t i = 0;
	for(uint64_t c=0; c < nbContigs; c++){
		string seq = contigs.get_seqContig(c);
		for(size_t j=0; j < seq.size(); j++){
				ultraSequence[i] = seq[j];
				i++;
		}
	}
	cl::Buffer buf_ultraseq (m_context, CL_MEM_READ_ONLY, buf_size);
	m_commandqueue.enqueueWriteBuffer(buf_ultraseq, CL_TRUE, 0, buf_size, ultraSequence);
	m_kernel.setArg(3, buf_ultraseq);
	bufferGlobalUsage += buf_size;

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

			The localisation (global or local) of buffer with NULL value is decided by the kernel parameter
			('global long* intarray' or 'local long* intarray').
		*/
	m_log->write("Prepare work-items buffers");
	buf_size = longuest_contig_size*work_group_size*sizeof(cl_char);
	txt= "charBufferGroup = "+to_string(buf_size)+"B";
	m_log->write(txt);
	m_kernel.setArg(4, buf_size, NULL); //Declare only the space size so it can be on local
	bufferLocalUsage += buf_size;

		// Each work item have need an int array. As it is on global, the buffer size must have nbContigs^2 elements (the number of work items in total)
	buf_size = longuest_contig_size*nbGlobalElem*sizeof(cl_long);
	txt= "intBufferGroup = "+to_string(buf_size)+"B";
	m_log->write(txt);
	m_kernel.setArg(5, buf_size, NULL); //Declare only the space size so it can be on local (but this one is too big to be on local)
	bufferGlobalUsage += buf_size;

	//Memory usage from buffer
	txt = "buffer_global_usage = "+to_string(bufferGlobalUsage)+"B";
	m_log->write(txt);
	txt = "buffer_global_percent = "+to_string(100*bufferGlobalUsage/m_device_global_bytes)+"%";
	m_log->write(txt);
	txt = "buffer_local_usage = "+to_string(bufferLocalUsage)+"B";
	m_log->write(txt);
	txt = "buffer_local_percent = "+to_string(100*bufferLocalUsage/m_device_local_bytes)+"%";
	m_log->write(txt);

	//Run the kernel and wait the end (global ids (contig1_id, contig2_id) is 2D to 1D (nbContigs*nbContigs), local is 1D (work group item id)
	//So contig2_id = global_id/nbContigs and contig1_id=global_id-contig2_id*nbContigs
	m_log->write("Run kernel");
		//Launch kernel
	m_commandqueue.enqueueNDRangeKernel(m_kernel,cl::NullRange, cl::NDRange(nbGlobalElem), cl::NDRange(work_group_size), NULL, &ev);
	ev.wait();

	//Get the scores matrix: get the buffer into a 1D array then convert de 2D vectors array
	int8_t* scores_1D = new int8_t[nbGlobalElem];
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
	long score_2_seq(local char *seq1, unsigned long seq1_size, global char *seq2, unsigned long seq2_size, global long *intarray){
		//Test: if same seq then =1 else =0
		if(seq1_size == seq2_size){
			bool same=true;
			for(size_t i=0; i < seq1_size; i++){
				if(seq1[i] != seq2[i]){
					intarray[0]=i;
					same = false;
					return -i;
					i = seq1_size;
				}
			}
			if(same){return 1;}
		}
		return 0;

	}

	kernel void cmp_2_contigs(__global unsigned long *infos, __global char *scores, __global unsigned long *seqs_sizes, __global char *ultraseq, __local char *charbufloc, __global long *intbufloc){
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
			//Get the position of the first sequence inside ultraseq
		unsigned long start = 0;
		for(unsigned long i=0; i < seq1_id; i++){start += seqs_sizes[i];}
			//Copy the contig sequence in local buffer (copy each char one by one from global to local)
		for(unsigned long c=0; c < seq1_size; c++){
			seq1[c] = ultraseq[start+c];
		}
		
		unsigned long seq1_start = start; //DEBUG

		//Get the second contig sequence and its infos. It is read only once, so there is no need to copy it in local.
			//Get size
		unsigned long seq2_size = seqs_sizes[seq2_id];
			//Calculate the begin of this seq in ultraseq
		start = 0;
		for(unsigned long i=0; i < seq2_id; i++){start += seqs_sizes[i];}
			//Get seq
		global char* seq2 = &ultraseq[start];
		
		unsigned long seq2_start = start; //DEBUG

		//Get the global int array buffer
		//global long *inta = &intbufloc[infos[2]*gid];
		global long *inta = NULL;

		//Get match score of seq2 on seq1
		scores[seq1_id+nbContigs*seq2_id]=score_2_seq(seq1, seq1_size, seq2, seq2_size, inta);
	}
)CLCODE";
