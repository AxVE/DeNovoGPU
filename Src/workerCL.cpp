#include "workerCL.hpp"


#define __CL_DEFINE_EXCEPTIONS
#include <CL/cl.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "log.hpp"
#include "reads.hpp"

using namespace std;

//This function return the gcd between 2 numbers
size_t gcd(size_t a, size_t b){
	//Using binary method (so can be heavily optimise)
	//Note: the gcd() in algorithm std IS actually (05/2017) NOT official ! So It's not use here
	//C++17 should include a gcd() function in <numeric>
	size_t d = 0;
	while(a%2==0 && b%2==0){
		a=a/2;
		b=b/2;
		d++;
	}
	while(a!=b){
		if(a%2==0){a=a/2;}
		else if(b%2==0){b=b/2;}
		else if(a>b){a=(a-b)/2;}
		else{b=(b-a)/2;}
	}
	return a*pow(2,d);
}

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
	/*
	txt = "itemBufferSize = "+to_string(longuest_contig_size);
	m_log->write(txt);
	txt = "workGroupSize = "+to_string(work_group_size);
	m_log->write(txt);
	size_t bufSize = longuest_contig_size*work_group_size*sizeof(char);
	txt= "charBufferGroup = "+to_string(bufSize)+"B";
	m_log->write(txt);
	cl::Buffer buf_localChar (m_context, CL_MEM_READ_WRITE, bufSize);
	m_kernel.setArg(4, buf_localChar);
	bufSize = longuest_contig_size*work_group_size*sizeof(int64_t);
	txt= "intBufferGroup = "+to_string(bufSize)+"B";
	m_log->write(txt);
	cl::Buffer buf_localInt (m_context, CL_MEM_READ_WRITE, bufSize);
	m_kernel.setArg(5, buf_localChar);
	*/

	//Run the kernel and wait the end (global ids (contig1_id, contig2_id) is 2D to 1D (nbContigs*nbContigs), local is 1D (work group item id)
	//So contig2_id = global_id/nbContigs and contig1_id=global_id-contig2_id*nbContigs
	m_log->write("Run kernel");
	size_t nbGlobalElem = nbContigs*nbContigs;
	txt="global_nb_elems="+to_string(nbGlobalElem); m_log->write(txt);
		//Important: the number of global ids MUST be a multiplicative of the number of local ids. So = greatest_common_divisor(global_ids,local_ids)
	size_t nbLocalElem = gcd(nbGlobalElem, work_group_size);
	txt="local_nb_elems="+to_string(nbLocalElem); m_log->write(txt);
		//Launch kernel
	m_commandqueue.enqueueNDRangeKernel(m_kernel,cl::NullRange, cl::NDRange(nbGlobalElem), cl::NDRange(work_group_size), NULL, &ev);
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

Reminder:
	infos[0]: number of contigs
	infos[1]: size of ultrasequence
	infos[2]: size of longuest contig (so the size of a lot of elements of buffers arrays)
*/

string WorkerCL::kernel_cmp_2_contigs = R"CLCODE(
	//kernel void cmp_2_contigs(global unsigned long *infos, global char *scores, global unsigned long *seqs_sizes, global char *ultraseq, local char *charbufloc, local long *intbufloc){
	kernel void cmp_2_contigs(global unsigned long *infos, global char *scores, global unsigned long *seqs_sizes, global char *ultraseq){
		int seq2_id = get_global_id(0)/infos[0];
		int seq1_id = get_global_id(0) - seq2_id*infos[0];
		int work_id = get_local_id(0);

		scores[seq1_id+infos[0]*seq2_id] = seq1_id*10+seq2_id;

		/*

		//Get the first contig sequence and its infos. As it will be read multiple times, copy it in local-item buffer
			//Get size
		unsigned long seq1_size = seqs_sizes[seq1_id];
			//Prepare the buffer to use
		local char* seq1_ptr = &charbufloc[infos[2]*work_id];
			//Calculate the begin for each item (the for boucle is on all contigs (infos[0]) but it's stop before (seq1_id < infos[0]))
		unsigned long start = 0;
		for(unsigned long i=0; i < seq1_id; i++){start += seqs_sizes[i];}
			//Copy the contig sequence in local buffer (copy each char one by one from global to local)
		for(unsigned long c=start; c < start+seq1_size; c++){
			seq1_ptr[c-start] = ultraseq[c];
		}

		//Get the second contig sequence and its infos. It is read only once, so there is no need to copy it in local.
			//Get size
		unsigned long seq2_size = seqs_sizes[seq2_id];
			//Prepare the buffer to use
		global char* seq2_ptr = NULL;
			//Calculate the begin of this seq in ultraseq
		start = 0;
		for(unsigned long i = 0; i < seq2_id; i++){start += seqs_sizes[i];}
			//Get seq
		seq2_ptr = &ultraseq[start];

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
		*/

	}
)CLCODE";
