#include "workerCL.hpp"


#define __CL_DEFINE_EXCEPTIONS
#include <CL/cl.hpp>

#include <vector>
#include <string>
#include <iostream>

#include "log.hpp"
#include "reads.hpp"

using namespace std;

WorkerCL::WorkerCL(size_t platform_id, size_t device_id){
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

void WorkerCL::run(Contigs contigs){
	/*
	Create the string containing all contigs sequences
	and store the list of contigs size (in same orders)
	to be able to retrieve contigs. Size is store in
	an uint64_t (ulong) to be sure it is 64bits as cl_ulong
	*/

	//Get list of contigs size and the total length of the
	//contigs concatenation
	uint64_t nbContigs = contigs.get_nbContigs();
	uint64_t ultraSequenceSize = 0;
	vector<uint64_t> contigs_size (nbContigs, 0);
	for(uint64_t i=0; i < nbContigs; i++){
		contigs_size[i] = contigs.get_sizeContig(i);
		ultraSequenceSize += contigs_size[i];
	}

	//Create the ultraSequence
	char* ultraSequence = new char[ultraSequenceSize];
	uint64_t i = 0;
		//Get each contigs sequence and add it in ultraSequence
	for(uint64_t c=0; c < nbContigs; c++){
		string seq = contigs.get_seqContig(c);
		for(size_t j=0; j < seq.size();j++){
			ultraSequence[i] = seq[j];
			i++;
		}
	}
	cout << "UltraSequence:" << endl;
	cout << ultraSequence << endl;

	//Prepare GPU for the run
	cl::Event ev;
		//infos buffer (64bits): nbcontigs
		//buffer only accepts non-dynamics arrays (even of size 1)
	uint64_t infos[1] = {ultraSequenceSize};
	cl::Buffer buf_infos (m_context, CL_MEM_READ_ONLY, sizeof(uint64_t));
	m_commandqueue.enqueueWriteBuffer(buf_infos, CL_TRUE, 0, sizeof(uint64_t), &infos);
		//sequences sizes (array of 64bits) buffer
	cl::Buffer buf_sizes (m_context, CL_ME_READ_ONLY, sizeof(uint64_t)*nbContigs);
	m_commandqueue.enqueueWriteBuffer(buf_sizes, CL_TRUE, 0, sizeof(uint64_t)*nbContigs, &contigs_size[0]);

	//Update the kernel (gpu function)
	m_kernel.setArg(0, buf_infos);

	//Run the kernel and wait the end
	m_commandqueue.enqueueNDRangeKernel(m_kernel,cl::NullRange, cl::NDRange::NDRange(nbContigs, nbContigs), cl::NullRange, NULL, &ev);
	ev.wait();
	

	//Clean the memory
	delete ultraSequence;

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
	txt = "platform_id\tdevice_id\tdevice_name\tdevice_vendor\tdevice_profile\tdevice_version\tdriver_version\topencl_c_version";
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
	kernel void cmp_2_contigs(){
		int i = 0;
	}
)CLCODE";
