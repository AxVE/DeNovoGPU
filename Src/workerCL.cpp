#include "workerCL.hpp"


#define __CL_DEFINE_EXCEPTIONS
#include <CL/cl.hpp>

#include <vector>
#include <string>
#include <iostream>

#include "log.hpp"

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

	//CL code
	
	//Build kernel sources
	
}

WorkerCL::~WorkerCL(){
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

string WorkerCL::kernel_cmp_2_contigs = R"CLCODE(
	kernel void cmp_2_contigs(){
		int i = 0;
	}
)CLCODE";
