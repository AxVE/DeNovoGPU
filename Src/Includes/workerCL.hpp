#ifndef WORKERCL_HPP
#define WORKERCL_HPP

#define __CL_DEFINE_EXCEPTIONS //Enable OpenCL's exceptions
#include <CL/cl.hpp> //version 1.2

#include "log.hpp"


class WorkerCL {
	public:
		WorkerCL(size_t platform_id=0, size_t device_id=0);
		~WorkerCL();
		static void list_infos(Log& output);

	private:
		//CL bases
		cl::Platform m_platform;
		cl::Device m_device;
		cl::Context m_context;
		cl::Program::Sources m_sources;
		cl::Program m_program;
};


#endif
