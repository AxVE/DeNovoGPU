#ifndef WORKERCL_HPP
#define WORKERCL_HPP

#define __CL_DEFINE_EXCEPTIONS //Enable OpenCL's exceptions
#include <CL/cl.hpp> //version 1.2

#include "log.hpp"
#include "reads.hpp"


class WorkerCL {
	public:
		WorkerCL(size_t platform_id=0, size_t device_id=0);
		~WorkerCL();

		//Return the matrix of coupling scores of a contigs set
		void run(Contigs contigs);
		
		static void list_infos(Log& output);

	private:
		//CL bases
		cl::Platform m_platform;
		cl::Device m_device;
		cl::Context m_context;
		//cl::Program::Sources m_sources;
		cl::Program m_program;

		//CL kernel
		static std::string kernel_cmp_2_contigs;
};


#endif
