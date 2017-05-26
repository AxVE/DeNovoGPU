#ifndef WORKERCL_HPP
#define WORKERCL_HPP

#define __CL_DEFINE_EXCEPTIONS //Enable OpenCL's exceptions
#include <CL/cl.hpp> //version 1.2

#include "log.hpp"
#include "reads.hpp"


class WorkerCL {
	public:
		WorkerCL(size_t platform_id, size_t device_id, Log& log);
		~WorkerCL();

		//Return the matrix of coupling scores of a contigs set
		std::vector< std::vector<int8_t> > run(const Contigs& contigs, size_t work_group_size=1);
		
		static void list_infos(Log& output);

	private:
		//CL bases
		cl::Platform m_platform;
		cl::Device m_device;
		cl::Context m_context;
		cl::CommandQueue m_commandqueue; //To send order/demand to gpu
		//cl::Program::Sources m_sources;
		cl::Program m_program;
		cl::Kernel m_kernel;

		//device infos
		std::string m_device_name="Unknown";
		size_t m_device_global_bytes=1;
		size_t m_device_local_bytes=1;
		size_t m_device_max_workgroupsize=1;
		size_t m_device_max_mem_alloc_size=1;
		
		//CL kernel
		static std::string kernel_cmp_2_contigs;

		//Ouputs infos manager
		Log* m_log = nullptr;
};


#endif
