#ifndef CMDPARSER_HPP
#define CMDPARSER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <thread>

#include "optionparser.h" // http://www.optionparser.sourceforge.net

using namespace std;

//Parameters storage
struct Params{
	string readsFile="required";
	string logs="cout";
	int8_t requiredScore=60;
	unsigned nbthreads=1;
	bool gpu_infos=false;
	bool gpu=false;
	size_t opencl_platform_id = 0;
	size_t opencl_device_id = 0;

	
	/*
	string outPrefix="Results/";
	size_t nb_iteration=1000;
	size_t size_x=10;
	size_t size_y=10;
	size_t size_z=10;
	size_t write_every=1;
	int parallel_nb_threads=1;
	size_t jobs_groupsize=1;
	*/

	string infosToString() const {
		string s;
		s += "readsFile="+readsFile;
		s += "\nlogs="+logs;
		s += "\nrequiredScore="+to_string(requiredScore);
		s += "\nnbthreads="+to_string(nbthreads);
		s += "\ngpu="+to_string(gpu);
		if(gpu){
			s+= "\nopencl_platform_id="+to_string(opencl_platform_id);
			s+= "\nopencl_device_id="+to_string(opencl_device_id);
		}
		return s;
	}
};

//Define arguments types (getting the default options)
struct Arg: public option::Arg
{
	static option::ArgStatus Unknown(const option::Option& option, bool msg){
		if(msg) cerr << "Unknown option '" << option.name << "'." << endl;
		return option::ARG_ILLEGAL;
	}
	static option::ArgStatus Required(const option::Option& option, bool msg){
		if (option.arg != 0) return option::ARG_OK;
		if(msg) cerr << "Option '" << option.name << "' requires an argument." << endl;
		return option::ARG_ILLEGAL;
	}
	
	static option::ArgStatus Numeric(const option::Option& option, bool msg){
		char* endptr = 0;
		if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
		if (endptr != option.arg && *endptr == 0) return option::ARG_OK;

		if (msg) cout << "Option '" << option.name << "' requires a numeric argument." << endl;
		return option::ARG_ILLEGAL;
	}
	static option::ArgStatus NonEmpty(const option::Option& option, bool msg){
		if (option.arg != 0 && option.arg[0] != 0) return option::ARG_OK;

		if (msg) cerr << "Option '" << option.name << "' requires a non-empty argument" << endl;;
			return option::ARG_ILLEGAL;
	}
};

//Define options list
//enum optionIndex { UNKNOWN, HELP, VERSION, POLYFILE, OUTPUT, ENVX, ENVY, ENVZ, NBITER, WRITEACH, NB_THREAD, PACKAGE_SIZE, GPU };
enum optionIndex { UNKNOWN, HELP, VERSION, READSFILE, LOGS, REQUIREDSCORE, NBTHREADS, GPU, GPUINFOS, CLPLATFORMID, CLDEVICEID};

//Define parser
const option::Descriptor usage[] = {
	{UNKNOWN, 0, "", "", Arg::Unknown,		"DeNovo reads alignment using GPU.\nDeveloper: Verdier Axel (axel.l.verdier@free.fr)\n"},
	{UNKNOWN, 0, "", "", Arg::Unknown, 		"Usage: ./denovoGPU -f <path> [options]"},
	{UNKNOWN, 0, "", "", Arg::Unknown, 		"\nBasic options:"},
	{HELP, 0, "h", "help", option::Arg::None, 	"     -h --help  \tPrint usage and exit."},
	{VERSION, 0, "V", "version", option::Arg::None,	"     -V --version  \tPrint version."},
	{UNKNOWN, 0, "", "", Arg::Unknown, 		"\nSimulation options:"},
	{READSFILE, 0, "f", "file", Arg::Required,	"     -f <path> --file <path>  \tFile of the reads. Required."},
	{LOGS, 0, "l", "log", Arg::NonEmpty,	"     -l <path> --logs <path>  \tFile of the logs. Note: 'cout' or 'stdout' will make logs to be send to the standard output. Default: cout."},
	{REQUIREDSCORE, 0, "s", "required-score", Arg::Numeric,	"     -s <numeric> --required-score <numeric> \tThe absolute value of the minimal required score of allignment between 2 contigs to merge them. Must be ]0;100]. Default: 60."}, 
	{NBTHREADS, 0, "t", "nbThreads", Arg::Numeric,	"     -t <numeric> --nbThreads <numeric> \tThe number of threads (main+additionals) to use. Must be > 0. Default: 1."},
	{UNKNOWN, 0, "", "", Arg::Unknown, 		"\nOpenCL options:"},
	{GPUINFOS, 0, "i", "opencl_infos", Arg::None, "     -i --opencl_infos  \tOutput OpenCL informations."},
	{GPU, 0, "g", "gpu_use", Arg::None, "     -g --gpu_use  \tSwitch on OpenCL (GPU) matrix calculation instead on CPU."},
	{CLPLATFORMID, 0, "p", "opencl_platform_id", Arg::Numeric, "     -p <numeric> --opencl_platform_id <numeric> \tThe id of the platform to use. Get list with --opencl_infos. Default: 0."},
	{CLDEVICEID, 0, "d", "opencl_device_id", Arg::Numeric, "     -d <numeric> --opencl_device_id <numeric> \tThe id of the device to use. Get list with --opencl_infos. Default: 0."},
	/*
	{OUTPUT, 0, "o", "output", Arg::NonEmpty,	"     -o <path>/ --output <path>/  \tThe output prefix. A directory (with '/' at the end) is advice. Default: \"Results/\"."},
	{ENVX, 0, "x", "envx", Arg::Numeric,		"     -x <int> --envx <int>  \tWidth of the environment. Default: 10."},
	{ENVY, 0, "y", "envy", Arg::Numeric,		"     -y <int> --envy <int>  \tHeight of the environment. Default: 10."},
	{ENVZ, 0, "z", "envz", Arg::Numeric,		"     -z <int> --envz <int>  \tDepth of the environment. Default: 10."},
	{NBITER, 0, "i", "nb_iteration", Arg::Numeric,	"     -i <int> --nb_iteration <int>  \tThe number of iterations = the number of time a plan is choose (half of monomers). Default: 1000."},
	{WRITEACH, 0, "w", "write_every", Arg::Numeric,"     -w <int> --write_every <int>  \tWrite the polymer state every X iterations. Default: 1."},
	{UNKNOWN, 0, "", "", Arg::Unknown, 		"\nParallelism options:"},
	{NB_THREAD, 0, "t", "nb_thread", Arg::Numeric,	"     -t <int> --nb_thread <int>  \tThe number of threads to use to move monomers. So there will be N+3 threads running as there is always 3 running (not for moves): the main process + 2 writers. Default: 1"},
	{PACKAGE_SIZE, 0, "p", "jobs_pack_size", Arg::Numeric,		"     -p <int> --jobs_pack_size <int>  \tThe size of a package. For each iteration, monomers to move are regroup to a package of a certain size, this package is put to the concurrent list for threads which work on this package. Impact the simulation speed. Default: 1."},
	{GPU, 0, "g", "gpu_use", Arg::None, "     -g --gpu_use  \tSwitch to opencl (gpu) mover instead of the cpu mover."},
	{UNKNOWN, 0, "", "", Arg::Unknown,			"\nExamples:\n"
							"     ./mst -f path/to/polyfile -x 100 -y 100 -z 100\n"
							"     ./mst -f path/to/polyfile -o output/directory/ -x 100 -y 100 -z 100 -i 10000 -w 100 -t 4 -p 10\n"
							"     ./mst -f path/to/polyfile -o output/directory/ -x 100 -y 100 -z 100 -i 10000 -w 100 -g\n"
	},
	*/
	{0, 0, 0 ,0, 0, 0} //Necessary
};

Params parse(int argc, char* argv[]){
	Params p;	
	
	//Parse
	argc-=(argc>0); argv+=(argc>0); // Skip program name argv[0] if present
	option::Stats stats(usage, argc, argv);
	vector<option::Option> options(stats.options_max);
	vector<option::Option> buffer(stats.buffer_max);
	option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

	//Stopping program checks
	if (parse.error()) exit(1);
	if (options[HELP] || argc == 0){
		option::printUsage(cout, usage);
		exit(0);
	}
	if(options[VERSION]){
		cout << PROG_VERSION << endl;
		exit(0);
	}

	//Set the params
		//Input read files
	if(!options[READSFILE]){ //Test the file is set
		cerr << "ERROR: reads file must be set (-f | --file)." << endl;
		exit(1);
	}
	else p.readsFile=options[READSFILE].arg;
		//Output information destination
	if(options[LOGS]){p.logs = options[LOGS].arg;}
		//Required score
	if(options[REQUIREDSCORE]){
		int8_t n = stoi(options[REQUIREDSCORE].arg);
		//Expected the number to be > 0 and <= 100
		if(n <= 0 or n > 100){
			cerr << "ERROR: expected a required score number ]0;100]. Value asked: " << to_string(n) << endl;
			exit(1);
		}
		p.requiredScore = n;
	}
		//Number of threads
	if(options[NBTHREADS]){
		unsigned n = stoul(options[NBTHREADS].arg);
		//Expected the number to be > 0
		if(n < 1){
			cerr << "ERROR: expected a number of threads >= 1. Value asked: " << to_string(n) << endl;
			exit(1);
		}
		else{
			//check the number of threads <= maximal number of threads supported by the hardware
			unsigned max = thread::hardware_concurrency(); // 0 => Can't know it
			if(max == 0){
				cerr << "WARNING: Impossible to know the hardware limit about the number of threads." << endl;
			}
			else if (max < n){
				cerr << "ERROR: Number of threads set (" << to_string(n)
					<< ") is superior of the number of threads supported ("
					<< to_string(max) << ")." << endl;
				exit(1);
			}
		}
		p.nbthreads = n;
	}

		//Get OpenCL params
	if(options[GPUINFOS]){p.gpu_infos = true;}
	if(options[GPU]){
		p.gpu = true;
		if(options[CLPLATFORMID]){p.opencl_platform_id = stoul(options[CLPLATFORMID].arg);}
		if(options[CLDEVICEID]){p.opencl_device_id = stoul(options[CLDEVICEID].arg);}
	}
	/*
	if(options[OUTPUT]){p.outPrefix = options[OUTPUT].arg;}
	if(options[ENVX]){p.size_x = stoull(options[ENVX].arg);}
	if(options[ENVY]){p.size_y = stoull(options[ENVY].arg);}
	if(options[ENVZ]){p.size_z = stoull(options[ENVZ].arg);}
	if(options[NBITER]){p.nb_iteration = stoull(options[NBITER].arg);}
	if(options[WRITEACH]){p.write_every = stoull(options[WRITEACH].arg);}
	if(options[NB_THREAD]){p.parallel_nb_threads = stoi(options[NB_THREAD].arg);}
	if(options[PACKAGE_SIZE]){p.jobs_groupsize = stoull(options[PACKAGE_SIZE].arg);}
	*/

	return p;
}



#endif
