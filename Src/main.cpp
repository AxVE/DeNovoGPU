#define PROG_VERSION "v0.1.0-alpha"

//Owns files
#include "cmdparser.hpp"
#include "log.hpp"
#include "reads.hpp"
#include "readsTools.hpp"
#include "workerCL.hpp"

//Standards libs
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <thread>


using namespace std;

/*
Function to get the scores matrix between all contigs of a same set. Used in mono and multithreading.
*/

void contigs_cmp(size_t id_begin, size_t id_end, vector< vector< int8_t> >& scores, const Contigs& contigs){
	ReadsTools readsTools; //Toolsbox containing the function to cmp reads
	for(size_t c = 0; c < contigs.get_nbContigs(); c++){
		string seq_c = contigs.get_seqContig(c);
		for(size_t i=id_begin; i <= id_end; i++){
			//if same read, score = 0
			//Else compare
			if(c == i){scores[c][i] = 0;}
			else{
				string seq_i = contigs.get_seqContig(i);
				
				scores[c][i] = readsTools.scoreCmpReads(seq_c, seq_i);
			}
		}
	}
}


/*
Heart of the program
*/

int main(int argc, char* argv[]){
	/* Initialisation :
		1/ Check and parse parameters and options
		2/ Display infos
		3/ Initialize OpenCL
	*/
	//Get parameters and options (Includes/cmdparser.hpp using optionsparser.h)
	Params params = parse(argc, argv);
	Log log(params.logs);

	//Display infos
	log.write( "[STEP]\tInitialisation" );

	if(!params.gpu_infos){ //normal infos
		log.write( "=== Infos ===" );
		log.write(params.infosToString());
	}
	else{ //Asked GPU infos
		log.write( "=== OpenCL infos ===" );
		WorkerCL::list_infos(log);
		return 0; //End prog
	}

	//Initialize gpu (used only if gpu is asked)
	WorkerCL* workerCL = nullptr;
	try{
		if(params.gpu){
			log.write("=== OpenCL initialisation ===");
			workerCL = new WorkerCL(params.opencl_platform_id, params.opencl_device_id, log);
		}
	}
	catch(string opencl_err){
		log.write(opencl_err);
	}

	/* Data gathering :
		1/ read file (get reads) (encode and store them)
		2/ Create the contigs list. It creates a "contig" for each read. A contig contain the ids of its reads and its sequence
		3/ Initialize the 2D array (vector) to store each reads-reads score alignment
	*/
	log.write( "\n[STEP]\tData gathering" );
	
	log.write( "=== File reading ===" );
	Reads reads;
	reads.readFile_fasta(params.readsFile);
	{
		string s = "nbReads = " + to_string(reads.get_nbReads());
		log.write(s);
	}
	log.write( "=== Create contigs ===" );
	Contigs contigs( &reads );
	size_t nbContigs = contigs.get_nbContigs();
	{
		string s = "number of contigs = " + to_string(nbContigs);
		log.write(s);
	}
	
	
	// Cycles
	log.write( "\n[STEP]\tCycles" );

	bool in_progress = true; //There is still progress ?
	size_t cycles = 0;

	while(in_progress){
		cycles++;
		{
			string s = "===== Cycle "+to_string(cycles)+" start =====";
			log.write(s);
			log.write( "=== Calculate score between each contig ===" );
		}

		/*
		 Create a 2D array of size of nb of contigs storing for each pair of contigs their alignment score.
		 0 means a bad alignment.
		 The second element is match on the first element.
		 >0 gives the score of a normal alignment
		 <0 gives the score (the absolute value) of the complementary alignment
		 The highest absolute value is stored (for example: for two read A et B if the score of the normal is 10 and
		 the score of the complementary alignment is 25, -25 will be stored
		*/
		nbContigs = contigs.get_nbContigs();
		vector< vector< int8_t > > scores = vector< vector<int8_t> >(nbContigs, vector<int8_t>(nbContigs, 0));
		
			// GPU ?
		if(params.gpu){
			scores = workerCL->run(contigs, params.opencl_work_group_size);
		}
			// Multithreading ?
		else if(params.nbthreads > 1){
			vector<thread> workers;
			size_t id_begin=0;

			//the main is a thread; So the number of worker = (number_of_threads_asked - 1).
			//Launch threads
			for(size_t w=0; w < params.nbthreads-1; w++){
				//Number of contigs to cmp in this thread
				size_t nb_contigs_cmp = nbContigs/params.nbthreads + (( w < (nbContigs % params.nbthreads) )?1:0);
				//cout << "Thread " << w << ": " << nb_contigs_cmp << " contigs." << endl;

				//If there is contigs to cmp, we run the thread
				if(nb_contigs_cmp){
					size_t id_end=id_begin+nb_contigs_cmp-1;
					workers.emplace_back(contigs_cmp, id_begin, id_end, ref(scores), ref(contigs));
					//cout << "Thread " << w << " on " << id_begin << " to " << id_end << endl;
					id_begin = id_end+1; //Prepare next
				}

			}
			
			//Run remainings contigs on main
			if(id_begin < nbContigs){
				//cout << "Thread main on " << id_begin << " to " << nbContigs-1 << endl;
				contigs_cmp(id_begin, nbContigs-1, scores, contigs);
			}

			//Wait the end of threads
			for(size_t w = 0; w < workers.size(); w++){
				workers[w].join();
			}

		}
			// Monothreading
		else{
			contigs_cmp(0, nbContigs-1, scores, contigs);
		}

		log.write("---- matrix of scores ----");
		string out;
		for(size_t i=0; i < nbContigs; i++){out += "\t" + to_string(i);}
		log.write(out);

		for(size_t i=0; i < nbContigs; i++){
			out = to_string(i);
			for(size_t j=0; j < nbContigs; j++){
				out += "\t"+to_string(scores[i][j]);
			}
			log.write(out);
		}
		log.write("---- end of matrix ----");

		//Merge best couples (when score >= requirement)
		log.write("=== Merging best contigs couples ===");
		int8_t b; //best couples score
		do {
			b=0;
			size_t b_id1=0;
			size_t b_id2=0;
			// Get best couple
			for(size_t id1=0; id1 < nbContigs; id1++){
				for(size_t id2=id1+1; id2 < nbContigs; id2++){
					//Get score, is it better ?
					if(abs(scores[id1][id2]) > b){
						//keep this couple
						b = abs(scores[id1][id2]);
						b_id1 = id1;
						b_id2 = id2;
					}
				}
			}
			//If score > requirement : merge this 2 
			if(b >= params.requiredScore){
				//Merge
				contigs.merge_contigs(b_id1, b_id2, scores[b_id1][b_id2]);

				//clean the done and not possible couples after this merge
					//This couple and its opposite
				scores[b_id1][b_id2]=0;
				scores[b_id2][b_id1]=0;
					//The 'end' of contig1 is used
						//All couple with contig1 first
				for(size_t i=0; i < nbContigs; i++){
					scores[b_id1][i]=0;
				}
						//All couple with contig1 reversed as second 
				for(size_t i=0; i < nbContigs; i++){
					if(scores[i][b_id1]<0){scores[i][b_id1]=0;}
				}

					//contig2 normal ? The begin is used
				if(scores[b_id1][b_id2] > 0){
					for(size_t i=0; i < nbContigs; i++){
						if(scores[i][b_id2]>0){scores[i][b_id2]=0;}
					}
				}
					//else contig2 is reversed: the end is used
				else{
						//All couples beginning with contig2
					for(size_t i=0; i < nbContigs; i++){
						scores[b_id2][i]=0;
					}
						//All couples with contig2 as second and reversed
					for(size_t i=0; i < nbContigs; i++){
						if(scores[i][b_id2] < 0){scores[i][b_id2]=0;}
					}
				}

			}
		} while (b >= params.requiredScore);


		//Synchronize the contigs list
		log.write("=== Synchronisation ===");
		contigs.sync_contigs_list();
		
		//Get the evolution. None: stop it.
		size_t new_nbContigs = contigs.get_nbContigs();

		string s = "nb contigs merged = " + to_string(nbContigs - new_nbContigs);
		s += "\nnb contigs remaining = " + to_string(new_nbContigs);
		if(new_nbContigs == nbContigs){
			s += "\nNo evolution. Stopping.";
			in_progress=false;
		}
		log.write(s);

		//End of cycle
		log.write("===== Cycle "+to_string(cycles)+" end =====\n");
	}

	//Free workerCL (if created)
	if(params.gpu){
		cout << "Deleting gpu worker" << endl;
		delete workerCL;
		workerCL = nullptr;
		cout << "\t- Done" << endl;
	}

	//End of program
	return 0;
}
