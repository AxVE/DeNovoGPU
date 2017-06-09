R"CLCODE(
	#define BASE_A 65
	#define BASE_T 84
	#define BASE_G 71
	#define BASE_C 67
	#define BASE_N 78
	
	//This function return the complet of a given char
	char complement(char c){
		switch(c){
			case BASE_A:
				c = BASE_T;
				break;
			case BASE_T:
				c = BASE_A;
				break;
			case BASE_G:
				c = BASE_C;
				break;
			case BASE_C:
				c = BASE_G;
				break;
			default:
				c = BASE_N;
		}
		return c;
	}

	//This function return the match score of seq2 on seq1. The array buffer intarray must be (at least) of the size of seq1 (so seq1_size).
	char score_2_seq(local char *seq1, unsigned long seq1_size, local char *seq2, unsigned long seq2_size, global long *intarray){
		/*
		 * Using the need1a algorithm (needleman but with only a 1D int array of seq1 size instead of a 2D array (seq1*seq2 sizes))
		*/
		long best=0; //Use to find the best score
		long previous_align_score = 0;//Use to keep in memory the previous "(mis)match" score
		//Doing first char of seq2 (can start everywhere in seq1)
		for(size_t i=0; i < seq1_size; i++){
			if(seq1[i]==seq2[0]){intarray[i]=1;}
			else{intarray[i]=-1;}
		}
		best = intarray[seq1_size-1];
		//Doing the rest of seq2
		for(size_t j=1; j < seq2_size; j++){
			previous_align_score = intarray[0];

			//Against first char of seq1: is an indel so it's -1*number_of_missing_nuc
			intarray[0] = -j;
			
			//Do others nucs
			for(size_t i=1; i < seq1_size; i++){
				//best score with (mis)match or indel ?
				long s = previous_align_score + ((seq1[i]==seq2[j])?1:-1);
				if(s < intarray[i-1]){s = intarray[i-1]-1;}
				if(s < intarray[i]){s = intarray[i]-1;}
				//Get 'current' align score
				previous_align_score = intarray[i];
				//Put new score
				intarray[i] = s;
			}

			//Is it the best align ?
			if(intarray[seq1_size-1]>best){best=intarray[seq1_size-1];}
		}

		//Check best case is not when seq2 is inside seq1 (intarray is now with the last nuc of seq2)
		for(size_t i=0; i < seq1_size; i ++){
			if(intarray[i] > best){best = intarray[i];}
		}
	
		//Send best score in %
		if(best < 0){best=0;}
		unsigned long min_size = (seq1_size < seq2_size)?seq1_size:seq2_size;
		return 100*best/min_size;
	}

	kernel void cmp_2_contigs(__global unsigned long *infos, __global char *scores, __global unsigned long *seqs_sizes, __global char *ultraseq, __global long *intbufloc, __local char *charbufloc){
		size_t gid = get_global_id(0);
		size_t seq2_id = gid/infos[0];
		size_t seq1_id = gid - seq2_id*infos[0];
		size_t work_id = get_local_id(0);
		unsigned long nbContigs = infos[0];

		//Only calculate score if it's not the same sequence
		if(seq1_id != seq2_id){
			//Get the first contig sequence and its infos. As it will be read multiple times, copy it in local-item buffer
				//Get size
			unsigned long seq1_size = seqs_sizes[seq1_id];
				//Prepare the buffer to use
			local char* seq1 = &charbufloc[infos[2]*work_id*2];
				//Get the position of the first sequence inside ultraseq
			unsigned long start = 0;
			for(unsigned long i=0; i < seq1_id; i++){start += seqs_sizes[i];}
				//Copy the contig sequence in local buffer (copy each char one by one from global to local)
			for(unsigned long c=0; c < seq1_size; c++){
				seq1[c] = ultraseq[start+c];
			}

			//Get the second contig sequence and its infos. As it will be reversed-complement, copy it in local-item-2 buffer
				//Get size
			unsigned long seq2_size = seqs_sizes[seq2_id];
				//Prepare the buffer to use
			local char* seq2 = &charbufloc[infos[2]*work_id*2+infos[2]];
				//Get the position of the first sequence inside ultraseq
			start = 0;
			for(unsigned long i=0; i < seq2_id; i++){start += seqs_sizes[i];}
				//Copy the contig sequence in local buffer (copy each char one by one from global to local)
			for(unsigned long c=0; c < seq2_size; c++){
				seq2[c] = ultraseq[start+c];
			}
			
			//Get the global int array buffer (it's an array of nbElem*longuestContig of long),
			//so the sub_array to this work item is at the position longuestContig*global_id
			global long *inta = &intbufloc[infos[2]*gid];

			//Get match score of seq2 on seq1
			//scores[seq1_id+nbContigs*seq2_id]=score_2_seq(seq1, seq1_size, seq2, seq2_size, inta);
				//normal
			char s_normal = score_2_seq(seq1, seq1_size, seq2, seq2_size, inta);
				//reverse-complement
			for(size_t i=0; i < seq2_size/2; i++){
				char c = complement(seq2[i]);
				seq2[i] = complement(seq2[seq2_size-1-i]);
				seq2[seq2_size-1-i] = c;
			}
				//Don't forget the middle nuc if seq2_size is odd
			if(seq2_size%2){seq2[seq2_size/2]=complement(seq2[seq2_size/2]);}
				//get score with reverse complement
			char s_reverse = score_2_seq(seq1, seq1_size, seq2, seq2_size, inta);

				//Get best score
			if(s_reverse > s_normal){scores[seq1_id+nbContigs*seq2_id] = s_reverse * -1;}
			else{scores[seq1_id+nbContigs*seq2_id] = s_normal;}

		}
		else{ //same sequences
			scores[seq1_id+nbContigs*seq2_id]=0;
		}
	}
)CLCODE";
