#include "readsTools.hpp"

#include <string>
#include <vector>
#include <iostream>

using namespace std;
	
//Static members/attributes;
const std::map<char,uint8_t> ReadsTools::nuc = {{'A',0},{'T',1},{'G',2},{'C',3},{'X',4}};

const std::vector< std::vector< int8_t > > ReadsTools::scoreMatrix =
	std::vector< std::vector< int8_t > > (5,
		std::vector<int8_t>(5,0)) = {
			//A T G C X(2nd id)
			{1,-1,-1,-1,-1}, //A (first id)
			{-1,1,-1,-1,-1}, //T (first id)
			{-1,-1,1,-1,-1}, //G (first id)
			{-1,-1,-1,1,-1}, //C (first id)
			{-1,-1,-1,-1,-1} //X (first id)
		};

	//function to get the best value in a matrix of int8_t (the scoreMatrix)
inline const int8_t bestValueInMatrix(vector< vector< int8_t> > matrix){
	if(!matrix.size() || !matrix[0].size()){
		cerr << "Error: matrix empty." << endl;
		return 0;
	}
	
	int8_t b = matrix[0][0];
	for(vector<vector<int8_t>>::const_iterator vi = matrix.cbegin(); vi != matrix.cend(); vi ++){
		for(vector<int8_t>::const_iterator vvi = vi->cbegin(); vvi != vi->cend(); vvi++){
			if(*vvi > b){b=*vvi;}
		}
	}
	return b;
}

const int8_t ReadsTools::bestScoreMatrix = bestValueInMatrix(scoreMatrix);
const int ReadsTools::gap_score = -1;

//Tool Public function
const int ReadsTools::scoreCmpReads(const std::string& read1, const string& read2){
	/* normal mapping */
	int best = scoreMatchingReads(read1,read2);

	/* read1 normal, read2 complementary (r2r) */
	string r2r = read2;
	const size_t r2size = read2.size();
	for(size_t i=0; i < r2size; i++){
		size_t r2id = r2size-1-i;
		switch(read2[i]){
			case 'A':
				r2r[r2id]='T';
				break;
			case 'T':
				r2r[r2id]='A';
				break;
			case 'G':
				r2r[r2id]='C';
				break;
			case 'C':
				r2r[r2id]='G';
				break;
			default:
				//cerr << "Unknown nuc: " << read2[i] << endl;
				r2r[r2id]='A';
				break;
		}
	}
	int score = scoreMatchingReads(read1, r2r);

	if(best < score){best=score*(-1);}

	return best;
}

const string ReadsTools::allignSeq(const string& seq1, const string& seq2){
	//Reminder: a needlemap has a size of seq1.size()+1 and seq1.size()+1
	//because the first line and first column must include an "empty" possibility
	struct Node {
		Node* previous = nullptr;
		size_t s1 = 0; //seq1 pos
		size_t s2 = 0; //seq2 pos
		int score=0;
	};
	vector< vector<Node> > needlemap (seq1.size()+1, vector<Node>(seq2.size()+1, Node()));

	/*
	cerr << "Merging following sequences" << endl;
	cerr << "Seq1:\n" << seq1 << endl;
	cerr << "Seq2:\n" << seq2 << endl;
	*/

	//Needlemap completion
		//first nuc (0,0)
		needlemap[0][0].score=scoreMatrix[nuc.at(seq1[0])][nuc.at(seq2[0])];
		//first line (we can begin anywere on the seq1 (reference) so
		//initial score is 0)
	for(size_t i=1; i < seq1.size()+1; i++){
		needlemap[i][0].score = 0;
		needlemap[i][0].previous = &needlemap[i-1][0];
		needlemap[i][0].s1=i-1;
	}
		//first column (we must match the beginning of the seq2)
	for(size_t i=1; i < seq2.size()+1; i++){
		needlemap[0][i].score = needlemap[0][i-1].score + gap_score;
		needlemap[0][i].previous = &needlemap[0][i-1];
		needlemap[0][i].s2=i-1;
	}

		//Do all other nodes (except last line and column)
	for(size_t i=1; i < seq1.size(); i++){
		for(size_t j=1; j < seq2.size(); j++){
				//Set position
			needlemap[i][j].s1=i-1;
			needlemap[i][j].s2=j-1;
				//(mis)match ?
			needlemap[i][j].score = needlemap[i-1][j-1].score +
				scoreMatrix[nuc.at(seq1[i-1])][nuc.at(seq2[j-1])];
			needlemap[i][j].previous = &needlemap[i-1][j-1];
				//seq1 gap ?
			int s = needlemap[i][j-1].score + gap_score;
			if(s > needlemap[i][j].score){
				needlemap[i][j].score = s;
				needlemap[i][j].previous = &needlemap[i][j-1];
			}
				//seq2 gap ?
			s = needlemap[i-1][j].score + gap_score;
			if(s > needlemap[i][j].score){
				needlemap[i][j].score = s;
				needlemap[i][j].previous = &needlemap[i-1][j];
			}
		}
	}

	//Last line and last column
	size_t n2 = seq2.size();
	for(size_t n1=1; n1 < seq1.size(); n1++){
		//Set pos
		needlemap[n1][n2].s1=n1-1;
		needlemap[n1][n2].s2=n2-1;

		//Match/mismatch ?
		needlemap[n1][n2].score = needlemap[n1-1][n2-1].score +
			scoreMatrix[nuc.at(seq1[n1-1])][nuc.at(seq2[n2-1])];
		needlemap[n1][n2].previous = &needlemap[n1-1][n2-1];

		//seq1 gap ? 
		int s = needlemap[n1][n2-1].score+gap_score;
		if(s > needlemap[n1][n2].score){
			needlemap[n1][n2].score = s;
			needlemap[n1][n2].previous = &needlemap[n1][n2-1];
		}
		//seq2 gap ? (There is no gap_score as we can accept "gap" at the end )
		s = needlemap[n1-1][n2].score;
		if(s > needlemap[n1][n2].score){
			needlemap[n1][n2].score = s;
			needlemap[n1][n2].previous = &needlemap[n1-1][n2];
		}
	}
	size_t n1 = seq1.size();
	for(size_t n2=1; n2 < seq2.size(); n2++){
		//Set pos
		needlemap[n1][n2].s1=n1-1;
		needlemap[n1][n2].s2=n2-1;

		//Match/mismatch ?
		needlemap[n1][n2].score = needlemap[n1-1][n2-1].score +
			scoreMatrix[nuc.at(seq1[n1-1])][nuc.at(seq2[n2-1])];
		needlemap[n1][n2].previous = &needlemap[n1-1][n2-1];

		//seq1 gap ? (There is no gap_score as we can accept "gap" at the end)
		int s = needlemap[n1][n2-1].score;
		if(s > needlemap[n1][n2].score){
			needlemap[n1][n2].score = s;
			needlemap[n1][n2].previous = &needlemap[n1][n2-1];
		}
		//seq2 gap ?
		s = needlemap[n1-1][n2].score+gap_score;
		if(s > needlemap[n1][n2].score){
			needlemap[n1][n2].score = s;
			needlemap[n1][n2].previous = &needlemap[n1-1][n2];
		}
	}

	//Do last row
		//Set pos
	needlemap[n1][n2].s1=n1-1;
	needlemap[n1][n2].s2=n2-1;
		
		//match/mismatch
	needlemap[n1][n2].score = needlemap[n1-1][n2-1].score +
		scoreMatrix[nuc.at(seq1[n1-1])][nuc.at(seq2[n2-1])];
	needlemap[n1][n2].previous = &needlemap[n1-1][n2-1];
	
	//Gap (no gap_score as it is the end of reads)
	if(needlemap[n1][n2].score < needlemap[n1-1][n2].score){
		needlemap[n1][n2].score = needlemap[n1-1][n2].score;
		needlemap[n1][n2].previous = &needlemap[n1-1][n2];
	}
	if(needlemap[n1][n2].score < needlemap[n1][n2-1].score){
		needlemap[n1][n2].score = needlemap[n1][n2-1].score;
		needlemap[n1][n2].previous = &needlemap[n1][n2-1];
	}
		
	//Needlemap lecture (begin by the end : sequence is reversed)
	Node* n = &needlemap[n1][n2];
	if(n == nullptr || n->previous == nullptr){
		cerr << "Error in needlemap reading" << endl;
		return "error while creating the allignment";
	}
	Node* p = n->previous;
	string s = "";
		
		//Reconstruct all nucleotides
	while(p != nullptr){
		//Do the last nucleotide
		//If previous is s1-1 and s2-1 so it's a (mis)match
		//Else it is a gap
		if(p->s1 < n->s1 && p->s2 < n->s2){
			if(seq1[n->s1] == seq2[n->s2]){
				s += seq1[n->s1];
			}
			else { s += 'X';}
		}
			//seq1 gap
		else if (p->s1 == n->s1){
			s += seq2[n->s2];
		}
			//seq2 gap
		else {
			s += seq1[n->s1];
		}

		//Do previous
		n = p;
		p = n->previous;
	}

	//Return the consensus seq (reverse the needlemap result)
	return string(s.rbegin(), s.rend());

}

//Tool Private function
const int ReadsTools::scoreMatchingReads(const string& read1, const string& read2){
	/*
	Calculing Comparaison score between the 2 read2s using
	needleman&wunsch algorithm without memoring the best path
	as we only need the score
	*/

	//Initialize mapping array
	vector< vector<int> > needlemap = vector< vector<int> >(read1.size(),
		vector<int>(read2.size(),0));
	
	/*
	Fill the array (the best score is the best value inside
	the last column and the last line).
	*/
	
	//Do the first value (row[0][0])
	needlemap[0][0]=scoreMatrix[nuc.at(read1[0])][nuc.at(read2[0])];
	//Do all the first line and the first column
	//It's all a read2 again the first nuc. Vice-versa.
		//read2 against first nuc of read1
	for(size_t n=1; n < read2.size(); n++){
		//Gap or (mis)match, keep the best
		int s=scoreMatrix[nuc.at(read1[0])][nuc.at(read2[n])];
		if(needlemap[0][n-1] + gap_score > s){s=needlemap[0][n-1]+gap_score;}
		needlemap[0][n]=s;
	}
		//read1 against first nuc of read2
	for(size_t n=1; n < read1.size(); n++){
		//Gap or (mis)match, keep the best
		int s=scoreMatrix[nuc.at(read1[n])][nuc.at(read2[0])];
		if(needlemap[n-1][0] + gap_score > s){s=needlemap[n-1][0]+gap_score;}
		needlemap[n][0]=s;
	}

	//Complete the array : do all other cases (except last column and last line)
	for(size_t n1=1; n1 < read1.size()-1; n1++){
		for(size_t n2=1; n2 < read2.size()-1; n2++){
			//Is it best to include a match or have a (mis)match ?
			//Get the mismatch score
			int s = needlemap[n1-1][n2-1]
				+ scoreMatrix[nuc.at(read1[n1])][nuc.at(read2[n2])];
			//Include a gap (2 possibilities) is better ?
			if(needlemap[n1][n2-1]+gap_score > s){s=needlemap[n1][n2-1]+gap_score;}
			if(needlemap[n1-1][n2]+gap_score > s){s=needlemap[n1-1][n2]+gap_score;}
			//Keep the best score in the current row
			needlemap[n1][n2]=s;
		}
	}

	//Do last line (permissiv to not finish read2) and search best score
	int best=0;
	const size_t n1=read1.size()-1;
	for(size_t n=1; n < read2.size(); n++){
		//Gap or mismatch ?
		int s=scoreMatrix[nuc.at(read1[n1])][nuc.at(read2[n])]+needlemap[n1-1][n-1];
		if(needlemap[n1-1][n]+gap_score > s){s=needlemap[n1][n]+gap_score;}
		if(needlemap[n1][n-1] > s){s=needlemap[n1][n-1];}
		needlemap[n1][n]=s;
		
		//Best ?
		if(s>best){best=s;}
	}

	//Do last column (permissiv to not finish read1) and search best score
	const size_t n2 = read2.size()-1;
	for(size_t n=1; n < read1.size(); n++){
		//Gap or mismatch ?
		int s=scoreMatrix[nuc.at(read1[n])][nuc.at(read2[n2])]+needlemap[n-1][n2-1];
		if(needlemap[n][n2-1]+gap_score > s){s=needlemap[n][n2-1]+gap_score;}
		if(needlemap[n-1][n2] > s){s=needlemap[n-1][n2];}
		needlemap[n][n2]=s;
		
		//Best ?
		if(s>best){best=s;}
	}

	//Do last row
	int s=scoreMatrix[nuc.at(read1[read1.size()-1])][nuc.at(read2[read2.size()-1])];
		//best gap ?
	if(needlemap[read1.size()-1][read2.size()-2] > s){needlemap[read1.size()-1][read2.size()-2];}
	if(needlemap[read1.size()-2][read2.size()-1] > s){needlemap[read1.size()-2][read2.size()-1];}
		//best score ?
	if(s>best){best=s;}

	//Calculate the maximum possible best score
	int maxbest = min(read1.size(),read2.size())*bestScoreMatrix;
	//output the best
	return 100*best/maxbest;
}
