#include "reads.hpp"
#include "readsTools.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <utility> //pair
#include <queue>
#include <cctype> //toupper()
#include <random>

using namespace std;

/* SEQ */
//Seq constructor
Seq::Seq(const std::string& name, const std::string& nucleotides_seq){
	//Store basic info
	m_name=name;
	m_size=nucleotides_seq.size();
	
	//Prepare the array of encoded nucleotides
	seqLength nbC = (m_size)/4 + (m_size%4 ? 1 : 0);
	m_encNuc=new unsigned char[nbC];

	/* Encode :
		Left to rights nucleotides
		Input on right and shift to left
	*/
	//Encode the completes 4based elements
	for(seqLength i=0; i < nbC-(m_size%4 ? 1 : 0); i++){ 
		//Set to null the element
		//*(m_encNuc+i) = 0b00000000;
		unsigned char e = 0b00000000;
		
		//Encode nucleotides
		for(unsigned char c = 0; c < 4; c++){
			//Shift the elem
			e <<= 2;
			//Get the nuc
			e |= encode_nuc.at(nucleotides_seq[i*4+c]);
		}
		*(m_encNuc+i) = e;
	}
	//Encode the last ones if they are [0;3]
	if(m_size%4){
		seqLength i = nbC-1;
		//Set to null the element
		*(m_encNuc+i) = 0b00000000;
	
		//Encode nucleotides
		for(unsigned char c=0; c < m_size%4; c++){
			//Shift the elem
			*(m_encNuc+i) <<= 2;
			//Get the nuc
			*(m_encNuc+i) |= encode_nuc.at(nucleotides_seq[i*4+c]);
		}
		//Fill "empty" spaces
		*(m_encNuc+i) <<= 2*(4-m_size%4);
	}
}

//Copy constructor
Seq::Seq(const Seq& other){
	//get infos
	m_name = other.get_name();
	m_size = other.get_size();

	//Get encoded sequence
	seqLength nbC = m_size/4 + (m_size%4 ? 1 : 0);
	m_encNuc = new unsigned char[nbC];
	unsigned char * otherEnc = other.get_encNuc();
	for(seqLength i=0; i < nbC; i++){
		*(m_encNuc + i) = *(otherEnc + i);
	}
}

//Move constructor
Seq::Seq(Seq& other) noexcept{
	m_size = other.get_size();
	m_name = other.get_name();
	m_encNuc = other.get_encNuc();
	other.reset();
}

//Seq destructor
Seq::~Seq(){
	delete[] m_encNuc;
}

//Copy assignment operator
Seq& Seq::operator= (const Seq& other){
	Seq tmp(other); //Re-use copy-constructor
	*this =std::move(tmp); //Re-use move-assignment
	return *this;
}

//Move assignment operator
Seq& Seq::operator= (Seq&& other) noexcept{
	m_size=other.get_size();
	m_name=other.get_name();
	delete[] m_encNuc;
	m_encNuc=other.get_encNuc();
	other.reset();
	return *this;
}

//Accessor of the decoded sequence
void Seq::get_seq(string& out_seq) const{
	//Prepapre the string
	out_seq.clear();
	out_seq.resize(m_size);
	
	//Get all 4 bases until the last char (a char = 4 bases except the last char)
	for(seqLength i=0; i < (m_size)/4; i++){
		unsigned char buf = *(m_encNuc+i);
		unsigned char mask = 0b00000011;

		unsigned char n=4;
		while(n>0){
			n--;
			char c = decode_nuc.at(buf&mask);
			buf >>= 2;
			out_seq[i*4+n] = c;
		}

	}

	//Get the nuc in the last 4 ([1;4]) Probleme here?
	if(m_size%4){
		seqLength i = m_size/4;
		unsigned char buf = *(m_encNuc+i);
		unsigned char mask = 0b00000011;

		//Remove the "unused" nuc
		buf >>= 2*(4-m_size%4);

		//Get the last 3||2||1 nuc

		unsigned char n=m_size%4;
		while(n > 0){
			n--;
			char c = decode_nuc.at(buf&mask);
			buf >>= 2;
			out_seq[m_size-m_size%4+n] = c;
		}
	}
}


/* READS */
int Reads::readFile_fasta(const std::string& file){
	//Prepare list of reads
	queue< pair< string, string > > toStoreReads; //Queue of reads : name + sequence
	unsigned int countReads=0;
	
	//Prepare random (for N, W, B, ...)
	default_random_engine generator(random_device{}()); //Random generator with a stochastic random seed
	uniform_int_distribution<unsigned short int> distrib2 (0,1);
	uniform_int_distribution<unsigned short int> distrib3 (0,2);
	uniform_int_distribution<unsigned short int> distrib4 (0,3);


	//Init file reader
	ifstream fileIn (file.c_str());

	//Check the file
	if(!fileIn.good()){
		m_cerr << "Error opening '" << file << "'. Bailing out." << endl;
		return -1;
	}
	
	//File ok, get content (reads and their name)
	else{
		string line, name, content;
		while(getline(fileIn, line).good()){
			//new read ? (name of read begins with a '>'
			if(line.empty() || line[0] == '>'){
				if(!name.empty()){
					//Increment the counter of reads
					++countReads;
					//Add this reads to the queue
					toStoreReads.push(make_pair(name, content));
					name.clear();
				}
				if(!line.empty()){
					//Remove the '>'
					name = line.substr(1);
				}
				content.clear();
			}
			// Or get the read content
			else if (!name.empty()) {
				//No space allowed in the middle of the seq
				if(line.find(' ') != string::npos){
					name.clear();
					content.clear();
				}
				else {
					//Clean and correct the line
					size_t i=0;
					while(i < line.size()){
						//min to maj
						//line[i]=toupper(i);
						char c = line[i];
						if(c=='A' || c=='T' || c=='G' || c=='C'){
							i++;
						}
						//Remove \r and \n 
						else if(line[i]=='\r' || line[i]=='\n'){
							line.erase(i, 1);
						}
						//Convert other letters
						else {
							string p = "ATGC";
							switch(c){
								case 'N': //A, T, G or C
									p = "ATGC";
									line[i] = p[distrib4(generator)];
									break;
								case 'R': //A or G
									p = "AG";
									line[i] = p[distrib2(generator)];
									break;  
								case 'Y': //T or C
									p = "TC";
									line[i] = p[distrib2(generator)];
									break;  
								case 'K': //T or G
									p = "TG";
									line[i] = p[distrib2(generator)];
									break;  
								case 'M': //A or C
									p = "AC";
									line[i] = p[distrib2(generator)];
									break;  
								case 'S': //G or C
									p = "GC";
									line[i] = p[distrib2(generator)];
									break;  
								case 'W': //A or T
									p = "AT";
									line[i] = p[distrib2(generator)];
									break;  
								case 'B': //T, G or C
									p = "TGC";
									line[i] = p[distrib3(generator)];
									break;  
								case 'D': //A, T or G
									p = "ATG";
									line[i] = p[distrib3(generator)];
									break;  
								case 'H': //A, T or C
									p = "ATC";
									line[i] = p[distrib3(generator)];
									break;  
								case 'V': //A, G or C
									p = "AGC";
									line[i] = p[distrib3(generator)];
									break; 
								default:
									cerr << "Unknown letter: " << c << ". Random between A, T, G and C." << endl;
									p = "ATGC";
									line[i] = p[distrib4(generator)];
							}
							i++;
						}
					}
					//Add the line to the content
					content += line;
				}
			}
		}
	
		//If there is a last reads
		if(!name.empty() && !content.empty()){
			//Increment the counter of reads
			++countReads;
			//Add this reads to the queue
			toStoreReads.push(make_pair(name, content));
		}
	}

	//Encode queue
	unsigned int i =0;
	while(!toStoreReads.empty()){
		//Display progress
		i++;
		//m_cout << "\rencoding: " << 100 - (toStoreReads.size()*100/countReads) << " % (" <<  i << "/" << countReads << ")" << flush;
		//Encode and store the current read
		m_seqs.emplace_back(toStoreReads.front().first, toStoreReads.front().second);
		//Remove it from the queue
		toStoreReads.pop();
	}
	//m_cout << "\rencoding: 100 % (" << i << "/" << countReads << ")" << endl;


	return 0;
}

bool Reads::display_read(const size_t id, ostream& os) const{
	if(id >= m_seqs.size()){
		cerr << "ERROR: sequence id asked > number of sequences" << endl;
		return 1;
	}
	os << ">" << m_seqs[id].get_name() << endl;
	string seq;
	m_seqs[id].get_seq(seq);
	os << seq << endl;
	return 0;
}

bool Reads::get_readName(const size_t id, string& str) const {
	if(id >= m_seqs.size()){
		cerr << "ERROR: sequence id asked > number of sequences" << endl;
		str = "No seq";
		return 1;
	}
	str = m_seqs[id].get_name();
	return 0;
}

bool Reads::get_readSeq(const size_t id, string& str) const {
	if(id >= m_seqs.size()){
		cerr << "ERROR: sequence id asked > number of sequences" << endl;
		str = "No seq";
		return 1;
	}
	m_seqs[id].get_seq(str);
	return 0;
}

/* Contigs */
Contigs::Contigs(const Reads* reads_set){
	//Keep memory where is this reads set
	m_reads_set = reads_set;

	//The max number of contig (the number of contigs will reduce after) is the number of reads so initialize the Contigs* arrays at this size
	size_t nbReads = reads_set->get_nbReads();
	m_real_contigs_list.resize(nbReads, nullptr);

	//For each read, create a contig, get its sequence, link correctly the "read" to the contig.
	for(size_t i = 0; i < nbReads; i++){
		Contig* c = new Contig; //Create contig
		c -> c_reads_ids.insert(i); //Store the read id
		reads_set->get_readSeq(i, c -> c_seq); //The sequence of the contig = read's sequence
		//Keep this contigs link
		m_real_contigs_list[i] = c;
	}

	//Sync the usable contigs list to be up-to-date
	sync_contigs_list();
}

Contigs::~Contigs(){
	//Desalloc all contigs
	for(size_t i=0; i < m_contigs_list.size(); i++){
		delete m_real_contigs_list[i];
		m_real_contigs_list[i] = nullptr;
	}
	sync_contigs_list();
}

void Contigs::merge_contigs(size_t contigID1, size_t contigID2, const int score){
	//Check values
	size_t nbC = m_contigs_list.size();
	if(contigID1 >= nbC){
		cerr << "Error: first id (" << contigID1 <<
			") is > number of known contigs (" << nbC << ")"
		<< endl;
	}
	if(contigID2 >= nbC){
		cerr << "Error: second id (" << contigID2 <<
			") is > number of known contigs (" << nbC << ")"
		<< endl;
	}
	if(contigID1 == contigID2){
		cerr << "Error: first and second id have the same value ("
			<< contigID1 << ")"
		<< endl;
	}

	//Get the addresses
	Contig* contig1 = m_contigs_list[contigID1];
	Contig* contig2 = m_contigs_list[contigID2];

	//Alligned the 2 sequences and keep this result in the contig1 sequence
	if(score >= 0){
	}
	else {
		//It is the reverse complement
		string s = ReadsTools::allignSeq(contig1->c_seq, contig2->c_seq);
		string rs = s; //reverse-complement
		size_t l = s.size()-1;
		for(size_t i=0; i<s.size(); i++){
			switch (s[i]){
				case 'A':
					rs[l-i] = 'T';
					break;
				case 'T':
					rs[l-i] = 'A';
					break;
				case 'G':
					rs[l-i] = 'C';
					break;
				case 'C':
					rs[l-i] = 'G';
					break;
				default:
					rs[l-i] = s[i];
					break;
			}
		}
		contig1->c_seq = ReadsTools::allignSeq(contig1->c_seq, rs);
	}
	
	//Add contig2 reads id in contig1 ones
	contig1->c_reads_ids.insert(contig2->c_reads_ids.begin(), contig2->c_reads_ids.end());
	
	//Search all "to contig2 link" and replace them by contig 1
	for(size_t i=0; i < nbC; i++){
		if(m_contigs_list[i] == contig2){m_contigs_list[i] = contig1;}
	}

	//Free and remove contig2
	for(size_t i=0; i < m_real_contigs_list.size(); i++){
		if(m_real_contigs_list[i] == contig2){
			//Free the contig
			delete contig2;
			m_real_contigs_list[i] = nullptr;
			//Erease element from the real list and stop the search
			m_real_contigs_list.erase(m_real_contigs_list.begin()+i);
			i = m_real_contigs_list.size();
		}
	}
}

void Contigs::sync_contigs_list(){m_contigs_list = m_real_contigs_list;}

