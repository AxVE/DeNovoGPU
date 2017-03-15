#ifndef READS_HPP
#define READS_HPP

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

#define BASE_A 0B00000000
#define BASE_T 0B00000011
#define BASE_G 0B00000010
#define BASE_C 0B00000001

typedef u_int32_t seqLength; //The type to store the size of sequences

const std::map<char, unsigned char> encode_nuc = {
	{'A', BASE_A},
	{'T', BASE_T},
	{'G', BASE_G},
	{'C', BASE_C}
};
const std::map<unsigned char, char> decode_nuc = {
	{BASE_A, 'A'},
	{BASE_T, 'T'},
	{BASE_G, 'G'},
	{BASE_C, 'C'}
};

class Seq {
	public:
		//Constructor
		Seq(const std::string& name, const std::string& nucleotides_seq);
		//Copy constructor
		Seq(const Seq& other);
		//Move constructor
		Seq(Seq& other) noexcept;
		//Destructor
		~Seq() noexcept;
		//Copy assignment operator
		Seq& operator= (const Seq& other);
		//Move assignment operator
		Seq& operator= (Seq&& other) noexcept;

		//Accessors
		seqLength get_size() const {return m_size;}
		std::string get_name() const {return m_name;}
		unsigned char* get_encNuc() const {return m_encNuc;}
		void get_seq(std::string& out_seq) const;

		//Setters
			//Reset put all values to 0 or null without deleting
		void reset(){m_name.clear(); m_size=0; m_encNuc=nullptr;}
		
	private:
		std::string m_name;
		unsigned char * m_encNuc=nullptr; //The encoded nucleotides sequence
		seqLength m_size;
};
	
class Reads {
	public:
		//Setters
		int readFile_fasta(const std::string& file);

		//Accessors
		unsigned int get_nbReads() const{return m_seqs.size();}
		bool display_read(const size_t id, std::ostream& os=std::cout) const;
		bool get_readName(const size_t id, std::string& str) const;
		//Get read (decoded)
		bool get_readSeq(const size_t id, std::string& str) const;


	protected:
		std::ostream& m_cout = std::cout;
		std::ostream& m_cerr = std::cerr;
		std::vector< Seq > m_seqs; //Sequences array

	private:

};

class Contigs {
	public:
		//Constructor need to have the set of all reads
		Contigs(const Reads* reads_set);

		//The destructor erased all contigs memory allowed
		~Contigs();

		//Accessors
			//Get the number of know contigs
		size_t get_nbContigs() const {return m_contigs_list.size();};
			//Get the real number of contigs (a sync will do "know" = "real")
		size_t get_nbContigs_real() const {return m_real_contigs_list.size();};
			//Get the sequence of a contig
		std::string get_seqContig(size_t contigID) const {
			return m_contigs_list[contigID]->c_seq;
		}; 

		//Setters
			//Merge 2 contigs (the 2nd is merged in the first)
			//The score is used to know if it is a normal merge (score >= 0)
			//or a reverse-complementary merge (score < 0)
		void merge_contigs(size_t contigID1, size_t contigID2, const int score);
			//Synchronize the "used" contigs list to the "real" list
			//To used on a "secure" state.
		void sync_contigs_list();


	private:
		//The set of the reads
		const Reads* m_reads_set;

		// A contig contains the ids of its reads and its sequence 
		struct Contig {
			std::set<size_t> c_reads_ids;
			std::string c_seq;
		};

		//This is the real list of actual contigs
		std::vector<Contig*> m_real_contigs_list;
		
		//This is the "buffer" list of contigs
		//The vector is sync to the real list when asked
		//If it's not sync, multiple elements can be linked to the same contig
		//It's used when editing an array (it should be sync when the calcs on the array is done)
		//Never doing a sync (when in "secure" state) after some contigs merge could occure a heavy ram load
		std::vector<Contig*> m_contigs_list;
};

#endif
