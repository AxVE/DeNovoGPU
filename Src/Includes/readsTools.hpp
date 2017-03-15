#ifndef READSTOOLS_HPP
#define READSTOOLS_HPP

#include <string>
#include <map>
#include <vector>


struct ReadsTools{
	public:
	/* The nucleotides list assigned to a specific id */
	const static std::map<char,uint8_t> nuc;
	
	/* The scoring matrix of the matching (and unmatching)
	between each nucleotide. It's a 2D array of size 5 (*2)
	The id (0,1,2,3,4)  represents A,T,G,C,X respectivly
	*/
	const static std::vector< std::vector< int8_t > > scoreMatrix;
	const static int8_t bestScoreMatrix;
	
	/* The cost of a gap */
	const static int gap_score;



	/* Calculate the score of matching between 2 reads.
	The value of a score is absolute. A negative score
	indicates than a read is reversed. */
	static const int scoreCmpReads(const std::string& read1, const std::string& read2);

	/* Return the best allignment between 2 sequences (without reversing, it's direct comparaison).
	It's used needleman / smith-watermann */
	static const std::string allignSeq(const std::string& seq1, const std::string& seq2);

	private:
	/* Calculate score between 2 reads (without reversing on of them). The score is centered using the mean of the size of the reads 
	((read1.size()+read2.size())/2). */
	static const int scoreMatchingReads(const std::string& read1, const std::string& read2);
};

#endif
