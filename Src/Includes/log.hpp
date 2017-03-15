#ifndef LOG_HPP
#define LOG_HPP

#include <fstream>
#include <string>

class Log {
	public:
		//Possibles modes of output
		enum Mode { STDOUT, FILE };

		//Constructor need the output destination of logs. 'cout' or 'stdout' will make use the other outputs. Either, it's a file.
		Log(const std::string& log_output); //Construtor need the out

		//Destructor must be defined to close the output file (if one was used
		~Log();

		//Function to send message to the output
		void write(const std::string& message);
	private:
		Mode m_mode; //Actual mode of output usedprivate:
		//Output destination parameters
		std::ofstream m_out;

};
#endif
