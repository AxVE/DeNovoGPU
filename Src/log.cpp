#include "log.hpp"

#include <string>
#include <fstream>
#include <iostream>

using namespace std;

Log::Log(const string& log_output){
	if(log_output == "" || log_output == "cout" || log_output == "stdout"){
		m_mode=STDOUT;
		m_out.basic_ios<char>::rdbuf(cout.rdbuf());
	}
	else{
		m_mode=FILE;
		m_out.open(log_output);
	}
}

Log::~Log(){
	if(m_mode==FILE){
		m_out.close();
	}
}

void Log::write(const string& message){
	m_out << message << endl;
}
