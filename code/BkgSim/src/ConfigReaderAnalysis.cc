#include "ConfigReaderAnalysis.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

ConfigReaderAnalysis::ConfigReaderAnalysis(std::string filename) : 
  ConfigReader(filename),
  _configName(""),
  _doEventLooper(false),
  _maxEvents(-1)
{}

void ConfigReaderAnalysis::setUp(){
				std::ifstream conffile;
				conffile.open(_filename.c_str(),std::ifstream::in);
				if(!conffile.is_open()){
								std::cerr << "Failed to open the config file for analysis: " << _filename;
								exit(1);
				}
				//Check that the config file isn't empty
				std::streampos fsize = 0;
				fsize = conffile.tellg();
				conffile.seekg( 0, std::ios::end );
				fsize = conffile.tellg() - fsize;
				if(fsize == 0){
								std::cerr << "Empty analysis config file specified";
								exit(2);
				}

				//Go back to the beginning, and make the all important 'temp' char
				conffile.seekg(0);
				char temp;

				while(!conffile.eof()){ //continue until we reach the end of the file
								bool newline(false);
								std::stringstream var,val;
								while(conffile.get(temp) && !newline){ //This is for each line of the file
												if(temp == '#'){ //Skip everything that appears after a #, as this is a comment
																while(conffile.get(temp)){
																				if(temp == '\n' || temp == ' '){
																								newline = true;
																								break;
																				}
																}
												}
												if(temp == '='){ //We have reached the end of a variable, figure out which one it is, and then fill the appropriate variable in the config file
																while(conffile.get(temp)){
																				if(temp == '\n' || temp == ' '){
																								newline = true;
																								break;
																				}
																				val << temp;
																}
												} else{
																var << temp;
												}
								}
								if(!conffile.eof()) conffile.unget(); //Otherwise we skip the first letter of new lines
								//Now fill variables as appropriate
								if(var.str() == "ConfigName"){
												setConfigName(val.str());
								}
								if(var.str() == "DoEventLooper"){
												setDoEventLooper(stoi(val.str()) != 0);
								}
								if(var.str() == "MaxEvents"){
												setMaxEvents(stoi(val.str()));
								}
								if(var.str() == "OutputFile"){
												setOutputFile(val.str());
								}

				}
}

void ConfigReaderAnalysis::setConfigName(const std::string configname){
  _configName = configname;
}

void ConfigReaderAnalysis::setOutputFile(const std::string outputfile){
  _outputFile = outputfile;
}

void ConfigReaderAnalysis::setDoEventLooper(const bool loop){
  _doEventLooper = loop;
}

void ConfigReaderAnalysis::setMaxEvents(const long int i){
  _maxEvents = i;
}

std::string ConfigReaderAnalysis::getConfigName() const{
  return _configName;
}

bool ConfigReaderAnalysis::getDoEventLooper() const{
  return _doEventLooper;
}

long int ConfigReaderAnalysis::getMaxEvents() const{
  return _maxEvents;
}
