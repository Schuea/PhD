#include "ConfigReader.h"

#include <iostream>

ConfigReader::ConfigReader(std::string filename) :
_filename(filename){
}

void ConfigReader::setFilename(const std::string filename){
  _filename = filename;
}

std::string ConfigReader::getFilename() const{
  return _filename;
}
