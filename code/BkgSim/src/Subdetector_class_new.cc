#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>

#include "Subdetector_class_new.h"

using namespace std;

Subdetector::Subdetector(){
}

Subdetector::Subdetector(string const subdetector_config_file){
  std::ifstream file(subdetector_config_file, std::ifstream::in);
  if(!file.is_open()){
    cerr << "Error opening subdetector config file: " << subdetector_config_file << endl;
    exit(2);
  }
  char variable[256], value[256];
  while(!file.eof()){
    file.getline(variable, 256, '=');
    file.getline(value, 256);
    if(string(variable) == "name") setName(string(value));
    if(string(variable) == "numberOfLayers") setNumberOfLayers(stoi(string(value)));
    if(string(variable) == "startBitLayer") setStartBitLayer(stoi(string(value)));
    if(string(variable) == "lengthBitLayer") setLengthBitLayer(stoi(string(value)));
    if(string(variable) == "rMin") setRMin(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "rMax") setRMax(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "zHalf") setZHalf(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "length") setLength(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "cellSizeX") setCellSizeX(stof(string(value)));
    if(string(variable) == "cellSizeY") setCellSizeY(stof(string(value)));

    if(string(variable) == "shape"){
      //Depending on the shape, this will have multiple layers which have to be taken into account differently
      for (int Layer = 0; Layer < getNumberOfLayers(); ++Layer){
        if(string(value) == "octagon-barrel"){//MuonBarrel
          getRLayer().push_back( sin(2*M_PI/8)*(getRMin().at(Layer)+(getRMax().at(Layer)-getRMin().at(Layer))/getNumberOfLayers())*(Layer+1)/sin((M_PI-2*M_PI/8)/2) );
          getArea().push_back( getRLayer().back()*8*getLength().at(Layer) );
          getNumberOfCells().push_back( (int)(getArea().back()/getCellSizeArea()) );
        }
        if(string(value) == "octagon-endcap"){//MuonEndcap
          getRLayer().push_back( getRMax().at(Layer) );
          getArea().push_back( pow(getRLayer().back(),2)*2*sqrt(2) - pow(getRMin().at(Layer),2)*8*tan(M_PI/8) );
          getNumberOfCells().push_back( (int)(getArea().back()/getCellSizeArea()) );
        }
        if(string(value) == "dodecagon-barrel"){//Ecal+HcalBarrel
          getRLayer().push_back( sin(2*M_PI/12)*(getRMin().at(Layer)+(getRMax().at(Layer)-getRMin().at(Layer))/getNumberOfLayers())*(Layer+1)/sin((M_PI-2*M_PI/12)/2) );
          getArea().push_back( getRLayer().back()*12*getLength().at(Layer) );
          getNumberOfCells().push_back( (int)(getArea().back()/getCellSizeArea()) );
        }
        if(string(value) == "dodecagon-endcap"){//Ecal+HcalEndcaps
          getRLayer().push_back( getRMax().at(Layer) );
          getArea().push_back( 3*pow(getRLayer().back(),2) - pow(getRMin().at(Layer),2)*12*(2-sqrt(3)) );
          getNumberOfCells().push_back( (int)(getArea().back()/getCellSizeArea()) );
        }
        if(string(value) == "circle-barrel"){//SiVertex+TrackerBarrel
          getRLayer().push_back( getRMax().at(Layer) );
          getArea().push_back( 2*M_PI*getRMin().at(Layer)*2*getLength().at(Layer) );
          getNumberOfCells().push_back( (int)(getArea().back()/getCellSizeArea()) );
        }
        if(string(value) == "circle-endcap"){//LumiCal,BeamCal,SiVertex+TrackerEndcap
          getRLayer().push_back( getRMax().at(Layer) );
          getArea().push_back( M_PI*(pow(getRLayer().back(),2) - pow(getRMin().at(Layer),2)) );
          getNumberOfCells().push_back( (int)(getArea().back()/getCellSizeArea()) );
        }
      }
    }
  }
  //Check the input parameters for zHalf and length:
  if (getLength().size() == 0 && getZHalf().size() != 0){
    for (size_t i=0; i < getZHalf().size(); ++i){
      getLength().at(i) = getZHalf().at(i)*2.0;
    }
  }
  else if (getLength().size() != 0 && getZHalf().size() == 0){
    for (size_t i=0; i < getLength().size(); ++i){
      getZHalf().at(i) = getLength().at(i)/2.0;
    }
  }
  else{
    if (getLength().size() == getZHalf().size()){
      for (size_t i=0; i < getLength().size(); ++i){
        if (getZHalf().at(i) < getLength().at(i)/2.0 + numeric_limits<double>::epsilon() &&
            getZHalf().at(i) > getLength().at(i)/2.0 - numeric_limits<double>::epsilon() ){
          cerr << "The given numbers for zHalf and length are not in agreement!" << endl; 
          cerr << "Exited in file " << __FILE__ << " on line " << __LINE__ << endl; 
          exit(-1);
        }
      }
    }
    else{
          cerr << "The number of input parameters for zHalf and length is not the same!" << endl; 
          cerr << "Exited in file " << __FILE__ << " on line " << __LINE__ << endl; 
          exit(-1);
    }
  }
  
}

std::string Subdetector::getName() const{
  return _name;
}

int Subdetector::getNumberOfLayers() const{
  return _numberOfLayers;
}

int Subdetector::getStartBitLayer() const{
  return _startBitLayer;
}

int Subdetector::getLengthBitLayer() const{
  return _lengthBitLayer;
}

vector< double > Subdetector::getRMin() const{
  return _rMin;
}

vector< double > Subdetector::getRMax() const{
  return _rMax;
}

vector< double > Subdetector::getZHalf() const{
  return _zHalf;
}

double Subdetector::getCellSizeX() const{
  return _cellSizeX;
}

double Subdetector::getCellSizeY() const{
  return _cellSizeY;
}

double Subdetector::getCellSizeArea() const{
  return _cellSizeArea;
}

vector< double > Subdetector::getLength() const{
  return _length;
}

std::vector< double > Subdetector::getRLayer() const{
  return _rLayer;
}

std::vector< double > Subdetector::getArea() const{
  return _area;
}

std::vector< double > Subdetector::getNumberOfCells() const{
  return _numberOfCells;
}

void Subdetector::setName(std::string const name){
  _name = name;
}

void Subdetector::setNumberOfLayers(int const number_of_layers){
  _numberOfLayers = number_of_layers;
}

void Subdetector::setStartBitLayer(int const start_bit_layer){
  _startBitLayer = start_bit_layer;
}

void Subdetector::setLengthBitLayer(int const length_bit_layer){
  _lengthBitLayer = length_bit_layer;
}

void Subdetector::setRMin(vector< double > const r_min){
  _rMin = r_min;
}

void Subdetector::setRMax(vector< double > const r_max){
  _rMax = r_max;
}

void Subdetector::setZHalf(vector< double > const z_half){
  _zHalf = z_half;
}

void Subdetector::setCellSizeX(double const cell_size_x){
  _cellSizeX = cell_size_x;
}

void Subdetector::setCellSizeY(double const cell_size_y){
  _cellSizeY = cell_size_y;
}

void Subdetector::setCellSizeArea(double const cell_size_area){
  _cellSizeArea = cell_size_area;
}

void Subdetector::setLength(std::vector< double > const length){
  _length = length;
}

void Subdetector::setRLayer(vector< double > const r_layer){
  _rLayer = r_layer;
}

void Subdetector::setArea(vector< double > const area){
  _area = area;
}

void Subdetector::setNumberOfCells(vector< double > const number_of_cells){
  _numberOfCells = number_of_cells;
}

vector< double > Subdetector::ConvertCSVToVectorDoubles(string const csv){
  int startpos(0);
  vector< double > result;
  while(csv.find(',',startpos) != string::npos){
    result.push_back(stof(csv.substr(startpos,csv.find(',',startpos) - startpos)));
    startpos = csv.find(',',startpos) + 1;
  }
  result.push_back(stof(csv.substr(startpos,csv.find(',',startpos) - startpos)));
  startpos = csv.find(',',startpos) + 1;
  
  FillUpVector( result );

  return result;
}

void Subdetector::FillUpVector( std::vector< double > vec ){
    if (vec.size() == 1){//if only one number is given, fill the vector with the same number for the rest of the layers
      for (int l=0; l < getNumberOfLayers()-1; ++l){
        vec.push_back( vec.at(0) );
      }
    }
}
