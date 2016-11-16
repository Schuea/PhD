#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Subdetector.h"

using namespace std;

vector< double > ConvertCSVToVectorDoubles(string const csv);

long long int MakeNewCellID(double const x, double const y, Subdetector const & component);

Subdetector::Subdetector(){
}

Subdetector::Subdetector(string const subdetector_config_file){
  std::ifstream file(subdetector_config_file, std::ifstream::in);
  if(!file.is_open()){
    cerr << "Error opening subdetector config file: " << subdetector_config_file << endl;
    exit(2);
  }
  char variable[256], value[256];
  double const PI(3.1415926535);
  while(!file.eof()){
    file.getline(variable, 256, '=');
    file.getline(value, 256);
    if(string(variable) == "numberOfLayers") setNumberOfLayers(stoi(string(value)));
    if(string(variable) == "startBitLayer") setStartBitLayer(stoi(string(value)));
    if(string(variable) == "lengthBitLayer") setLengthBitLayer(stoi(string(value)));
    if(string(variable) == "rMin") setRMin(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "rMax") setRMax(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "zHalf") setZHalf(ConvertCSVToVectorDoubles(string(value)));
    if(string(variable) == "length") setLength(value);
    //Add something that checks values of zHalf and length, if one is set and the other isn't, set the other based on the first
    if(string(variable) == "shape"){
      //Depending on the shape, this will have multiple layers which have to be taken into account differently
      if(string(value) == "barrel") getRLayer.push_back(2*PI*getRMin()*length);
      //Finish this for the other shapes
      //Set the area using rLayer*shape_number_of_sides*length or equivilent formula
      //Number of cells is just area divided by cellSizeArea
    }
  }
  
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

double Subdetector::getLength() const{
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

void Subdetector::setLength(double const length){
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

long long int MakeNewCellID(double const x, double const y, Subdetector const & component){
  int newX = static_cast<int>(x/component.getCellSizeArea()); //Check if Cell Size Area is the same as Cell Dimension
  int newY = static_cast<int>(y/component.getCellSizeArea()); //Check if Cell Size Area is the same as Cell Dimension
  if(x >= 0) ++newX;
  if(y >= 0) ++newY;
  bitset<32> bitY = newY;
  newY = 0;
  for(int i = 0; i < 31; ++i){
    newY += bitY[i]*pow(2,i);
  }
  return newX << 32 | newY;
}

vector< double > ConvertCSVToVectorDoubles(string const csv){
  int startpos(0);
  vector< double > result;
  while(csv.find(',',startpos) != string::npos){
    result.push_back(stof(csv.substr(startpos,csv.find(',',startpos) - startpos)));
    startpos = csv.find(',',startpos) + 1;
  }
  result.push_back(stof(csv.substr(startpos,csv.find(',',startpos) - startpos)));
  startpos = csv.find(',',startpos) + 1;
  return result;
}

int main(){
  Subdetector component("SiTrackerBarrel.cfg");
  return 0;
}


