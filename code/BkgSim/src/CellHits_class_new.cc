/*
 * CellHits.cpp
 *
 *  Created on: Jan 11, 2016
 *      Author: schuea
 */

#include <iostream>

#include <bitset>
#include <vector>
#include <map>
#include <cmath>
#include <math.h>       /* atan2 */

#include "CellHits_class_new.h"

std::vector<long long int> CellHits::Get_CellID() const {
	return CellID;
}
std::vector<int> CellHits::Get_HitCount() const {
	return HitCount;
}
std::vector< int > CellHits::Get_Layer(Subdetector det, long long int const id) {
  std::bitset<64> cellidbit(id);
  std::string CellID_ = cellidbit.to_string();
  Layer = -1;
  std::stringstream LayerID;

  //This for loop calculates the layer id
  //From a sring of 0's and 1's, e.g. 00001011010010
  //The StartBin is the first bin in the string we are interested in (when reading from right to left)
  //The LengthBin is the length of the string we are interested in
  //We read from left to right, but we specify the start position from right to left
  //There is a magic +1 in there because strings start at element 0.
  for (int i = CellID_.size() - (det.getStartBitLayer + det.getLengthBitLayers); i <= CellID_.size() - (det.getStartBitLayers + 1);
      ++i) {
    LayerID << CellID_.at(i);
  }

  std::bitset<64> LayerIDbit(LayerID.str());
  Layer = LayerIDbit.to_ulong();

  return Layer;
}
int CellHits::Get_NumberHitsPerLayer(int LayerNumber) {
	return Calculate_NumberHitsPerLayer(LayerNumber);
}

void CellHits::Check_CellID(long long int const id, float const x, float const y, float const z) {
	bool cell_exists(false);
	int vector_element(-1);
	for (size_t i = 0; i < CellID.size(); ++i) {
		if (CellID.at(i) == id) {
			cell_exists = true;
			vector_element = i; //Check at which position in vector the ID is stored
			break;
		}
	}
	if (cell_exists) {
		HitCount.at(vector_element) += 1;
	} else {
		CellID.push_back(id);
		HitCount.push_back(1);
		HitPosition_x.push_back(x);
		HitPosition_y.push_back(y);
		HitPosition_z.push_back(z);
		Layer.push_back(GetLayer(subdetector, id));
	}
}

int CellHits::Calculate_NumberHitsPerLayer(int LayerNumber) {
	int NumberHitsPerLayer(0);
	for (size_t i = 0; i < Layer.size(); ++i){
		if (Layer.at(i) == LayerNumber) NumberHitsPerLayer += 1;
	}
	return NumberHitsPerLayer;
}

