/*
 * CellHits.cpp
 *
 *  Created on: Jan 11, 2016
 *      Author: schuea
 */

#include <iostream>
#include <sstream>

#include <bitset>
#include <vector>
#include <map>
#include <cmath>
#include <math.h>       /* atan2 */

#include "CellHits_class_new.h"

std::vector<uint64_t> CellHits::Get_CellID() const {
	return CellID;
}
std::vector<int> CellHits::Get_HitCount() const {
	return HitCount;
}
std::vector< float > CellHits::Get_HitPosition(char xyz) const {
	if (xyz != 'x' && xyz != 'y' && xyz != 'z') {
		std::cerr << "Input not correct! Has to be 'x', 'y' or 'z'!" << std::endl;
		exit(1);
	}
	else if (xyz == 'x') return HitPosition_x;
	else if (xyz == 'y') return HitPosition_y;
	else if (xyz == 'z') return HitPosition_z;
	else
	{
		std::cout << "Something weird happening with x y and z" << std::endl;
		exit(1);
	}
}
std::vector<int> CellHits::Get_Layer() const {
	return Layer;
}
uint64_t CellHits::CalculateLayer(uint64_t const id, Subdetector const & subdetector) {
  uint64_t LayerID64;

  LayerID64 = id << (64 - subdetector.getLengthBitLayer() - subdetector.getStartBitLayer());
  LayerID64 = LayerID64 >> (64 - subdetector.getLengthBitLayer());

  return LayerID64;
}
int CellHits::Get_NumberHitsPerLayer(int LayerNumber) {
	return Calculate_NumberHitsPerLayer(LayerNumber);
}

int CellHits::Check_CellID(uint64_t const id, float const x, float const y, float const z, Subdetector const & subdetector) {
	bool cell_exists(false);
	int vector_element(-1);
        //std::cout << id << ", " << x << ", " << y << ", " << z << std::endl;
	for (size_t i = 0; i < CellID.size(); ++i) {
		if (CellID.at(i) == id) {
			cell_exists = true;
			vector_element = i; //Check at which position in vector the ID is stored
			break;
		}
	}
	if (cell_exists) {
		HitCount.at(vector_element) += 1;
    return vector_element;
	} else {
		CellID.push_back(id);
		HitCount.push_back(1);
		HitPosition_x.push_back(x);
		HitPosition_y.push_back(y);
		HitPosition_z.push_back(z);
		Layer.push_back(CalculateLayer(id, subdetector));
    return CellID.size()-1;
	}
}

int CellHits::Calculate_NumberHitsPerLayer(int LayerNumber) {
	int NumberHitsPerLayer(0);
	for (size_t i = 0; i < Layer.size(); ++i){
		if (Layer.at(i) == LayerNumber) NumberHitsPerLayer += 1;
	}
	return NumberHitsPerLayer;
}

