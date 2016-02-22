#include "LayerCodeInCellID_class.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <sstream>

int LayerCodeInCellID::FindLayer(std::string const CellID_, int const StartBin_layers, int const LengthBin_layers) {

	std::stringstream LayerID;

	//This for loop calculates the layer id
	//From a sring of 0's and 1's, e.g. 00001011010010
	//The StartBin is the first bin in the string we are interested in (when reading from right to left)
	//The LengthBin is the length of the string we are interested in
	//We read from left to right, but we specify the start position from right to left
	//There is a magic +1 in there because strings start at element 0.
	for (int i = CellID_.size() - (StartBin_layers + LengthBin_layers); i <= CellID_.size() - (StartBin_layers + 1);
			++i) {
		LayerID << CellID_.at(i);
	}

	std::bitset<64> LayerIDbit(LayerID.str());
	int tempLayer = -1;
	tempLayer = LayerIDbit.to_ulong();
	return tempLayer;
}

int LayerCodeInCellID::GetLayer(std::string const CellID, int const StartBin_layers, int const LengthBin_layers) {
	Layer = -1;
	Layer = FindLayer(CellID, StartBin_layers, LengthBin_layers);
	if (Layer >= 0) {
		return Layer;
	} else
		throw std::runtime_error("You didn't initialize the class with an ID!");
}

int LayerCodeInCellID::GetLayer(long long const cellid, int const StartBin_layers, int const LengthBin_layers) {

	std::bitset<64> cellidbit(cellid);
	Layer = -1;
	Layer = FindLayer(cellidbit.to_string(), StartBin_layers, LengthBin_layers);
	if (Layer >= 0) {
		return Layer;
	} else
		throw std::runtime_error("You didn't initialize the class with an ID!");
}

