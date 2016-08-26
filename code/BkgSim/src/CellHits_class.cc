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

#include "CellHits_class.h"

std::vector<long long int> CellHits::Get_CellID() const {
	return CellID;
}
std::vector<int> CellHits::Get_HitCount() const {
	return HitCount;
}
std::vector<int> CellHits::Get_Particle() const {
	return Particle;
}
std::vector<int> CellHits::Get_ParticleCount() const {
	return ParticleCount;
}
std::vector<float> CellHits::Get_HitTime() const {
	return HitTime;
}
std::vector< float > CellHits::Get_HitPosition(char xyz) const {
	if (xyz != 'x' && xyz != 'y' && xyz != 'z') {
		std::cerr << "Input not correct! Has to be 'x', 'y' or 'z'!" << std::endl;
		exit(1);
	}
	else if (xyz == 'x') return HitPosition_x;
	else if (xyz == 'y') return HitPosition_y;
	else if (xyz == 'z') return HitPosition_z;
}
std::vector< float > CellHits::Get_HitMomentum(char zT) const {
	if (zT != 'T' && zT != 'z') {
		std::cerr << "Input not correct! Has to be 'T' or 'z'!" << std::endl;
		exit(1);
	}
	else if (zT == 'T') return HitMomentum_T;
	else if (zT == 'z') return HitMomentum_z;
}
std::vector< int > CellHits::Get_Layer() const {
	return Layer;
}
int CellHits::Get_NumberHitsPerLayer(int LayerNumber) {
	return Calculate_NumberHitsPerLayer(LayerNumber);
}
std::vector<float> CellHits::Get_Position_Radius() const {
	return Position_Radius;
}
std::vector<float> CellHits::Get_Position_Phi() const {
	return Position_Phi;
}
std::map<int, std::pair<std::vector<int>, std::pair<float, float> > > CellHits::Get_AverageOccupancy_Rad() const {
	return AverageOccupancy_Rad;
}
std::map<int, std::pair<std::vector<int>, std::pair<float, float> > > CellHits::Get_AverageOccupancy_Phi() const {
	return AverageOccupancy_Phi;
}
int CellHits::Get_BunchNumber() const {
	return BunchNumber;
}
void CellHits::Set_BunchNumber(int const bunchnumber) {
	BunchNumber = bunchnumber;
}

void CellHits::Check_ParticleID(long long int const id, int const particle_id, float const time, float const x, float const y, float const z, float const p_T, float const p_z ) {
	bool particleID_exists(false);
	int vector_element(-1);
	for (size_t i = 0; i < Particle.size(); ++i) {
		if (Particle.at(i) == particle_id) {
			particleID_exists = true;
			vector_element = i; //Check at which position in vector the ID is stored
			break;
		}
	}
	if (particleID_exists) {
		ParticleCount.at(vector_element) += 1;
	} else {
		Particle.push_back(particle_id);
		ParticleCount.push_back(1);
		CellID.push_back(id);
		HitTime.push_back(time);
		HitMomentum_z.push_back(p_z);
		HitMomentum_T.push_back(p_T);
		HitPosition_x.push_back(x);
		HitPosition_y.push_back(y);
		HitPosition_z.push_back(z);
		Layer.push_back(SubDetector->GetLayer(id));
	}
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
		Layer.push_back(SubDetector->GetLayer(id));
		Position_Radius.push_back(sqrt(pow(x, 2) + pow(y, 2)));
		Position_Phi.push_back(atan2(y, x));
	}
}

int CellHits::Calculate_NumberHitsPerLayer(int LayerNumber) {
	int NumberHitsPerLayer(0);
	for (size_t i = 0; i < Layer.size(); ++i){
		if (Layer.at(i) == LayerNumber) NumberHitsPerLayer += 1;
	}
	return NumberHitsPerLayer;
}

void CellHits::Check_Rad_Position() {

	for (size_t j = 0; j < Position_Radius.size(); ++j) {

		if (SubDetector->GetName() == "SiVertexEndcap" || SubDetector->GetName() == "SiVertexBarrel") {
			if (Position_Radius.at(j) < 10)
				AverageOccupancy_Rad[5].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 10 && Position_Radius.at(j) < 20)
				AverageOccupancy_Rad[15].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 20 && Position_Radius.at(j) < 30)
				AverageOccupancy_Rad[25].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 30 && Position_Radius.at(j) < 40)
				AverageOccupancy_Rad[35].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 40 && Position_Radius.at(j) < 50)
				AverageOccupancy_Rad[45].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 50 && Position_Radius.at(j) < 60)
				AverageOccupancy_Rad[55].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 60 && Position_Radius.at(j) < 70)
				AverageOccupancy_Rad[65].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 70 && Position_Radius.at(j) < 80)
				AverageOccupancy_Rad[75].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 80)
				AverageOccupancy_Rad[85].first.push_back(HitCount.at(j));
		} else {
			//if (found == std::string::npos) {
			if (Position_Radius.at(j) < 200)
				AverageOccupancy_Rad[100].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 200 && Position_Radius.at(j) < 400)
				AverageOccupancy_Rad[300].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 400 && Position_Radius.at(j) < 600)
				AverageOccupancy_Rad[500].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 600 && Position_Radius.at(j) < 800)
				AverageOccupancy_Rad[700].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 800 && Position_Radius.at(j) < 1000)
				AverageOccupancy_Rad[900].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 1000 && Position_Radius.at(j) < 1200)
				AverageOccupancy_Rad[1100].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 1200 && Position_Radius.at(j) < 1400)
				AverageOccupancy_Rad[1300].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 1400 && Position_Radius.at(j) < 1600)
				AverageOccupancy_Rad[1500].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 1600 && Position_Radius.at(j) < 1800)
				AverageOccupancy_Rad[1700].first.push_back(HitCount.at(j));
			if (Position_Radius.at(j) >= 1800 && Position_Radius.at(j) < 2000)
				AverageOccupancy_Rad[1900].first.push_back(HitCount.at(j));
		}
	}
	Calculate_Average(AverageOccupancy_Rad);
}
void CellHits::Check_Phi_Position() {
	for (size_t j = 0; j < Position_Phi.size(); ++j) {

		 if (Position_Phi.at(j) < 0 && Position_Phi.at(j) > -0.4)
		 AverageOccupancy_Phi[-2].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -0.4 && Position_Phi.at(j) > -0.8)
		 AverageOccupancy_Phi[-6].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -0.8 && Position_Phi.at(j) > -1.2)
		 AverageOccupancy_Phi[-10].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -1.2 && Position_Phi.at(j) > -1.6)
		 AverageOccupancy_Phi[-14].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -1.6 && Position_Phi.at(j) > -2.0)
		 AverageOccupancy_Phi[-18].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -2.0 && Position_Phi.at(j) > -2.4)
		 AverageOccupancy_Phi[-22].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -2.4 && Position_Phi.at(j) > -2.8)
		 AverageOccupancy_Phi[-26].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) <= -2.8 && Position_Phi.at(j) > -3.2)
		 AverageOccupancy_Phi[-30].first.push_back(HitCount.at(j));

		 if (Position_Phi.at(j) >= 0 && Position_Phi.at(j) < 0.4)
		 AverageOccupancy_Phi[2].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 0.4 && Position_Phi.at(j) < 0.8)
		 AverageOccupancy_Phi[6].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 0.8 && Position_Phi.at(j) < 1.2)
		 AverageOccupancy_Phi[10].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 1.2 && Position_Phi.at(j) < 1.6)
		 AverageOccupancy_Phi[14].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 1.6 && Position_Phi.at(j) < 2.0)
		 AverageOccupancy_Phi[18].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 2.0 && Position_Phi.at(j) < 2.4)
		 AverageOccupancy_Phi[22].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 2.4 && Position_Phi.at(j) < 2.8)
		 AverageOccupancy_Phi[26].first.push_back(HitCount.at(j));
		 if (Position_Phi.at(j) >= 2.8 && Position_Phi.at(j) < 3.2)
		 AverageOccupancy_Phi[30].first.push_back(HitCount.at(j));
	}
	Calculate_Average(AverageOccupancy_Phi);
}
void CellHits::Calculate_Average(std::map<int, std::pair<std::vector<int>, std::pair<float, float> > > & AverageMap) {
	for (auto &iterator : AverageMap) {
		float average = 0;
		float stddev = 0;

		for (size_t i = 0; i < iterator.second.first.size(); ++i) {
			average += iterator.second.first.at(i);
		}
		if (iterator.second.first.size() != 0) average /= float(iterator.second.first.size());
		else average = 0;

		for (size_t i = 0; i < iterator.second.first.size(); ++i) {
			stddev += pow(iterator.second.first.at(i) - average,2);
		}
		if (iterator.second.first.size() != 0){
			stddev /= float(iterator.second.first.size());
			stddev = sqrt(stddev);
		}
		else stddev = 0;

		iterator.second.second.first = average;
		iterator.second.second.second = stddev;
	}
}
