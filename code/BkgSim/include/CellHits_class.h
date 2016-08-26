/*
 * CellHits_class.h
 *
 *  Created on: Jan 11, 2016
 *      Author: schuea
 */

#ifndef CELLHITSCLASS_H_
#define CELLHITSCLASS_H_

#include "Subdetector_class.h"

class CellHits {
public:
  CellHits(){}
	CellHits(Subdetector* SubDetector_) :
			SubDetector(SubDetector_), 
      CellID(), 
      HitCount(), 
			Particle(),
      ParticleCount(), 
			HitTime(),
			HitMomentum_z(),
			HitMomentum_T(),
			HitPosition_x(), 
			HitPosition_y(), 
			HitPosition_z(), 
			Layer(),
      Position_Radius(),
      Position_Phi(),
      AverageOccupancy_Rad(),
      AverageOccupancy_Phi(),
      BunchNumber(0) {
	}
	~CellHits() {
	}

	std::vector< long long int > Get_CellID() const;
	std::vector< int > Get_HitCount() const;
	std::vector< int > Get_Particle() const;
	std::vector< int > Get_ParticleCount() const;
	std::vector< float > Get_HitTime() const;
	std::vector< float > Get_HitPosition(char xyz) const;
	std::vector< float > Get_HitMomentum(char zT) const;
	std::vector< int > Get_Layer() const;
	int Get_NumberHitsPerLayer(int LayerNumber);
	std::vector< float > Get_Position_Radius() const;
	std::vector< float > Get_Position_Phi() const;
	std::map<int, std::pair<std::vector<int>, std::pair<float, float> > > Get_AverageOccupancy_Rad() const;
	std::map<int, std::pair<std::vector<int>, std::pair<float, float> > > Get_AverageOccupancy_Phi() const;
	int Get_BunchNumber() const;
	void Set_BunchNumber(int const bunchnumber);

	int Calculate_NumberHitsPerLayer(int LayerNumber);
	void Check_ParticleID(long long int const id, int const particle_id, float const time, float const x, float const y, float const z, float const p_T, float const p_z ); 
	void Check_CellID(long long int const id, float const x, float const y, float const z);
	void Check_Rad_Position();
	void Check_Phi_Position();
	void Calculate_Average(std::map<int, std::pair<std::vector<int>, std::pair<float, float> > > & AverageMap);

protected:
	Subdetector* SubDetector;
	std::vector<long long int> CellID;
	std::vector<int> HitCount;
	std::vector<int> Particle;
	std::vector<int> ParticleCount;
	std::vector<float> HitTime;
	std::vector<float> HitMomentum_z;
	std::vector<float> HitMomentum_T;
	std::vector<float> HitPosition_x;
	std::vector<float> HitPosition_y;
	std::vector<float> HitPosition_z;
	std::vector<int> Layer;
	std::vector<float> Position_Radius;
	std::vector<float> Position_Phi;
	std::map<int, std::pair<std::vector<int>, std::pair< float, float > > > AverageOccupancy_Rad;
	std::map<int, std::pair<std::vector<int>, std::pair< float, float > > > AverageOccupancy_Phi;
	int BunchNumber;
};

#endif /* CELLHITSCLASS_H_ */
