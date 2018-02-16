/*
 * CellHits_class.h
 *
 *  Created on: Jan 11, 2016
 *      Author: schuea
 */

#ifndef CELLHITSCLASS_NEW_H_
#define CELLHITSCLASS_NEW_H_

#include "Subdetector_class_new.h"

class CellHits {
public:
	CellHits(Subdetector *subdetector) :
      SubDetector(subdetector),
      CellID(), 
      HitCount(), 
			HitPosition_x(), 
			HitPosition_y(), 
			HitPosition_z(), 
			Layer()
  {
	}
	~CellHits() {
	}

	std::vector< uint64_t > Get_CellID() const;
	std::vector< int > Get_HitCount() const;
	std::vector< float > Get_HitPosition(char xyz) const;
	std::vector< int > Get_Layer() const;
  uint64_t CalculateLayer(uint64_t const id, Subdetector const & subdetector);
	int Get_NumberHitsPerLayer(int LayerNumber);

	int Calculate_NumberHitsPerLayer(int LayerNumber);
	void Check_CellID(uint64_t const id, float const x, float const y, float const z, Subdetector const & subdetector);

protected:
	Subdetector* SubDetector;
	std::vector<uint64_t> CellID;
	std::vector<int> HitCount;
	std::vector<float> HitPosition_x;
	std::vector<float> HitPosition_y;
	std::vector<float> HitPosition_z;
	std::vector<int> Layer;
};

#endif /* CELLHITSCLASS_H_ */
