//HELIX_H_
#ifndef HELIX_H_
#define HELIX_H_

#include <math.h>

#include <iostream>
#include <vector>

class Helix {

				public: 
								Helix(){
								}
								Helix(double momx, double momy, double momz, int particle_charge, std::vector< float > origin, float Field):
												px(momx), py(momy), pz(momz),
												charge(particle_charge),
												x0(origin.at(0)), y0(origin.at(1)), z0(origin.at(2)),
												B(Field){
								}

								~Helix(){
								}

								std::vector< double > Get_position(float pos_z);

				private:
								double px, py, pz;
								int charge;
								float x0, y0, z0;
								float z;
								float B;
								double length;
								double number_turn;
								double xi; 
								std::vector< double > position;

								void Calculate_number_turn();
								void Calculate_xi();
								void Calculate_position();
};

#endif /* HELIX_H_*/
