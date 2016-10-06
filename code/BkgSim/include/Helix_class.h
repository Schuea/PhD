//HELIX_H_
#ifndef HELIX_H_
#define HELIX_H_

#include <cmath>

#include <iostream>
#include <vector>

class Helix {

				public: 
								Helix():
												B(0),
												px(0), py(0), pz(0),
												x0(0), y0(0), z0(0),
												charge(0), z(0),
												radius(0), number_turn(0),
												xi(0), cx(0), cy(0),
												alpha(0){
												//position(){
																std::cout << "This should not be printed!" << std::endl;
								}
								Helix(float Field) :	
												B(Field),
												px(0), py(0), pz(0),
												x0(0), y0(0), z0(0),
												charge(0), z(0),
												radius(0), number_turn(0),
												xi(0), cx(0), cy(0),
												alpha(0){
												//position(){
																std::cout << "B(Field) = " << "B(" << Field <<") = " << B << std::endl;
								}

								~Helix(){
								}

								std::vector< double > Get_position(std::vector< double > mom, float particle_charge, std::vector< double > origin, double pos_z);

				private:
								float B;
								double px, py, pz;
								double x0, y0, z0;
								float charge;
								double z;
								double radius;
								double number_turn;
								double xi; 
				        double cx, cy;
								double alpha;

								//std::vector< double > position;

								void Calculate_radius();
								void Calculate_number_turn();
								void Calculate_xi();
								void Calculate_circlecenter();
								std::vector<double> Calculate_position();
};

#endif /* HELIX_H_*/
