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
												radius(0), 
												cx(0), cy(0){
																std::cout << "This should not be printed!" << std::endl;
								}
								Helix(float const Field) :	
												B(Field),
												px(0), py(0), pz(0),
												x0(0), y0(0), z0(0),
												charge(0), z(0),
												radius(0),
												cx(0), cy(0){
																std::cout << "B(Field) = " << "B(" << Field <<") = " << B << std::endl;
								}

								~Helix(){
								}

								std::vector< double > Get_position(double const pos_z);
                void Set_particlevalues(std::vector< double > const mom, float const particle_charge, std::vector< double > const origin);

				private:
								float B;
								double px, py, pz;
								double x0, y0, z0;
								float charge;
								double z;
								double radius;
				        double cx, cy;

								void Calculate_radius();
								double Calculate_number_turn() const;
								double Calculate_xi() const;
								void Calculate_circlecenter();
								std::vector<double> Calculate_position();
				
								double alpha;
								double sinalpha;
								double cosalpha;
};

#endif /* HELIX_H_*/
