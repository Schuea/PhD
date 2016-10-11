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
												cx(0), cy(0),
												xi(0),
												beta(0), alpha(0),
												sinalpha(0), cosalpha(0),
												position_prime_x(0), position_prime_y(0),
												position()
												{
																std::cout << "This should not be printed!" << std::endl;
								}
								Helix(float const Field) :	
												B(Field),
												px(0), py(0), pz(0),
												x0(0), y0(0), z0(0),
												charge(0), z(0),
												radius(0),
												cx(0), cy(0),
												xi(0),
												beta(0), alpha(0),
												sinalpha(0), cosalpha(0),
												position_prime_x(0), position_prime_y(0),
												position()
												{
																std::cout << "B(Field) = " << "B(" << Field <<") = " << B << std::endl;
								}

								~Helix(){
								}

								double* Get_position(double const pos_z);
                void Set_particlevalues(std::vector< double > const mom, float const particle_charge, std::vector< double > const origin);

				private:
								float B;
								double px, py, pz;
								double x0, y0, z0;
								float charge;
								double z;
								double radius;
				        double cx, cy;
				        double xi;
				
								double beta;
								double alpha;
								double sinalpha;
								double cosalpha;
								double position_prime_x;
								double position_prime_y;

								double position [3];

								void Calculate_radius();
								double Calculate_number_turn() const;
								double Calculate_xi() const;
								void Calculate_circlecenter();
};

#endif /* HELIX_H_*/
