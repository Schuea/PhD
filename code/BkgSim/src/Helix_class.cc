#include "Helix_class.h"

void Helix::Calculate_radius(){
								std::cout << __FILE__ << " " << __LINE__ << std::endl;
								std::cout << "B = " << B << std::endl;
								std::cout << "px = " << px << std::endl;
								std::cout << "py = " << py << std::endl;
				radius = sqrt(px*px+py*py)/(0.3*B);
				std::cout << "radius = " << radius << std::endl;
}

void Helix::Calculate_number_turn(){
								std::cout << __FILE__ << " " << __LINE__ << std::endl;
								std::cout << "pz = " << pz << std::endl;
				number_turn = (z-z0)/(pz*2*M_PI/(0.3*B));
			  std::cout << "number_turn = " << number_turn << std::endl;
}

void Helix::Calculate_xi(){
				Calculate_number_turn();
				xi = number_turn*(2*M_PI);
			  std::cout << "xi = " << xi << std::endl;
}

void Helix::Calculate_circlecenter(){
				double alpha, beta = 0.0;
				alpha = atan2(py,px);
				beta = alpha + charge*M_PI/2;
				cx = radius*cos(beta);
				cy = radius*sin(beta);
}

std::vector< double > Helix::Calculate_position(){

				Calculate_radius();
				radius = radius * charge;
				Calculate_xi();

				Calculate_circlecenter();
				std::cout << "cx = " << cx << std::endl;
				std::cout << "cy = " << cy << std::endl;

				std::vector< double > position;
				position.push_back( radius*cos(xi) + x0 + cx);
				position.push_back( radius*sin(xi) + y0 + cy);
				position.push_back( z + z0 );
				return position;
}

std::vector< double > Helix::Get_position(std::vector< double > mom, float particle_charge, std::vector< double > origin, double pos_z){
				px = mom.at(0);
				py = mom.at(1);
				pz = mom.at(2);
				charge = particle_charge;
				x0 = origin.at(0);
				y0 = origin.at(1);
				z0 = origin.at(2);
				z = pos_z/1000.0; //to convert mm into m
			
				//if (x0>0 || y0>0 ||z0>0){
								std::cout << __FILE__ << " " << __LINE__ << std::endl;
								std::cout << "mom.at(0) = " << mom.at(0) << std::endl;
								std::cout << "mom.at(1) = " << mom.at(1) << std::endl;
								std::cout << "mom.at(2) = " << mom.at(2) << std::endl;
				//}
				
				return Calculate_position();

}
