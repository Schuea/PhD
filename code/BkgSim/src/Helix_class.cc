#include "Helix_class.h"

void Helix::Calculate_radius(){
				radius = sqrt(px*px+py*py)/(0.3*B);
				//std::cout << "radius = " << radius << std::endl;
}

double Helix::Calculate_number_turn() const{
			  //std::cout << "number_turn = " << number_turn << std::endl;
        return (z-z0)/(pz*2.0*M_PI/(0.3*B));
}

double Helix::Calculate_xi() const{
			  //std::cout << "xi = " << xi << std::endl;
				if (py/std::abs(py)*charge < 0 ){ //if the particle's circle is on the right hand side of the p_T vector, the angle starts off from pi and not 0.
          return charge*Calculate_number_turn()*(2.0*M_PI) + M_PI;
        } else{
          return charge*Calculate_number_turn()*(2.0*M_PI);
        }
}

void Helix::Calculate_circlecenter(){
				beta = M_PI*0.5 + py/std::abs(py)*charge*M_PI*0.5; //beta is the angle between the x-axis and the axis perpendicular to p_T in the xy-plane
			  //std::cout << "beta = " << beta << std::endl;
				cx = radius*cos(beta);//cx and cy are the x- and y-coordinates of the center of the circle that the helix performs in the xy-plane
				cy = radius*sin(beta);
}

std::vector< double > Helix::Calculate_position(){

				xi = Calculate_xi();

				position_prime_x = radius*cos(xi) + x0 + cx;
				position_prime_y = radius*sin(xi) + y0 + cy;

				//The actual position can be gained by rotating the coordinate system with the rotation matrix:
				position.push_back( position_prime_x*cosalpha - position_prime_y*sinalpha );
				position.push_back( position_prime_x*sinalpha + position_prime_y*cosalpha );
				position.push_back( z + z0 );

			  //if (std::abs(position.at(0)) > 4.0/1000.0 || std::abs(position.at(1)) > 4.0/1000.0){ //Only for the positions that are 4 mm away from the IP (0.0.0)
			  //				std::cout << "charge = " << charge << std::endl;
			  //				std::cout << "px = " << px << std::endl;
			  //				std::cout << "py = " << py << std::endl;
				//				std::cout << "radius = " << radius << std::endl;
			  //				std::cout << "cx = " << cx << std::endl;
			  //				std::cout << "cy = " << cy << std::endl;
			  //				std::cout << "number_turn = " << number_turn << std::endl;
			  //				std::cout << "xi = " << xi << std::endl;
			  //				std::cout << "position.at(0) = " << position.at(0) << std::endl;
			  //				std::cout << "position.at(1) = " << position.at(1) << std::endl;
			  //}
				return position;
}

std::vector< double > Helix::Get_position(double const pos_z){
				z = pos_z*0.001; //to convert mm into m (the histogramm will be drawn in mm, therefore z is given in mm)
				return Calculate_position();
}

void Helix::Set_particlevalues(std::vector< double > const mom, float const particle_charge, std::vector< double > const origin){
				px = mom.at(0);
				py = mom.at(1);
				pz = mom.at(2);
				charge = particle_charge;
				x0 = origin.at(0);
				y0 = origin.at(1);
				z0 = origin.at(2);
				
				//Everything that depends only on px and py can be calculated already here:
				Calculate_radius();
				Calculate_circlecenter();
				alpha = atan2(py,px); // alpha is the angle between x-axis and p_T vector in the xy-plane
        sinalpha = sin(alpha);
        cosalpha = cos(alpha);
}
