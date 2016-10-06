#include "Helix_class.h"

void Helix::Calculate_radius(){
				radius = sqrt(px*px+py*py)/(0.3*B);
				//std::cout << "radius = " << radius << std::endl;
}

void Helix::Calculate_number_turn(){
				number_turn = (z-z0)/(pz*2*M_PI/(0.3*B));
			  //std::cout << "number_turn = " << number_turn << std::endl;
}

void Helix::Calculate_xi(){
				Calculate_number_turn();
				xi = number_turn*(2*M_PI);
				xi = xi*charge;
			  //std::cout << "xi = " << xi << std::endl;
}

void Helix::Calculate_circlecenter(){
				double beta = 0.0;
				//alpha = atan2(py,px); // alpha is the angle between x-axis and p_T vector in the xy-plane
			  //std::cout << "alpha = " << alpha << std::endl;
				beta = M_PI/2 + py/std::abs(py)*charge*M_PI/2; //beta is the angle between the x-axis and the axis perpendicular to p_T in the xy-plane
				//beta = alpha + py/std::abs(py)*charge*M_PI/2; //beta is the angle between the x-axis and the axis perpendicular to p_T in the xy-plane
			  //std::cout << "beta = " << beta << std::endl;
				cx = radius*cos(beta);//cx and cy are the x- and y-coordinates of the center of the circle that the helix performs in the xy-plane
				cy = radius*sin(beta);
}

std::vector< double > Helix::Calculate_position(){

				Calculate_radius();
				Calculate_xi();
				Calculate_circlecenter();

				//radius = radius * charge;//electrons rotate clockwise, positrons anticlockwise

				std::vector< double > position_prime;
				if (py/std::abs(py)*charge < 0 ) xi += M_PI;//if the particle's circle is on the right hand side of the p_T vector, the angle starts off from pi and not 0.
				position_prime.push_back( radius*cos(xi) + x0 + cx);
				position_prime.push_back( radius*sin(xi) + y0 + cy);

				std::vector< double > position;
				//The actual position can be gained by rotating the coordinate system with the rotation matrix:
				alpha = atan2(py,px); // alpha is the angle between x-axis and p_T vector in the xy-plane
				position.push_back( position_prime.at(0)*cos(alpha) - position_prime.at(1)*sin(alpha) );
				position.push_back( position_prime.at(0)*sin(alpha) + position_prime.at(1)*cos(alpha) );
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

std::vector< double > Helix::Get_position(std::vector< double > mom, float particle_charge, std::vector< double > origin, double pos_z){
				px = mom.at(0);
				py = mom.at(1);
				pz = mom.at(2);
				charge = particle_charge;
				x0 = origin.at(0);
				y0 = origin.at(1);
				z0 = origin.at(2);
				z = pos_z/1000.0; //to convert mm into m (the histogramm will be drawn in mm, therefore z is given in mm)
			
				return Calculate_position();

}
