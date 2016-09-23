#include "Helix_class.h"

void Helix::Calculate_radius(){
				radius = sqrt(px*px+py*py)/(0.3*B);
}

void Helix::Calculate_number_turn(){
				number_turn = (z-z0)/(pz*2*M_PI/(0.3*B));
}

void Helix::Calculate_xi(){
				Calculate_number_turn();
				xi = number_turn/(2*M_PI);
}

void Helix::Calculate_position(){

				Calculate_radius();
				radius = radius * charge;
				Calculate_xi();

				position.clear();
				position.push_back( radius*cos(xi) + x0 );
				position.push_back( radius*sin(xi) + y0 );
				position.push_back( z + z0 );
}

std::vector< double > Helix::Get_position(float pos_z){
				z = pos_z;
				Calculate_position();

				return position;
}
