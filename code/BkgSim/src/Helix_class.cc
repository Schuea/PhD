#include "Helix_class.h"

void Helix::Calculate_number_turn(){
				number_turn = (z-z0)/(pz*2*M_PI/(0.3*B));
}

void Helix::Calculate_xi(){
				xi = number_turn/(2*M_PI);
}

void Helix::Calculate_position(){

				Calculate_number_turn();
				Calculate_xi();

				position.clear();
				position.push_back( cos(xi) + x0 );
				position.push_back( sin(xi) + y0 );
				position.push_back( z + z0 );
}

std::vector< double > Helix::Get_position(float pos_z){
				z = pos_z;
				Calculate_position();

				return position;
}
