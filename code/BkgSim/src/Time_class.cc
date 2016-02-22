/*
 * Time_class.cpp
 *
 *  Created on: Jan 20, 2016
 *      Author: schuea
 */

#include "Time_class.h"
#include <iostream>
#include <map>

Time::Time(): bunch_spacing(554.0), train_spacing(200000000.0),
	number_of_bunch(0), number_of_train(0), passedbytime(0.0) {
}

Time::~Time() {
	// TODO Auto-generated destructor stub
}

void Time::Set_number_of_train(int train){
	number_of_train = train;
}
void Time::Set_number_of_bunch(int bunch){
	number_of_bunch = bunch;
}
void Time::Update_time_map(){
	passedbytimes[std::pair< int, int >(number_of_train, number_of_bunch)] = passedbytime;
}
void Time::Calculate_passedbytime(int train, int bunch){
	Set_number_of_train(train);
	Set_number_of_bunch(bunch);

	if (number_of_train == 1 && number_of_bunch == 1) passedbytime = 0.0;
	else if (number_of_train == 1 && number_of_bunch > 1) passedbytime = passedbytimes[std::pair< int, int >(number_of_train, number_of_bunch-1)] + bunch_spacing;
	else if (number_of_train > 1 && number_of_bunch == 1) passedbytime = passedbytimes[std::pair< int, int >(number_of_train-1, 1312)] + train_spacing;
	else if (number_of_train > 1 && number_of_bunch > 1) passedbytime = passedbytimes[std::pair< int, int >(number_of_train-1, number_of_bunch-1)] + bunch_spacing + train_spacing;

	Update_time_map();
}

float Time::Get_passedbytime() const{
	return passedbytime;
}
