/*
 * Time.h
 *
 *  Created on: Jan 20, 2016
 *      Author: schuea
 */

#ifndef TIMECLASS_H_
#define TIMECLASS_H_

#include <map>

class Time {
public:
	Time();
	virtual ~Time();

	void Calculate_passedbytime(int train, int bunch);
	void Set_number_of_train(int train);
	void Set_number_of_bunch(int bunch);
	void Update_time_map();

	float Get_passedbytime()  const;

private:
	float bunch_spacing; //ns
	float train_spacing; //ns = 200 ms
	int number_of_bunch;
	int number_of_train;
	float passedbytime; //ns

	std::map< std::pair< int, int >, float > passedbytimes;
};

#endif /* TIMECLASS_H_ */
