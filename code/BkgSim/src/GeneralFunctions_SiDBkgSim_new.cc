#include "GeneralFunctions_SiDBkgSim_new.h"


TTree* Get_TTree(TFile* inputfile, std::string subdetector_name) {
	std::stringstream temp;
	temp << "Tree_" << subdetector_name;

	TTree* Tree = nullptr;
	inputfile->GetObject(temp.str().c_str(), Tree);
	if (!Tree) {
		throw std::exception();
	}
	return Tree;
}

std::pair<int, int> Set_train_bunch_number(int number_of_file) {
  int number_of_train = 0;
  int number_of_bunch = 0;
  if ((number_of_file + 1) <= 1312) { //number of bunches starts with 1, not 0
    number_of_train = 1;
    number_of_bunch = number_of_file + 1;
  } else if ((number_of_file + 1) > 1312 && (number_of_file + 1) <= 2624) {
    number_of_train = 2;
    number_of_bunch = (number_of_file + 1) - 1312;
  }
  return std::pair<int, int>(number_of_train, number_of_bunch);
}
