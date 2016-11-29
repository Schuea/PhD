#ifndef GENERALFUNCTIONSFORSIDBKGSIM_NEW_H_
#define GENERALFUNCTIONSFORSIDBKGSIM_NEW_H_

#include "Subdetector_class_new.h"
#include "TFile.h"
#include "TTree.h"
#include <sstream>
#include <iostream>

TTree* Get_TTree(TFile* inputfile, std::string subdetector_name);

std::pair<int, int> Set_train_bunch_number(int number_of_file);

#endif //GENERALFUNCTIONSFORSIDBKGSIM_H_
