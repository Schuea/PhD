#include "TFile.h"
#include "TTree.h"

#include <sstream>
#include <iostream>

#include "Subdetector_class.h"


void InitializeAllCaloSubdetectors(std::vector<Subdetector*> * SubDetectors) {
	SubDetectors->push_back(new EcalBarrel());
	SubDetectors->push_back(new EcalEndcap());
	SubDetectors->push_back(new HcalBarrel());
	SubDetectors->push_back(new HcalEndcap());
	SubDetectors->push_back(new MuonBarrel());
	SubDetectors->push_back(new MuonEndcap());
	SubDetectors->push_back(new BeamCal());
	SubDetectors->push_back(new LumiCal());
}
void InitializeAllTrackerSubdetectors(std::vector<Subdetector*> * SubDetectors) {
	SubDetectors->push_back(new SiVertexBarrel());
	SubDetectors->push_back(new SiVertexEndcap());
	SubDetectors->push_back(new SiTrackerBarrel());
	SubDetectors->push_back(new SiTrackerEndcap());
	SubDetectors->push_back(new SiTrackerForward());
}
void InitializeWhichSubdetector(std::string SubdetectorName, std::vector<Subdetector*> * SubDetectors) {
	std::cout << "SubdetectorName = " << SubdetectorName << std::endl;
	if (SubdetectorName == std::string("EcalBarrel"))
		SubDetectors->push_back(new EcalBarrel());
	else if (SubdetectorName == std::string("EcalEndcap"))
		SubDetectors->push_back(new EcalEndcap());
	else if (SubdetectorName == std::string("HcalBarrel"))
		SubDetectors->push_back(new HcalBarrel());
	else if (SubdetectorName == std::string("HcalEndcap"))
		SubDetectors->push_back(new HcalEndcap());
	else if (SubdetectorName == std::string("MuonBarrel"))
		SubDetectors->push_back(new MuonBarrel());
	else if (SubdetectorName == std::string("MuonEndcap"))
		SubDetectors->push_back(new MuonEndcap());
	else if (SubdetectorName == std::string("BeamCal"))
		SubDetectors->push_back(new BeamCal());
	else if (SubdetectorName == std::string("LumiCal"))
		SubDetectors->push_back(new LumiCal());
	else if (SubdetectorName == std::string("SiVertexBarrel"))
		SubDetectors->push_back(new SiVertexBarrel());
	else if (SubdetectorName == std::string("SiVertexEndcap"))
		SubDetectors->push_back(new SiVertexEndcap());
	else if (SubdetectorName == std::string("SiTrackerBarrel"))
		SubDetectors->push_back(new SiTrackerBarrel());
	else if (SubdetectorName == std::string("SiTrackerEndcap"))
		SubDetectors->push_back(new SiTrackerEndcap());
	else if (SubdetectorName == std::string("SiTrackerForward"))
		SubDetectors->push_back(new SiTrackerForward());
	else {
		std::cerr << "The given name doesn't match any subdetector!" << std::endl;
		std::terminate();
	}
}
void InitializeAllSubdetectors(std::vector<Subdetector*> * SubDetectors) {
	InitializeAllCaloSubdetectors(SubDetectors);
	InitializeAllTrackerSubdetectors(SubDetectors);
}

void SetupSubDetectorsVector(std::vector<Subdetector*> * SubDetectors, std::string *several_subdetector_names,
		std::vector<std::string> argument_subdetectors) {
	for (int s = 0; s < argument_subdetectors.size(); ++s) {

		if (argument_subdetectors.size() == 1) {
			(*several_subdetector_names) = argument_subdetectors.at(s);
			if (argument_subdetectors.at(s) == std::string("allCalo")) {

				SubDetectors->clear();
				InitializeAllCaloSubdetectors(SubDetectors);
			}

			else if (argument_subdetectors.at(s) == std::string("allTracker")) {
				(*several_subdetector_names) = argument_subdetectors.at(s);

				SubDetectors->clear();
				InitializeAllTrackerSubdetectors(SubDetectors);
			}

			else if (argument_subdetectors.at(s) == std::string("all")) {
				(*several_subdetector_names) = argument_subdetectors.at(s);

				SubDetectors->clear();
				InitializeAllSubdetectors(SubDetectors);
			} else {
				InitializeWhichSubdetector(argument_subdetectors.at(s), SubDetectors);
			}
		} else {
			if (argument_subdetectors.at(s) == std::string("all")) {
				(*several_subdetector_names) = argument_subdetectors.at(s);

				SubDetectors->clear();
				InitializeAllSubdetectors(SubDetectors);
				break;
			} else if (argument_subdetectors.at(s) == std::string("allCalo")) {
				(*several_subdetector_names) = argument_subdetectors.at(s);

				SubDetectors->clear();
				InitializeAllCaloSubdetectors(SubDetectors);

				if (argument_subdetectors.at(s + 1) == std::string("allTracker")) {
					std::cerr
							<< "If you want to compare all subdetectors, just use the argument 'all' after the flag -s!"
							<< std::endl;
					std::terminate();
				} else
					break;
			} else if (argument_subdetectors.at(s) == std::string("allTracker")) {
				(*several_subdetector_names) = argument_subdetectors.at(s);

				SubDetectors->clear();
				InitializeAllTrackerSubdetectors(SubDetectors);

				if (argument_subdetectors.at(s + 1) == std::string("allCalo")) {
					std::cerr
							<< "If you want to compare all subdetectors, just use the argument 'all' after the flag -s!"
							<< std::endl;
					std::terminate();
				} else
					break;

			} else {
				(*several_subdetector_names) += argument_subdetectors.at(s);
				InitializeWhichSubdetector(argument_subdetectors.at(s), SubDetectors);
			}
		}
	}
	for (size_t i = 0; i < SubDetectors->size(); ++i) {
		SubDetectors->at(i)->Initialize();
	}
}


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
