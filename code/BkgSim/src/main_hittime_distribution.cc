#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include "ConfigReaderAnalysis.h"
#include "UsefulFunctions.h"
#include "CellHits_class.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "Time_class.h"

#include "Style.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
	UsePhDStyle();

	//ConfigReaderAnalysis config(argv[1]);
	//config.setUp();
	//cout << config.getConfigName() << endl;
	//cout << config.getDoEventLooper() << endl;
	//cout << config.getMaxEvents() << endl;
	//exit(1);
	//Three occupancy plots
	//Two that are 1D histograms, y-axis is average occupancy and x-axis is radius/phi
	//The difficult plot is buffer depth plot: y-axis is probability of a specific cell occupancy occuring and x-axis is occupancy

	//The input is a TTree ROOT file(s)
	//The output is .pdf and .C files

	std::vector<std::string> *inputfilenames = new std::vector<std::string>();
	std::vector<std::string> argument_subdetectors;

	int NUMBER_OF_FILES = 0;
	bool NUMBER_OF_FILES_set = false;
	bool inputfile_set = false;
	bool subdetector_set = false;

	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-n")) {
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-s") && argv[i + 1] != std::string("-i")) {
				NUMBER_OF_FILES = std::stoi(argv[i + 1]);
				std::cout << "Number of input files = " << NUMBER_OF_FILES << std::endl;
				NUMBER_OF_FILES_set = true;
			} else {
				std::cerr << "You didn't give an argument for the number of files!" << std::endl;
			}
		}
	}
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")) {
			if (argv[i + 1] != NULL) {
				int j = 1;
				do {
					if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")) {
						inputfilenames->push_back(argv[i + j]);
						++j;
					} else {
						break;
					}
				} while (j <= NUMBER_OF_FILES);
				inputfile_set = true;
			} else {
				std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
			}
		}
		if (argv[i] == std::string("-s")) {
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-i")) {
				if (argv[i + 1] != NULL) {
					int j = 1;
					do {
						argument_subdetectors.push_back(argv[i + j]);
						++j;
					} while (argv[i + j] != NULL && argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s"));
					subdetector_set = true;
				} else {
					std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
				}
			}
		}
	}
	if (!inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set) {
		std::cerr
				<< "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
				<< std::endl;
		exit(1);
	}

	std::vector<Subdetector*> * SubDetectors = new std::vector<Subdetector*>();
	std::string* subdetector_name = new std::string("");
	SetupSubDetectorsVector(SubDetectors, subdetector_name, argument_subdetectors);
	std::string subdetectornames = (*subdetector_name);
	std::vector<float> range_array;

	//Make histogram for storing the information
	//std::vector< TH2D* > Hits_Time_rtime_2D_;
	range_array = SubDetectors->at(0)->GetROOTHisto_binning3D(); //SOMETHING HARD CODED!!
	float rmax = 0.7*sqrt(pow(range_array[5], 2) + pow(range_array[8], 2));//Make the plot a big smaller than the actual measurement of the subdetector
	float rmin = 0.;
	float rrange = rmax - rmin;
	float timemax = 0.0; //ns (one bunch spacing is 554 ns)
	if(NUMBER_OF_FILES>1) timemax = NUMBER_OF_FILES * 600.0; //ns (one bunch spacing is 554 ns)
	else timemax = 70.0; //ns (one bunch spacing is 554 ns)
	float timemin = 0.;
	float timerange = timemax - timemin;
	std::array<float, 6> axis_vector = { float(timerange / 10.0), timemin, timemax, float(rrange / 10.0), rmin, rmax };

	std::string const histo_title = "Radial position of hits over hit time for " + subdetectornames + ";Hit time [ns];r [mm]";
	std::vector<TH2D*> Hits_Time_rtime_;

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {
		Time PassedTime;

		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
			std::stringstream histoname;
			histoname << "HitTime_" << subdetectornames << "_bunch#" << file_iterator+1;
			std::string histo_name = histoname.str();
			Hits_Time_rtime_.emplace_back(new TH2D(histo_name.c_str(), histo_title.c_str(), 
									axis_vector.at(0), axis_vector.at(1), axis_vector.at(2),
									axis_vector.at(3), axis_vector.at(4), axis_vector.at(5)));
			TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			float actualtime = 0.0;
			float x_cal = 0.0;
			float y_cal = 0.0;
			double x_tracker = 0.0;
			double y_tracker = 0.0;
			tree->SetBranchStatus("*", 0);
			tree->SetBranchStatus("HitPosition_x", kTRUE); 
			tree->SetBranchStatus("HitPosition_y", kTRUE);
			
			if (tree->GetName() == std::string("Tree_EcalBarrel") 
											|| tree->GetName() == std::string("Tree_EcalEndcap")
											|| tree->GetName() == std::string("Tree_HcalBarrel") 
											|| tree->GetName() == std::string("Tree_HcalEndcap")
											|| tree->GetName() == std::string("Tree_MuonBarrel") 
											|| tree->GetName() == std::string("Tree_MuonEndcap")
											|| tree->GetName() == std::string("Tree_BeamCal") 
											|| tree->GetName() == std::string("Tree_LumiCal")) {
							tree->SetBranchStatus("HitContrTime", 1);
							tree->SetBranchAddress("HitContrTime", &actualtime);
							tree->SetBranchAddress("HitPosition_x", &x_cal);
							tree->SetBranchAddress("HitPosition_y", &y_cal);
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
							tree->SetBranchStatus("HitTime", 1);
							tree->SetBranchAddress("HitTime", &actualtime);
							tree->SetBranchAddress("HitPosition_x", &x_tracker);
							tree->SetBranchAddress("HitPosition_y", &y_tracker);
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}
			float x = 0.0;
			float y = 0.0;

			//Figure out which bunch in which train we are, in order to calculate the time that's already passed:
			std::pair<int, int> Number_train_bunch = Set_train_bunch_number(file_iterator);
			PassedTime.Calculate_passedbytime(Number_train_bunch.first, Number_train_bunch.second);//first: number of train, second: number of bunch
  		float absolutetime = 0.0;

			//Have the markers for different bunches coloured differently
			Hits_Time_rtime_.at(file_iterator)->SetMarkerColor(Number_train_bunch.second);

			//Now we loop through the tree

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				if (tree->GetName() == std::string("Tree_EcalBarrel") 
											|| tree->GetName() == std::string("Tree_EcalEndcap")
											|| tree->GetName() == std::string("Tree_HcalBarrel") 
											|| tree->GetName() == std::string("Tree_HcalEndcap")
											|| tree->GetName() == std::string("Tree_MuonBarrel") 
											|| tree->GetName() == std::string("Tree_MuonEndcap")
											|| tree->GetName() == std::string("Tree_BeamCal") 
											|| tree->GetName() == std::string("Tree_LumiCal")) {
								x = x_cal;
								y = y_cal;
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
								x = x_tracker;
								y = y_tracker;
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}
			//absolute time = time in respect to the current bunch interaction + time passed by since first bunch interaction
			  absolutetime = actualtime + PassedTime.Get_passedbytime();
				Hits_Time_rtime_.at(file_iterator)->Fill(absolutetime, sqrt(pow(x, 2) + pow(y, 2)));
			}
			file->Close();
		}
	}
	gStyle->SetOptStat(11);

	//Plot the histogram and save it
	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);

	canvas1->cd();
	//canvas1->SetLogz();
	//histo1->Draw("colz");
	Hits_Time_rtime_.at(0)->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec;
	st_vec.push_back(new TPaveStats());
	st_vec.at(0) = (TPaveStats*)Hits_Time_rtime_.at(0)->GetListOfFunctions()->FindObject("stats");
	st_vec.at(0)->SetLineColor(0+1);
	st_vec.at(0)->SetX1NDC(0.65); //new x start position
	st_vec.at(0)->SetX2NDC(0.85); //new x end position
	st_vec.at(0)->SetY1NDC(0.83); //new x start position
	st_vec.at(0)->SetY2NDC(0.9); //new x end position
	float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

	for (int file_iterator = 1; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
		Hits_Time_rtime_.at(file_iterator)->Draw("SAMES");
		canvas1->Update();
		st_vec.push_back(new TPaveStats());
		st_vec.at(file_iterator)= (TPaveStats*)Hits_Time_rtime_.at(file_iterator)->GetListOfFunctions()->FindObject("stats");
		st_vec.at(file_iterator)->SetLineColor(file_iterator+1);
		st_vec.at(file_iterator)->SetX1NDC(0.65); //new x start position
		st_vec.at(file_iterator)->SetX2NDC(0.85); //new x end position
		st_vec.at(file_iterator)->SetY2NDC(st_vec.at(file_iterator-1)->GetY1NDC()); //new x end position
		st_vec.at(file_iterator)->SetY1NDC(st_vec.at(file_iterator)->GetY2NDC()-boxsize); //new x start position
	}


	canvas1->Print(("output/hittime_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/hittime_"+subdetectornames+".cxx").c_str());

	return 0;
}

