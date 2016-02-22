#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "ConfigReaderAnalysis.h"
#include "UsefulFunctions.h"
#include "CellHits_class.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "Time_class.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
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
	Time passedbytime;

	std::vector<Subdetector*> * SubDetectors = new std::vector<Subdetector*>();
	std::string* subdetector_name = new std::string("");
	SetupSubDetectorsVector(SubDetectors, subdetector_name, argument_subdetectors);
	std::string subdetectornames = (*subdetector_name);
	std::vector<float> range_array;

	range_array = SubDetectors->at(0)->GetROOTHisto_binning3D(); //SOMETHING HARD CODED!!
	float rmax = sqrt(pow(range_array[5], 2) + pow(range_array[8], 2));
	float rmin = 0.;
	float rrange = rmax - rmin;
	float zmax = 3500;//range_array[2];
	float zmin = -zmax;
	float zrange = rmax - zmin;
	std::array<float, 6> axis_vector = { float(zrange / 8.0), zmin, zmax, float(rrange / 20.0), rmin, 400};

	//Make histogram for storing the information
	std::string const title1 = "Time < 10ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm]; r [mm]";
	TH2D* histo1 = new TH2D("histo1", title1.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);
	std::string const title2 = "10ns <= Time < 20ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm]; r [mm]";
	TH2D* histo2 = new TH2D("histo2", title2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);
	std::string const title3 = "20ns <= Time < 50ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm]; r [mm]";
	TH2D* histo3 = new TH2D("histo3", title3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {

		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
			TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			tree->SetBranchStatus("*", 0);
			tree->SetBranchStatus("HitContrTime", 1);
			tree->SetBranchStatus("HitMotherVertex_x", 1);
			tree->SetBranchStatus("HitMotherVertex_y", 1);
			tree->SetBranchStatus("HitMotherVertex_z", 1);

			float actualtime = 0.0;
			double vertex_x = 0.0;
			double vertex_y = 0.0;
			double vertex_z = 0.0;

			tree->SetBranchAddress("HitContrTime", &actualtime);
			tree->SetBranchAddress("HitMotherVertex_x", &vertex_x);
			tree->SetBranchAddress("HitMotherVertex_y", &vertex_y);
			tree->SetBranchAddress("HitMotherVertex_z", &vertex_z);

			std::array<double, 3> vertex;
			float absolutetime = 0.;
			std::pair<int, int> Number_train_bunch = Set_train_bunch_number(file_iterator);
			passedbytime.Calculate_passedbytime(Number_train_bunch.first, Number_train_bunch.second);

			//Now we loop through the tree
			//Combine the two Cell ID's into a single new Cell ID
			//See how often the new Cell ID occurs in total, this is the occupancy

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				absolutetime = actualtime + passedbytime.Get_passedbytime();
				vertex = { vertex_x, vertex_y, vertex_z };
				if (absolutetime < 10.0) histo1->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
				else if (absolutetime >= 10.0 && absolutetime < 20.0) histo2->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
				else if (absolutetime >= 20.0 && absolutetime < 50.0) histo3->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
			}
			file->Close();
		}
	}

	gStyle->SetOptStat(11111);
	//Plot the histogram and save it
	TCanvas *canvas1 = new TCanvas("canvas", "canvas", 800, 600);
	histo1->Draw("colz");
	canvas1->Print("../output/hitmaps_particleorigins_time1.pdf");
	canvas1->Print("../output/hitmaps_particleorigins_time1.C");
	TCanvas *canvas2 = new TCanvas("canvas", "canvas", 800, 600);
	histo2->Draw("colz");
	canvas2->Print("../output/hitmaps_particleorigins_time2.pdf");
	canvas2->Print("../output/hitmaps_particleorigins_time2.C");
	TCanvas *canvas3 = new TCanvas("canvas", "canvas", 800, 600);
	histo3->Draw("colz");
	canvas3->Print("../output/hitmaps_particleorigins_time3.pdf");
	canvas3->Print("../output/hitmaps_particleorigins_time3.C");

	return 0;
}

