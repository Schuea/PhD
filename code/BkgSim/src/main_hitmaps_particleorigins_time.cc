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

	range_array = SubDetectors->at(0)->GetROOTHisto_binning3D(); //SOMETHING HARD CODED!!
	float rmax = 1.5*sqrt(pow(range_array[5], 2) + pow(range_array[8], 2));//Make the plot a big bigger than the data
	float rmin = 0.;
	float rrange = rmax - rmin;
	float zmax = 3500;
	//float zmax = range_array[2];
	float zmin = -zmax;
	float zrange = rmax - zmin;
	std::array<float, 6> axis_vector = { float(zrange / 8.0), zmin, zmax, float(rrange / 10.0), rmin, rmax};

	//Make histogram for storing the information
	std::string const title1 = "Time < 10ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm];r [mm];# of origins";
	TH2D* histo1 = new TH2D("histo1", title1.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);
	std::string const title2 = "10ns <= Time < 20ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm];r [mm];# of origins";
	TH2D* histo2 = new TH2D("histo2", title2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);
	std::string const title3 = "20ns <= Time < 30ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm];r [mm];# of origins";
	TH2D* histo3 = new TH2D("histo3", title3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);
	std::string const title4 = "30ns <= Time < 50ns, Hit maps of particle origins for those particles hitting subdetector "
			+ subdetectornames + ";z [mm];r [mm];# of origins";
	TH2D* histo4 = new TH2D("histo4", title4.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
			axis_vector[4], axis_vector[5]);

	std::stringstream subdetector_names;

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {
	  subdetector_names << SubDetectors->at(subdetector_iterator)->GetName();

		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
			TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			float actualtime = 0.0;
			double vertex_x = 0.0;
			double vertex_y = 0.0;
			double vertex_z = 0.0;
			tree->SetBranchStatus("*", 0);
			
			if (tree->GetName() == std::string("Tree_EcalBarrel") 
											|| tree->GetName() == std::string("Tree_EcalEndcap")
											|| tree->GetName() == std::string("Tree_HcalBarrel") 
											|| tree->GetName() == std::string("Tree_HcalEndcap")
											|| tree->GetName() == std::string("Tree_MuonBarrel") 
											|| tree->GetName() == std::string("Tree_MuonEndcap")
											|| tree->GetName() == std::string("Tree_BeamCal") 
											|| tree->GetName() == std::string("Tree_LumiCal")) {
							tree->SetBranchStatus("HitContrTime", 1);
							tree->SetBranchStatus("HitMotherVertex_x", 1);
							tree->SetBranchStatus("HitMotherVertex_y", 1);
							tree->SetBranchStatus("HitMotherVertex_z", 1);
							tree->SetBranchAddress("HitContrTime", &actualtime);
							tree->SetBranchAddress("HitMotherVertex_x", &vertex_x);
							tree->SetBranchAddress("HitMotherVertex_y", &vertex_y);
							tree->SetBranchAddress("HitMotherVertex_z", &vertex_z);
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
							tree->SetBranchStatus("HitTime", 1);
							tree->SetBranchStatus("HitParticleVertex_x", 1);
							tree->SetBranchStatus("HitParticleVertex_y", 1);
							tree->SetBranchStatus("HitParticleVertex_z", 1);
							tree->SetBranchAddress("HitTime", &actualtime);
							tree->SetBranchAddress("HitParticleVertex_x", &vertex_x);
							tree->SetBranchAddress("HitParticleVertex_y", &vertex_y);
							tree->SetBranchAddress("HitParticleVertex_z", &vertex_z);
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}

			std::array<double, 3> vertex;
			std::pair<int, int> Number_train_bunch = Set_train_bunch_number(file_iterator);

			//Now we loop through the tree
			//Combine the two Cell ID's into a single new Cell ID
			//See how often the new Cell ID occurs in total, this is the occupancy

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				vertex = { vertex_x, vertex_y, vertex_z };
				if (actualtime < 10.0) histo1->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
				else if (actualtime >= 10.0 && actualtime < 20.0) histo2->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
				else if (actualtime >= 20.0 && actualtime < 30.0) histo3->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
				else if (actualtime >= 30.0 && actualtime < 50.0) histo4->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
			}
			file->Close();
		}
	}
	gStyle->SetOptStat(111111);

	//Plot the histogram and save it
	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);
	TCanvas *canvas2 = new TCanvas("canvas2", "canvas", 800, 600);
	TCanvas *canvas3 = new TCanvas("canvas3", "canvas", 800, 600);
	TCanvas *canvas4 = new TCanvas("canvas4", "canvas", 800, 600);

	canvas1->cd();
	canvas1->SetLogz();
	histo1->Draw("colz");
	canvas1->Update();
	TPaveStats *st1 = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
	st1->SetX1NDC(0.65); //new x start position
	st1->SetX2NDC(0.85); //new x end position
	st1->SetY1NDC(0.6); //new x start position
	st1->SetY2NDC(0.9); //new x end position


	canvas1->Print(("output/hitmaps_particleorigins_time1_"+subdetector_names.str()+".pdf").c_str());
	canvas1->Print(("output/hitmaps_particleorigins_time1_"+subdetector_names.str()+".cxx").c_str());

	canvas2->cd();
	canvas2->SetLogz();
	histo2->Draw("colz");
	canvas2->Update();
	TPaveStats *st2 = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
	st2->SetX1NDC(0.65); //new x start position
	st2->SetX2NDC(0.85); //new x end position
	st2->SetY1NDC(0.6); //new x start position
	st2->SetY2NDC(0.9); //new x end position

	canvas2->Print(("output/hitmaps_particleorigins_time2_"+subdetector_names.str()+".pdf").c_str());
	canvas2->Print(("output/hitmaps_particleorigins_time2_"+subdetector_names.str()+".cxx").c_str());

	canvas3->cd();
	canvas3->SetLogz();
	histo3->Draw("colz");
	canvas3->Update();
	TPaveStats *st3 = (TPaveStats*)histo3->GetListOfFunctions()->FindObject("stats");
	st3->SetX1NDC(0.65); //new x start position
	st3->SetX2NDC(0.85); //new x end position
	st3->SetY1NDC(0.6); //new x start position
	st3->SetY2NDC(0.9); //new x end position

	canvas3->Print(("output/hitmaps_particleorigins_time3_"+subdetector_names.str()+".pdf").c_str());
	canvas3->Print(("output/hitmaps_particleorigins_time3_"+subdetector_names.str()+".cxx").c_str());

	canvas4->cd();
	canvas4->SetLogz();
	histo4->Draw("colz");
	canvas4->Update();
	TPaveStats *st4 = (TPaveStats*)histo4->GetListOfFunctions()->FindObject("stats");
	st4->SetX1NDC(0.65); //new x start position
	st4->SetX2NDC(0.85); //new x end position
	st4->SetY1NDC(0.6); //new x start position
	st4->SetY2NDC(0.9); //new x end position

	canvas4->Print(("output/hitmaps_particleorigins_time4_"+subdetector_names.str()+".pdf").c_str());
	canvas4->Print(("output/hitmaps_particleorigins_time4_"+subdetector_names.str()+".cxx").c_str());


	return 0;
}

