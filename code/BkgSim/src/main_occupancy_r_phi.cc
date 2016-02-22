#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>

#include "ConfigReaderAnalysis.h"
#include "UsefulFunctions.h"
#include "CellHits_class.h"

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

	//Open the tree for making the plot
	std::vector<std::string> *inputfilenames = new std::vector<std::string>();
	std::string subdetector = "";
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
				subdetector = argv[i + 1];
				std::cout << "Subdetectors = " << subdetector << std::endl;
				subdetector_set = true;

			} else {
				std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
			}
		}
	}
	if (!inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set) {
		std::cerr
				<< "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
				<< std::endl;
		exit(1);
	}

	string const tree_name = "Tree_" + subdetector;
	//map<string, int> hit_map;
	CellHits HitCount;

	for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
		TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
		TTree *tree;
		file->GetObject(tree_name.c_str(), tree);

		//Set the branches
		tree->SetBranchStatus("*", 0);
		tree->SetBranchStatus("HitCellID0", 1);
		tree->SetBranchStatus("HitCellID1", 1);
		tree->SetBranchStatus("HitPosition_x", 1);
		tree->SetBranchStatus("HitPosition_y", 1);

		int HitCellID0(0), HitCellID1(0);
		float HitPosition_x(0.0), HitPosition_y(0.0);

		tree->SetBranchAddress("HitCellID0", &HitCellID0);
		tree->SetBranchAddress("HitCellID1", &HitCellID1);
		tree->SetBranchAddress("HitPosition_x", &HitPosition_x);
		tree->SetBranchAddress("HitPosition_y", &HitPosition_y);

		//Now we loop through the tree
		//Combine the two Cell ID's into a single new Cell ID
		//See how often the new Cell ID occurs in total, this is the occupancy

		long long int const entries = tree->GetEntries();
		for (long long int i = 0; i < entries; ++i) {
			tree->GetEntry(i);
			//Make a combined cell ID
			long long int const combined_cell_id = (long long) HitCellID1 << 32 | HitCellID0;
			//Use the CellHits class for storing the hit cells and their hitcounts
			HitCount.Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y);
		}
		file->Close();
	}

	//Make histogram for storing the information
	std::string const title = "Normalized buffer depth for subdetector " + subdetector;
	TH1D *histo = new TH1D("histo", title.c_str(), 20, 0, 20);

	for (size_t hitcounts = 0; hitcounts < HitCount.Get_HitCount().size(); ++hitcounts) {
		histo->Fill(HitCount.Get_HitCount().at(hitcounts));
	}

	NormalizeHistogram(histo, 1.0);

	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	canvas->SetLogy(1);
	histo->SetMinimum(0.0000001);
	histo->Draw();
	canvas->Print("../output/buffer_depth.pdf");
	canvas->Print("../output/buffer_depth.C");

	return 0;
}

