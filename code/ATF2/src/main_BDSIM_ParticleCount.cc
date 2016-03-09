#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TROOT.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "UsefulFunctions.h"
#include "Style.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
	UsePhDStyle();
	//The input is a TTree ROOT file(s)
	//The output is .pdf and .C files

	std::vector < std::string > *inputfilenames =
			new std::vector<std::string>();
	std::string sample;

	int NUMBER_OF_FILES = 0;
	bool NUMBER_OF_FILES_set = false;
	bool inputfile_set = false;
	bool sample_set = false;

	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-n")) {
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-s")
					&& argv[i + 1] != std::string("-i")) {
				NUMBER_OF_FILES = std::stoi(argv[i + 1]);
				std::cout << "Number of input files = " << NUMBER_OF_FILES
						<< std::endl;
				NUMBER_OF_FILES_set = true;
			} else {
				std::cerr
						<< "You didn't give an argument for the number of files!"
						<< std::endl;
			}
		}
	}
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i") && argv[i + 1] != std::string("-n")
				&& argv[i + 1] != std::string("-s")) {
			if (argv[i + 1] != NULL) {
				int j = 1;
				do {
					if (argv[i + j] != std::string("-n")
							&& argv[i + j] != std::string("-s")) {
						inputfilenames->push_back(argv[i + j]);
						++j;
					} else {
						break;
					}
				} while (j <= NUMBER_OF_FILES);
				inputfile_set = true;
			} else {
				std::cerr << "You didn't give an argument for the inputfile(s)!"
						<< std::endl;
			}
		}
		if (argv[i] == std::string("-s")) {
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n")
					&& argv[i + 1] != std::string("-i")) {
				sample = argv[i + 1];
				sample_set = true;
			} else {
				std::cerr << "You didn't give an argument for the sample!"
						<< std::endl;
			}
		}
	}
	if (!inputfile_set || !sample_set || !NUMBER_OF_FILES_set) {
		std::cerr
				<< "You didn't give the name for the sample, the inputfiles or the amount of files. Please try again!"
				<< std::endl;
		exit(1);
	}

	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	//Make histogram for storing the information
	std::string const title =
			"Number of primary and secondary particles;PDG code;Count";
	TH1D *histo_primaries = new TH1D("Primaries", title.c_str(), 70, -35, 35);
	histo_primaries->SetLineColor(4);	//blue
	TH1D *histo_secondaries = new TH1D("Secondaries", title.c_str(), 70, -35,
			35);
	histo_secondaries->SetLineColor(2);	//red

	std::string Sampler_Tree_name = "Sampler_" + sample;
	std::cout << Sampler_Tree_name << std::endl;
	
	//For setting the y-axis of the plot
	int primary_number = 0;
	int secondary_number = 0;

	for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES;
			++file_iterator) {
		TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
		TTree *tree = nullptr;
		file->GetObject(Sampler_Tree_name.c_str(), tree);

		//Set the branches
		tree->SetBranchStatus("*", 0);
		tree->SetBranchStatus("partID", 1);
		tree->SetBranchStatus("parentID", 1);
		/*
		 tree->SetBranchStatus("x0", 1);
		 tree->SetBranchStatus("y0", 1);
		 tree->SetBranchStatus("z0", 1);
		 tree->SetBranchStatus("x_prod", 1);
		 tree->SetBranchStatus("y_prod", 1);
		 tree->SetBranchStatus("z_prod", 1);
		 */

		int partID(0), parentID(0);
		//float x0(0.0), y0(0.0), z0(0.0);
		//float x_prod(0.0), y_prod(0.0), z_prod(0.0);

		tree->SetBranchAddress("partID", &partID);
		tree->SetBranchAddress("parentID", &parentID);
		/*
		 tree->SetBranchAddress("x0", &x0);
		 tree->SetBranchAddress("y0", &y0);
		 tree->SetBranchAddress("z0", &z0);
		 tree->SetBranchAddress("x_prod", &x_prod);
		 tree->SetBranchAddress("y_prod", &y_prod);
		 tree->SetBranchAddress("z_prod", &z_prod);
		 */

		//Now we loop through the tree
		long long int const entries = tree->GetEntries();
		for (long long int i = 0; i < entries; ++i) {
			tree->GetEntry(i);
			if (parentID == 0) {
				histo_primaries->Fill(partID);
				primary_number++;
			} else {
				histo_secondaries->Fill(partID);
				secondary_number++;
			}
		} 
		file->Close();
	}
	
	int maxnumber = secondary_number;
	if (primary_number > secondary_number){
		maxnumber = primary_number;
	}
	histo_primaries->SetMaximum(maxnumber + 0.5*maxnumber);
	histo_secondaries->SetMaximum(maxnumber + 0.5*maxnumber);
	histo_primaries->SetMinimum(0.01);
	histo_secondaries->SetMinimum(0.01);

	//NormalizeHistogram(histo, 1.0);
	canvas->SetLogy();
	
	histo_primaries->Draw();
	histo_secondaries->Draw("SAMES");

	canvas->Update();
	TPaveStats *st1 =
			(TPaveStats*) histo_primaries->GetListOfFunctions()->FindObject(
					"stats");
	st1->SetTextColor(4);
	st1->SetX1NDC(0.75); //new x start position
	st1->SetX2NDC(0.9); //new x end position
	st1->SetY1NDC(0.7); //new x start position
	st1->SetY2NDC(0.9); //new x end position
	float boxsize = st1->GetY2NDC() - st1->GetY1NDC();
	TPaveStats *st2 =
			(TPaveStats*) histo_secondaries->GetListOfFunctions()->FindObject(
					"stats");
	st2->SetTextColor(2);
	st2->SetX1NDC(0.75); //new x start position
	st2->SetX2NDC(0.9); //new x end position
	st2->SetY2NDC(st1->GetY1NDC()); //new x start position
	st2->SetY1NDC(st2->GetY2NDC() - boxsize); //new x end position

	canvas->Print("output/Histo_PrimariesSecondaries_Halo.pdf");
	canvas->Print("output/Histo_PrimariesSecondaries_Halo.cxx");

	return 0;
}
