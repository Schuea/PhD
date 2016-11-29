#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"


#include <bitset>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

#include "UsefulFunctions.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "CellHits_class_new.h"
#include "Subdetector_class_new.h"

using namespace std;

long long int MakeNewCellID(double const x, double const y, Subdetector const & component);

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
	std::string argument_subdetectors;

	int NUMBER_OF_FILES = 0;
	bool NUMBER_OF_FILES_set = false;
	bool inputfile_set = false;
	bool subdetector_set = false;

	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-n")) {
			if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-s") && 
          argv[i + 1] != std::string("-i")) {
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
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-s")) {
        argument_subdetectors = argv[i + 1];
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

  Subdetector det( argument_subdetectors );

	CellHits * HitCount = new CellHits( &det );

		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
			TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, det.getName());

			//Set the branches
			tree->SetBranchStatus("*", 0);
			//tree->SetBranchStatus("HitCellID", 1);
			tree->SetBranchStatus("HitCellID0", 1);
			//tree->SetBranchStatus("HitCellID1", 1);
			tree->SetBranchStatus("HitPosition_x", 1);
			tree->SetBranchStatus("HitPosition_y", 1);
			tree->SetBranchStatus("HitPosition_z", 1);

			int HitCellID0(0);
      //int HitCellID1(0);
			float HitPosition_x(0.0), HitPosition_y(0.0), HitPosition_z(0.0);
			//double HitPosition_x(0.0), HitPosition_y(0.0);

			//tree->SetBranchAddress("HitCellID", &HitCellID0);
			tree->SetBranchAddress("HitCellID0", &HitCellID0);
			//tree->SetBranchAddress("HitCellID1", &HitCellID1);
			tree->SetBranchAddress("HitPosition_x", &HitPosition_x);
			tree->SetBranchAddress("HitPosition_y", &HitPosition_y);
			tree->SetBranchAddress("HitPosition_z", &HitPosition_z);

			//Now we loop through the tree
			//Combine the two Cell ID's into a single new Cell ID
			//See how often the new Cell ID occurs in total, this is the occupancy

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				//if (HitPosition_z < 0) continue;
				//Make a combined cell ID
        long long int HitCellID1 = MakeNewCellID(HitPosition_x,HitPosition_y,det);
				long long int const combined_cell_id = (long long) HitCellID1 << 32 | HitCellID0;
				//Use the CellHits class for storing the hit cells and their hitcounts
				HitCount->Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y, HitPosition_z);
			}
			file->Close();
		}

	//Make histogram for storing the information
  std::string subdetectorname = det.getName();
	std::string const title = "Occupancy for subdetector " + subdetectorname;
	//std::string const title = "Normalized buffer depth for subdetector " + subdetectorname;
	std::vector< TH1D* > histos;
	TH1D* All_Layers_histo = new TH1D("All layers", title.c_str(), 10, 0, 10);
	std::vector< TPaveStats* > stats;

	int tot_no_hits = 0;

	int max_num_layers = det.getNumberOfLayers();
	for (int number_layer = 0; number_layer <= max_num_layers; ++number_layer) {
		std::stringstream layername;
		layername << "Layer " << number_layer;
		histos.emplace_back(new TH1D(layername.str().c_str(), title.c_str(), 100000, 0, 100000));
	}
  for (size_t vecpos = 0; vecpos < HitCount->Get_HitCount().size(); ++vecpos) {
    if(HitCount->Get_HitCount().at(vecpos) > 0){
      histos.at(HitCount->Get_Layer().at(vecpos))->Fill(HitCount->Get_HitCount().at(vecpos));
      All_Layers_histo->Fill(HitCount->Get_HitCount().at(vecpos));
      tot_no_hits += HitCount->Get_HitCount().at(vecpos);
    }
  }

	std::cout<< "---------------" <<std::endl;
	std::cout<< "Total number of hits counted for this histogram: " << tot_no_hits <<std::endl;
	std::cout<< "---------------" <<std::endl;

	//NormalizeHistogram(histo, 1.0);
	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	canvas->SetLogy(1);
	float boxsize = 0.0;
	for (size_t number_histo = 0; number_histo< histos.size(); ++number_histo) {
					histos.at(number_histo)->SetMinimum(0.1);
					if(number_histo == 0){
									histos.at(number_histo)->SetLineColor(2);
									histos.at(number_histo)->Draw();
									canvas->Update();
									TPaveStats* st =  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats");
									st->SetX1NDC(0.65); //new x start position
									st->SetX2NDC(0.85); //new x end position
									st->SetY1NDC(0.8); //new y start position
									st->SetY2NDC(0.9); //new y end position
									stats.push_back(st);
									boxsize = stats.at(number_histo)->GetY2NDC() - stats.at(number_histo)->GetY1NDC();
					}
					if(number_histo > 0){
									histos.at(number_histo)->SetLineColor(2+number_histo);
									histos.at(number_histo)->Draw("SAMES");
									canvas->Update();
									TPaveStats* st =  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats");
									stats.push_back(st);
									stats.at(number_histo)->SetTextColor(2+number_histo);
									stats.at(number_histo)->SetX1NDC(0.65); //new x start position
									stats.at(number_histo)->SetX2NDC(0.85); //new x end position
									stats.at(number_histo)->SetY2NDC(stats.at(number_histo-1)->GetY1NDC()); //new y end position
									stats.at(number_histo)->SetY1NDC(stats.at(number_histo-1)->GetY2NDC()-boxsize); //new y start position
					}
	}
  std::stringstream output;
  output << "output/muon_occupancy_" << subdetectorname;
	canvas->Print((output.str() + ".pdf").c_str());
	canvas->Print((output.str() + ".cxx").c_str());

	All_Layers_histo->SetMinimum(0.1);
	All_Layers_histo->SetLineColor(2);
	All_Layers_histo->Draw();
	canvas->Update();
	TPaveStats* st =  (TPaveStats*)All_Layers_histo->GetListOfFunctions()->FindObject("stats");
	st->SetX1NDC(0.65); //new x start position
	st->SetX2NDC(0.85); //new x end position
	st->SetY1NDC(0.8); //new y start position
	st->SetY2NDC(0.9); //new y end position

  std::stringstream output2;
  output2 << "output/muon_occupancy_all_layers_" << subdetectorname;
	canvas->Print((output2.str() + ".pdf").c_str());
	canvas->Print((output2.str() + ".cxx").c_str());

	return 0;
}

long long int MakeNewCellID(double const x, double const y, Subdetector const & component){
  int newX = static_cast<int>(x/component.getCellSizeArea()); //Check if Cell Size Area is the same as Cell Dimension
  int newY = static_cast<int>(y/component.getCellSizeArea()); //Check if Cell Size Area is the same as Cell Dimension
  if(x >= 0) ++newX;
  if(y >= 0) ++newY;
  bitset<32> bitY = newY;
  newY = 0;
  for(int i = 0; i < 31; ++i){
    newY += bitY[i]*pow(2,i);
  }
  return newX << 32 | newY;
}
