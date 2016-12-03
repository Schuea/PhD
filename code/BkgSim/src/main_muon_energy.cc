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
#include "GeneralFunctions_SiDBkgSim_new.h"

using namespace std;

void Draw_All_plots_together ( std::vector< TH1D* > histo, TCanvas* canvas);

double weight = 0.08846; //The weight is the same for both scenarios, for both, the electron and the positron line, e.g.: 5sp+wall, elec: (4321/10155 * 898.34)/4321

int main(int const argc, char const * const * const argv) {
	std::vector<std::string> *inputfilenames = new std::vector<std::string>();

	int NUMBER_OF_FILES = 0;
	bool NUMBER_OF_FILES_set = false;
	bool inputfile_set = false;

	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-n")) {
			if (argv[i + 1] != NULL && 
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
    if (argv[i] == std::string("-i") && 
        argv[i + 1] != std::string("-n") ) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-n") ) {
            inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_FILES*2);//2 files per scenario
        inputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
      }
    }
  }
	if (!inputfile_set || !NUMBER_OF_FILES_set) {
		std::cerr
				<< "You didn't give the name for the inputfiles or the amount of files. Please try again!"
				<< std::endl;
		exit(1);
	}

  std::vector< TH1D* > Energy_histos;
  float max_energy = 260.0;
  std::string name1 = "Energy_sp";
  std::string name2 = "Energy_spwall";
  std::string title = "Energy distribution of the muons reaching the detector;Energy (GeV);Number of muons";
  Energy_histos.emplace_back(new TH1D(name1.c_str(),title.c_str(),(int)max_energy/3.0,0,(int)max_energy) );
  Energy_histos.emplace_back(new TH1D(name2.c_str(),title.c_str(),(int)max_energy/3.0,0,(int)max_energy) );

  for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES*2; ++file_iterator) {
    std::cout << inputfilenames->at(file_iterator) << std::endl;
    TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, "MCP");

    //Set the branches
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("Energy", 1);

    double energy(0.0);

    tree->SetBranchAddress("Energy", &energy);

    //Now we loop through the tree

    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      Energy_histos.at(file_iterator/NUMBER_OF_FILES)->Fill(energy, weight);
    }
    file->Close();
  }

	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	canvas->SetLogy(1);

   Draw_All_plots_together( Energy_histos, canvas); 
  std::stringstream output;
  output << "output/muon_energy";
  canvas->Print((output.str() + ".pdf").c_str());
  canvas->Print((output.str() + ".cxx").c_str());
  
 	return 0;
}
void Draw_All_plots_together ( std::vector< TH1D* > histo, TCanvas* canvas){
 int tot = 0;
  for(int bin = 1; bin <= histo.at(0)->GetNbinsX(); ++bin){
    tot += histo.at(0)->GetBinContent(bin);
  }
  int tot2 = 0;
  for(int bin = 1; bin <= histo.at(1)->GetNbinsX(); ++bin){
    tot2 += histo.at(1)->GetBinContent(bin);
  }
  for (size_t vec_entry = 0; vec_entry < histo.size(); ++vec_entry){
			histo.at(vec_entry)->Sumw2(1);
  }
	double max=GetMinMaxForMultipleOverlappingHistograms(histo,true).second;
  for (size_t vec_entry = 0; vec_entry < histo.size(); ++vec_entry){
    histo.at(vec_entry)->SetMinimum(max);
    histo.at(vec_entry)->SetMinimum(0.1);
    if(vec_entry == 0){
      histo.at(vec_entry)->SetLineColor(2);
      histo.at(vec_entry)->Draw("hist");
      canvas->Update();
    }
    else{
      histo.at(vec_entry)->SetLineColor(4);
      histo.at(vec_entry)->Draw("hist,SAMES");
      canvas->Update();
    }
  }
  TPaveText *text1 = new TPaveText(0.75,0.8,0.95,0.9,"brNDC");
  text1->SetTextFont(62);
  text1->SetTextColor(2);
  text1->SetFillColor(0);
  text1->AddText("5 spoilers");
  text1->AddLine(0,0.5,1,0.5);
  std::stringstream entries_sp;
  entries_sp << "Entries = " << tot;
  text1->AddText(entries_sp.str().c_str());
  TPaveText *text2 = new TPaveText(0.75,0.7,0.95,0.8,"brNDC");
  text2->SetTextFont(62);
  text2->SetTextColor(4);
  text2->SetFillColor(0);
  //text->AddLine(0,0.5,1,0.5);
  text2->AddText("5 spoilers + wall");
  text2->AddLine(0,0.5,1,0.5);
  std::stringstream entries_spwall;
  entries_spwall << "Entries = " << tot2;
  text2->AddText(entries_spwall.str().c_str());
  text1->Draw();
  text2->Draw();

}

