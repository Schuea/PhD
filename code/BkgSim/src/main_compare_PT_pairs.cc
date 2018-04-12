#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include "UsefulFunctions.h"
#include "GeneralFunctions_SiDBkgSim_new.h"
#include "Style.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();
				TH1::SetDefaultSumw2();

				std::vector<std::string> *inputfilenames_250 = new std::vector<std::string>();
				std::vector<std::string> *inputfilenames_500 = new std::vector<std::string>();

				int NUMBER_OF_FILES_250 = 0;
				bool NUMBER_OF_FILES_250_set = false;
				int NUMBER_OF_FILES_500 = 0;
				bool NUMBER_OF_FILES_500_set = false;

				bool inputfiles_250_set = false;
				bool inputfiles_500_set = false;

        for (int i = 1; i < argc; i++) {
          if (argv[i] == std::string("-n250")) {
            if (argv[i + 1] != NULL && argv[i + 1] != std::string("-i250") 
                && argv[i + 1] != std::string("-i500")
                && argv[i + 1] != std::string("-n500")) {
              NUMBER_OF_FILES_250 = std::stoi(argv[i + 1]);
              std::cout << "Number of pair input files = " << NUMBER_OF_FILES_250 << std::endl;
              NUMBER_OF_FILES_250_set = true;
            } else {
              std::cerr << "You didn't give an argument for the number of pair files!" << std::endl;
            }
          }
          if (argv[i] == std::string("-n500")) {
            if (argv[i + 1] != NULL && argv[i + 1] != std::string("-i250") 
                && argv[i + 1] != std::string("-i500")
                && argv[i + 1] != std::string("-n250")) {
              NUMBER_OF_FILES_500 = std::stoi(argv[i + 1]);
              std::cout << "Number of pair input files = " << NUMBER_OF_FILES_500 << std::endl;
              NUMBER_OF_FILES_500_set = true;
            } else {
              std::cerr << "You didn't give an argument for the number of pair files!" << std::endl;
            }
          }
        }  
        for (int i = 1; i < argc; i++) {
          if (argv[i] == std::string("-i250") 
              && argv[i + 1] != std::string("-n250") 
              && argv[i + 1] != std::string("-i500")
              && argv[i + 1] != std::string("-n500")) {
            if (argv[i + 1] != NULL) {
              std::string filelist = argv[i + 1];
              if (access(filelist.c_str(), F_OK) == -1 ) {
                std::cerr << "The text file does not exist!" << std::endl;
                exit(2);
              }
              else {
                std::ifstream inputfilelist( filelist );
                int line = 1;
                std::string file;
                do {
                  std::getline(inputfilelist, file);
                  inputfilenames_250->push_back(file);
                  ++line;
                } while (line <= NUMBER_OF_FILES_250);
                inputfiles_250_set = true;
              }

            } else {
              std::cerr << "You didn't give an argument for the pair inputfile(s)!" << std::endl;
            }
          }
          if (argv[i] == std::string("-i500") 
              && argv[i + 1] != std::string("-n250") 
              && argv[i + 1] != std::string("-i250")
              && argv[i + 1] != std::string("-n500")) {
            if (argv[i + 1] != NULL) {
              std::string filelist = argv[i + 1];
              if (access(filelist.c_str(), F_OK) == -1 ) {
                std::cerr << "The text file does not exist!" << std::endl;
                exit(2);
              }
              else {
                std::ifstream inputfilelist( filelist );
                int line = 1;
                std::string file;
                do {
                  std::getline(inputfilelist, file);
                  inputfilenames_500->push_back(file);
                  ++line;
                } while (line <= NUMBER_OF_FILES_500);
                inputfiles_500_set = true;
              }

            } else {
              std::cerr << "You didn't give an argument for the pair inputfile(s)!" << std::endl;
            }
          }
        }
				if (!inputfiles_250_set || !NUMBER_OF_FILES_250_set || !inputfiles_500_set || !NUMBER_OF_FILES_500_set) {
								std::cerr
												<< "You didn't give the name for the inputfiles, or the amount of pair inputfiles. Please try again!"
												<< std::endl;
								exit(1);
				}

  double weight(1);
  weight = (double)(1.0/1312.0);
  //Make histogram for storing the information
	std::string const histo_title = "P_T of particles from pair background for 250 and 500 GeV;P_{T} [GeV]; Number of particles";
	TH1D* Hits_P_T_250 = new TH1D("PT_250GeV",histo_title.c_str(),50,0,2);
	TH1D* Hits_P_T_500 = new TH1D("PT_500GeV",histo_title.c_str(),50,0,2);
	
	std::string const histo_title_z = "P_z of particles from pair background for 250 and 500 GeV;P_{z} [GeV]; Number of particles";
	TH1D* Hits_P_z_250 = new TH1D("Pz_250GeV",histo_title_z.c_str(),180,0,270);
	TH1D* Hits_P_z_500 = new TH1D("Pz_500GeV",histo_title_z.c_str(),180,0,270);

  for (size_t file_iterator = 0; file_iterator < inputfilenames_250->size(); ++file_iterator) {
    TFile *file = TFile::Open(inputfilenames_250->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, "MCP");

    //Set the branches
    double mom_x = 0.0;
    double mom_y = 0.0;
    double mom_z = 0.0;
    int particle_pdg = -999;
    tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("Particle_PDG", kTRUE);
    tree->SetBranchAddress("Particle_PDG", &particle_pdg);
    tree->SetBranchStatus("Momentumx", 1);
    tree->SetBranchAddress("Momentumx", &mom_x);
    tree->SetBranchStatus("Momentumy", 1);
    tree->SetBranchAddress("Momentumy", &mom_y);
    tree->SetBranchStatus("Momentumz", 1);
    tree->SetBranchAddress("Momentumz", &mom_z);

    //Now we loop through the tree
    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      if ( particle_pdg != 11 && particle_pdg != -11) continue;
      //if ( particle_pdg != 2112) continue;
      Hits_P_T_250->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight);
      Hits_P_z_250->Fill(abs(mom_z),weight);
    }
    file->Close();
  }
  for (size_t file_iterator = 0; file_iterator < inputfilenames_500->size(); ++file_iterator) {
    TFile *file = TFile::Open(inputfilenames_500->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, "MCP");

    //Set the branches
    double mom_x = 0.0;
    double mom_y = 0.0;
    double mom_z = 0.0;
    int particle_pdg = -999;
    tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("Particle_ID", kTRUE);
    tree->SetBranchAddress("Particle_ID", &particle_pdg);
    tree->SetBranchStatus("Momentumx", 1);
    tree->SetBranchAddress("Momentumx", &mom_x);
    tree->SetBranchStatus("Momentumy", 1);
    tree->SetBranchAddress("Momentumy", &mom_y);
    tree->SetBranchStatus("Momentumz", 1);
    tree->SetBranchAddress("Momentumz", &mom_z);

    //Now we loop through the tree
    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      if ( particle_pdg != 11 && particle_pdg != -11) continue;
      //if ( particle_pdg != 2112) continue;
      Hits_P_T_500->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight);
      Hits_P_z_500->Fill(abs(mom_z),weight);
    }
    file->Close();
  }

	std::vector< TH1D* > histos_PT;
	histos_PT.push_back(Hits_P_T_250);
	histos_PT.push_back(Hits_P_T_500);
	double max=GetMinMaxForMultipleOverlappingHistograms(histos_PT,true).second;
	for(size_t iterator = 0; iterator < histos_PT.size(); ++iterator){
					histos_PT.at(iterator)->SetMinimum(0.5);
					histos_PT.at(iterator)->SetMaximum(max);
	}
	histos_PT.at(0)->SetMarkerColor(kCyan+2);
	histos_PT.at(0)->SetLineColor(kCyan+2);
	histos_PT.at(1)->SetMarkerColor(kPink-5);
	histos_PT.at(1)->SetMarkerStyle(26);
	histos_PT.at(1)->SetMarkerSize(0.8);
	histos_PT.at(1)->SetLineColor(kPink-5);

  std::vector< TH1D* > histos_Pz;
	histos_Pz.push_back(Hits_P_z_250);
	histos_Pz.push_back(Hits_P_z_500);
	max=GetMinMaxForMultipleOverlappingHistograms(histos_Pz,true).second;
	for(size_t iterator = 0; iterator < histos_Pz.size(); ++iterator){
					histos_Pz.at(iterator)->SetMinimum(0.5);
					histos_Pz.at(iterator)->SetMaximum(max);
	}
	histos_Pz.at(0)->SetMarkerColor(kCyan+2);
	histos_Pz.at(0)->SetLineColor(kCyan+2);
	histos_Pz.at(1)->SetMarkerColor(kPink-5);
	histos_Pz.at(1)->SetMarkerStyle(26);
	histos_Pz.at(1)->SetMarkerSize(0.8);
	histos_Pz.at(1)->SetLineColor(kPink-5);

	
	gStyle->SetOptStat(0);
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);
	canvas1->SetLogy(1);

	//Plot the histogram and save it
	histos_PT.at(0)->Draw();
	canvas1->Update();
	//std::vector<TPaveStats*> st_vec;
	//st_vec.push_back(new TPaveStats());
	//st_vec.at(0) = (TPaveStats*)histos_PT.at(0)->GetListOfFunctions()->FindObject("stats");
	//st_vec.at(0)->SetLineColor(0+1);
	//st_vec.at(0)->SetX1NDC(0.5); //new x start position
	//st_vec.at(0)->SetX2NDC(0.8); //new x end position
	//st_vec.at(0)->SetY1NDC(0.8); //new x start position
	//st_vec.at(0)->SetY2NDC(0.88); //new x end position
	//float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

	histos_PT.at(1)->Draw("SAMES");
  canvas1->Update();
  //st_vec.push_back(new TPaveStats());
	//st_vec.at(1)= (TPaveStats*)histos_PT.at(1)->GetListOfFunctions()->FindObject("stats");
	//st_vec.at(1)->SetLineColor(kPink-5);
	//st_vec.at(1)->SetX1NDC(0.5); //new x start position
	//st_vec.at(1)->SetX2NDC(0.8); //new x end position
	//st_vec.at(1)->SetY2NDC(st_vec.at(0)->GetY1NDC()); //new x end position
	//st_vec.at(1)->SetY1NDC(st_vec.at(1)->GetY2NDC()-boxsize); //new x start position


	TPaveText *text1 = new TPaveText(0.55,0.80,0.8,0.88,"brNDC");
	text1->SetTextFont(62);
	text1->SetTextColor(kCyan+2);
	text1->SetFillColor(0);
	text1->AddText("ILC250");
	text1->AddLine(0,0.5,1,0.5);
	std::stringstream entries_PT_250;
	entries_PT_250 << "Entries = " << (int)Hits_P_T_250->Integral();
	text1->AddText(entries_PT_250.str().c_str());
	TPaveText *text2 = new TPaveText(0.55,0.72,0.8,0.80,"brNDC");
	text2->SetTextFont(62);
	text2->SetTextColor(kPink-5);
	text2->SetFillColor(0);
	//text->AddLine(0,0.5,1,0.5);
	text2->AddText("ILC500");
	text2->AddLine(0,0.5,1,0.5);
	std::stringstream entries_PT_500;
	entries_PT_500 << "Entries = " << (int)Hits_P_T_500->Integral();
	text2->AddText(entries_PT_500.str().c_str());
	text1->Draw();
	text2->Draw();

	canvas1->Print("output/250_500_pairs_comparison_PT.pdf");
	canvas1->Print("output/250_500_pairs_comparison_PT.cxx");

  //Plot the histogram and save it
	histos_Pz.at(0)->Draw();
	canvas1->Update();
	//std::vector<TPaveStats*> st_vec2;
	//st_vec2.push_back(new TPaveStats());
	//st_vec2.at(0) = (TPaveStats*)histos_Pz.at(0)->GetListOfFunctions()->FindObject("stats");
	//st_vec2.at(0)->SetLineColor(0+1);
	//st_vec2.at(0)->SetX1NDC(0.5); //new x start position
	//st_vec2.at(0)->SetX2NDC(0.8); //new x end position
	//st_vec2.at(0)->SetY1NDC(0.8); //new x start position
	//st_vec2.at(0)->SetY2NDC(0.88); //new x end position
	//float boxsize2 = st_vec2.at(0)->GetY2NDC()-st_vec2.at(0)->GetY1NDC();

	histos_Pz.at(1)->Draw("SAMES");
  canvas1->Update();
  //st_vec2.push_back(new TPaveStats());
	//st_vec2.at(1)= (TPaveStats*)histos_Pz.at(1)->GetListOfFunctions()->FindObject("stats");
	//st_vec2.at(1)->SetLineColor(kPink-5);
	//st_vec2.at(1)->SetX1NDC(0.5); //new x start position
	//st_vec2.at(1)->SetX2NDC(0.8); //new x end position
	//st_vec2.at(1)->SetY2NDC(st_vec2.at(0)->GetY1NDC()); //new x end position
	//st_vec2.at(1)->SetY1NDC(st_vec2.at(1)->GetY2NDC()-boxsize2); //new x start position

	TPaveText *text3 = new TPaveText(0.55,0.80,0.8,0.88,"brNDC");
	text3->SetTextFont(62);
	text3->SetTextColor(kCyan+2);
	text3->SetFillColor(0);
	text3->AddText("ILC250");
	text3->AddLine(0,0.5,1,0.5);
	std::stringstream entries_Pz_250;
	entries_Pz_250 << "Entries = " << (int)Hits_P_z_250->Integral();
	text3->AddText(entries_Pz_250.str().c_str());
	TPaveText *text4 = new TPaveText(0.55,0.72,0.8,0.80,"brNDC");
	text4->SetTextFont(62);
	text4->SetTextColor(kPink-5);
	text4->SetFillColor(0);
	//text->AddLine(0,0.5,1,0.5);
	text4->AddText("ILC500");
	text4->AddLine(0,0.5,1,0.5);
	std::stringstream entries_Pz_500;
	entries_Pz_500 << "Entries = " << (int)Hits_P_z_500->Integral();
	text4->AddText(entries_Pz_500.str().c_str());
	text3->Draw();
	text4->Draw();

	canvas1->Print("output/250_500_pairs_comparison_Pz.pdf");
	canvas1->Print("output/250_500_pairs_comparison_Pz.cxx");

	return 0;
}

