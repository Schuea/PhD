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
#include "Style.h"

using namespace std;

void BinLogX(TH1D*h);

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();
				TH1::SetDefaultSumw2();

				std::vector<std::string> *inputfilenames_neutrons = new std::vector<std::string>();

				int NUMBER_OF_FILES_neutrons = 0;
				bool NUMBER_OF_FILES_neutrons_set = false;

				bool inputfiles_neutrons_set = false;

        for (int i = 1; i < argc; i++) {
          if (argv[i] == std::string("-n")) {
            if (argv[i + 1] != NULL && argv[i + 1] != std::string("-i") 
                ) {
              NUMBER_OF_FILES_neutrons = std::stoi(argv[i + 1]);
              std::cout << "Number of pair input files = " << NUMBER_OF_FILES_neutrons << std::endl;
              NUMBER_OF_FILES_neutrons_set = true;
            } else {
              std::cerr << "You didn't give an argument for the number of pair files!" << std::endl;
            }
          }
        }  
        for (int i = 1; i < argc; i++) {
          if (argv[i] == std::string("-i") 
              && argv[i + 1] != std::string("-n")) {
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
                  inputfilenames_neutrons->push_back(file);
                  ++line;
                } while (line <= NUMBER_OF_FILES_neutrons);
                inputfiles_neutrons_set = true;
              }

            } else {
              std::cerr << "You didn't give an argument for the pair inputfile(s)!" << std::endl;
            }
          }
        }
				if (!inputfiles_neutrons_set || !NUMBER_OF_FILES_neutrons_set) {
								std::cerr
												<< "You didn't give the name for the inputfiles, or the amount of pair inputfiles. Please try again!"
												<< std::endl;
								exit(1);
				}

  double weight(1);
  //Make histogram for storing the information
	std::string const histo_title_x = "P_x of EXT line neutrons;P_{x} [MeV]; Number of particles";
	TH1D* Hits_P_x_neutrons = new TH1D("Px_neutrons",histo_title_x.c_str(),75,-50,50);

	std::string const histo_title_y = "P_y of EXT line neutrons;P_{y} [MeV]; Number of particles";
	TH1D* Hits_P_y_neutrons = new TH1D("Py_neutrons",histo_title_y.c_str(),75,-50,50);
	
	std::string const histo_title_z = "P_z of EXT line neutrons;P_{z} [MeV]; Number of particles";
	TH1D* Hits_P_z_neutrons = new TH1D("Pz_neutrons",histo_title_z.c_str(),30,0,70);

  for (size_t file_iterator = 0; file_iterator < inputfilenames_neutrons->size(); ++file_iterator) {
    TFile *file = TFile::Open(inputfilenames_neutrons->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, "MCP");
    //TTree *tree1 = Get_TTree(file, "MuonEndcap");
    //TTree *tree2 = Get_TTree(file, "MuonBarrel");
    //TTree *tree3 = Get_TTree(file, "BeamCal");

    //Set the branches
    double mom_x = 0.0;
    double mom_y = 0.0;
    double mom_z = 0.0;
    int particle_pdg = -999;
    bool InCalo = 0;
    tree->SetBranchStatus("*", 0);
    //tree1->SetBranchStatus("*", 0);
    //tree2->SetBranchStatus("*", 0);
    //tree3->SetBranchStatus("*", 0);

    tree->SetBranchStatus( "DecayedInCalo_Status", kTRUE);
    tree->SetBranchAddress("DecayedInCalo_Status", &InCalo);
    tree->SetBranchStatus("Particle_PDG", kTRUE);
    tree->SetBranchAddress("Particle_PDG", &particle_pdg);
    tree->SetBranchStatus("Momentumx", 1);
    tree->SetBranchAddress("Momentumx", &mom_x);
    tree->SetBranchStatus("Momentumy", 1);
    tree->SetBranchAddress("Momentumy", &mom_y);
    tree->SetBranchStatus("Momentumz", 1);
    tree->SetBranchAddress("Momentumz", &mom_z);
    
    //tree1->SetBranchStatus("HitMotherParticle_PDG", kTRUE);
    //tree1->SetBranchAddress("HitMotherParticle_PDG", &particle_pdg);
    //tree1->SetBranchStatus("HitMotherMomentum_x", 1);
    //tree1->SetBranchAddress("HitMotherMomentum_x", &mom_x);
    //tree1->SetBranchStatus("HitMotherMomentum_y", 1);
    //tree1->SetBranchAddress("HitMotherMomentum_y", &mom_y);
    //tree1->SetBranchStatus("HitMotherMomentum_z", 1);
    //tree1->SetBranchAddress("HitMotherMomentum_z", &mom_z);
    //tree2->SetBranchStatus("HitMotherParticle_PDG", kTRUE);
    //tree2->SetBranchAddress("HitMotherParticle_PDG", &particle_pdg);
    //tree2->SetBranchStatus("HitMotherMomentum_x", 1);
    //tree2->SetBranchAddress("HitMotherMomentum_x", &mom_x);
    //tree2->SetBranchStatus("HitMotherMomentum_y", 1);
    //tree2->SetBranchAddress("HitMotherMomentum_y", &mom_y);
    //tree2->SetBranchStatus("HitMotherMomentum_z", 1);
    //tree2->SetBranchAddress("HitMotherMomentum_z", &mom_z);
    //tree3->SetBranchStatus("HitMotherParticle_PDG", kTRUE);
    //tree3->SetBranchAddress("HitMotherParticle_PDG", &particle_pdg);
    //tree3->SetBranchStatus("HitMotherMomentum_x", 1);
    //tree3->SetBranchAddress("HitMotherMomentum_x", &mom_x);
    //tree3->SetBranchStatus("HitMotherMomentum_y", 1);
    //tree3->SetBranchAddress("HitMotherMomentum_y", &mom_y);
    //tree3->SetBranchStatus("HitMotherMomentum_z", 1);
    //tree3->SetBranchAddress("HitMotherMomentum_z", &mom_z);

    //Now we loop through the tree
    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      if ( particle_pdg != 2112 || !InCalo) continue;
      Hits_P_x_neutrons->Fill(mom_x*1000.0,weight);//Plot in MeV
      Hits_P_y_neutrons->Fill(mom_y*1000.0,weight);//Plot in MeV
      Hits_P_z_neutrons->Fill(abs(mom_z)*1000.0,weight);//Plot in MeV
    }
    ////Now we loop through the tree
    //long long int const entries1 = tree1->GetEntries();
    //for (long long int i = 0; i < entries1; ++i) {
    //  tree1->GetEntry(i);
    //  if ( particle_pdg != 2112) continue;
    //  Hits_P_x_neutrons->Fill(mom_x*1000.0,weight);//Plot in MeV
    //  Hits_P_y_neutrons->Fill(mom_y*1000.0,weight);//Plot in MeV
    //  Hits_P_z_neutrons->Fill(abs(mom_z)*1000.0,weight);//Plot in MeV
    //}
    ////Now we loop through the tree
    //long long int const entries2 = tree2->GetEntries();
    //for (long long int i = 0; i < entries2; ++i) {
    //  tree2->GetEntry(i);
    //  if ( particle_pdg != 2112) continue;
    //  Hits_P_x_neutrons->Fill(mom_x*1000.0,weight);//Plot in MeV
    //  Hits_P_y_neutrons->Fill(mom_y*1000.0,weight);//Plot in MeV
    //  Hits_P_z_neutrons->Fill(abs(mom_z)*1000.0,weight);//Plot in MeV
    //}
    ////Now we loop through the tree
    //long long int const entries3 = tree3->GetEntries();
    //for (long long int i = 0; i < entries3; ++i) {
    //  tree3->GetEntry(i);
    //  if ( particle_pdg != 2112) continue;
    //  Hits_P_x_neutrons->Fill(mom_x*1000.0,weight);//Plot in MeV
    //  Hits_P_y_neutrons->Fill(mom_y*1000.0,weight);//Plot in MeV
    //  Hits_P_z_neutrons->Fill(abs(mom_z)*1000.0,weight);//Plot in MeV
    //}


    file->Close();
  }
	std::vector< TH1D* > histos_px;
	histos_px.push_back(Hits_P_x_neutrons);
	double max=GetMinMaxForMultipleOverlappingHistograms(histos_px,true).second;
	for(size_t iterator = 0; iterator < histos_px.size(); ++iterator){
					//histos_px.at(iterator)->SetMinimum(0.5);
					histos_px.at(iterator)->SetMaximum(max);
	}
	histos_px.at(0)->SetMarkerColor(kCyan+2);
	histos_px.at(0)->SetLineColor(kCyan+2);

	std::vector< TH1D* > histos_py;
	histos_py.push_back(Hits_P_y_neutrons);
	max=GetMinMaxForMultipleOverlappingHistograms(histos_py,true).second;
	for(size_t iterator = 0; iterator < histos_py.size(); ++iterator){
					//histos_py.at(iterator)->SetMinimum(0.5);
					histos_py.at(iterator)->SetMaximum(max);
	}
	histos_py.at(0)->SetMarkerColor(kCyan+2);
	histos_py.at(0)->SetLineColor(kCyan+2);

  std::vector< TH1D* > histos_Pz;
	histos_Pz.push_back(Hits_P_z_neutrons);
	max=GetMinMaxForMultipleOverlappingHistograms(histos_Pz,true).second;
	for(size_t iterator = 0; iterator < histos_Pz.size(); ++iterator){
					histos_Pz.at(iterator)->SetMinimum(0.5);
					histos_Pz.at(iterator)->SetMaximum(max);
	}
	histos_Pz.at(0)->SetMarkerColor(kCyan+2);
	histos_Pz.at(0)->SetLineColor(kCyan+2);
	
	gStyle->SetOptStat(111);
	//gStyle->SetOptStat(0);
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);
	canvas1->SetLogy(1);
	canvas1->SetLogx(0);
	canvas1->SetGrid(1);

	//Plot the histogram and save it
	histos_px.at(0)->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec;
	st_vec.push_back(new TPaveStats());
	st_vec.at(0) = (TPaveStats*)histos_px.at(0)->GetListOfFunctions()->FindObject("stats");
	st_vec.at(0)->SetLineColor(0+1);
	st_vec.at(0)->SetX1NDC(0.6); //new x start position
	st_vec.at(0)->SetX2NDC(0.8); //new x end position
	st_vec.at(0)->SetY1NDC(0.8); //new x start position
	st_vec.at(0)->SetY2NDC(0.88); //new x end position
	//float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();


	//TPaveText *text1 = new TPaveText(0.55,0.80,0.8,0.88,"brNDC");
	//text1->SetTextFont(62);
	//text1->SetTextColor(kCyan+2);
	//text1->SetFillColor(0);
	//text1->AddText("ILCneutrons");
	//text1->AddLine(0,0.5,1,0.5);
	//std::stringstream entries_px_neutrons;
	//entries_px_neutrons << "Entries = " << (int)histos_px.at(0)->Integral();
	//text1->AddText(entries_px_neutrons.str().c_str());
	//text1->Draw();

	canvas1->Print("output/neutrons_Px.pdf");
	canvas1->Print("output/neutrons_Px.cxx");

  //Plot the histogram and save it
	histos_py.at(0)->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec2;
	st_vec2.push_back(new TPaveStats());
	st_vec2.at(0) = (TPaveStats*)histos_py.at(0)->GetListOfFunctions()->FindObject("stats");
	st_vec2.at(0)->SetLineColor(0+1);
	st_vec2.at(0)->SetX1NDC(0.6); //new x start position
	st_vec2.at(0)->SetX2NDC(0.8); //new x end position
	st_vec2.at(0)->SetY1NDC(0.8); //new x start position
	st_vec2.at(0)->SetY2NDC(0.88); //new x end position

	//TPaveText *text2 = new TPaveText(0.55,0.80,0.8,0.88,"brNDC");
	//text2->SetTextFont(62);
	//text2->SetTextColor(kCyan+2);
	//text2->SetFillColor(0);
	//text2->AddText("ILCneutrons");
	//text2->AddLine(0,0.5,1,0.5);
	//std::stringstream entries_py_neutrons;
	//entries_py_neutrons << "Entries = " << (int)histos_py.at(0)->Integral();
	//text2->AddText(entries_py_neutrons.str().c_str());
	//text2->Draw();

	canvas1->Print("output/neutrons_Py.pdf");
	canvas1->Print("output/neutrons_Py.cxx");

	//canvas1->SetLogx(1);
  //Plot the histogram and save it
	//BinLogX( histos_Pz.at(0) );
	histos_Pz.at(0)->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec3;
	st_vec3.push_back(new TPaveStats());
	st_vec3.at(0) = (TPaveStats*)histos_Pz.at(0)->GetListOfFunctions()->FindObject("stats");
	st_vec3.at(0)->SetLineColor(0+1);
	st_vec3.at(0)->SetX1NDC(0.6); //new x start position
	st_vec3.at(0)->SetX2NDC(0.8); //new x end position
	st_vec3.at(0)->SetY1NDC(0.8); //new x start position
	st_vec3.at(0)->SetY2NDC(0.88); //new x end position

	//TPaveText *text3 = new TPaveText(0.55,0.80,0.8,0.88,"brNDC");
	//text3->SetTextFont(62);
	//text3->SetTextColor(kCyan+2);
	//text3->SetFillColor(0);
	//text3->AddText("ILCneutrons");
	//text3->AddLine(0,0.5,1,0.5);
	//std::stringstream entries_Pz_neutrons;
	//entries_Pz_neutrons << "Entries = " << (int)Hits_P_z_neutrons->Integral();
	//text3->AddText(entries_Pz_neutrons.str().c_str());

	canvas1->Print("output/neutrons_Pz.pdf");
	canvas1->Print("output/neutrons_Pz.cxx");

	return 0;
}

void BinLogX(TH1D*h)
{

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = pow(10, from + i * width);

  }
  axis->Set(bins, new_bins);
  delete new_bins;
} 
