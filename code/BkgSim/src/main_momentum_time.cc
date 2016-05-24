#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TAxis.h"
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
  TH1D* ParticleCreationTime = new TH1D("Creation time", ("Creation time of particles hitting subdetector "+subdetectornames+";Particle creation time [ns]; #hits").c_str(), 4000,0, 2000);

  std::vector< TH1D* > Momentum_histos;
	std::string const title1 = "Particle momentum in certain time interval for particles created before 1ns, hitting subdetector "
			+ subdetectornames + ";Total momentum [GeV];#hits";
  Momentum_histos.emplace_back( new TH1D("0 < hittime < 10ns", title1.c_str(), 50, 0, 1) );
  Momentum_histos.emplace_back( new TH1D("10 < hittime < 20ns", title1.c_str(), 50, 0, 1) );
  Momentum_histos.emplace_back( new TH1D("20 < hittime < 30ns", title1.c_str(), 50, 0, 1) );
  Momentum_histos.emplace_back( new TH1D("30 < hittime < 50ns", title1.c_str(), 50, 0, 1) );
	
  std::string const title2 = "Particle momentum at hittime for particles created before 1ns, hitting subdetector "
			+ subdetectornames + ";hit time [ns];Total momentum [GeV]";
	TH2D* Momentum_time = new TH2D("Momentum_time", title2.c_str(), 70,0,100, 20, 0, 0.2);

  std::vector< TH1D* > BackScatter_Momentum_histos;
	std::string const BackScatter_title1 = "Particle momentum in certain time interval for particles created after 1ns, hitting subdetector "
			+ subdetectornames + ";Total momentum [GeV];#hits";
  BackScatter_Momentum_histos.emplace_back( new TH1D("0 < hittime < 10ns", title1.c_str(), 50, 0, 1) );
  BackScatter_Momentum_histos.emplace_back( new TH1D("10 < hittime < 20ns", title1.c_str(), 50, 0, 1) );
  BackScatter_Momentum_histos.emplace_back( new TH1D("20 < hittime < 30ns", title1.c_str(), 50, 0, 1) );
  BackScatter_Momentum_histos.emplace_back( new TH1D("30 < hittime < 50ns", title1.c_str(), 50, 0, 1) );
	
  std::string const BackScatter_title2 = "Particle momentum at hittime for particles created after 1ns, hitting subdetector "
			+ subdetectornames + ";hit time [ns];Total momentum [GeV]";
	TH2D* BackScatter_Momentum_time = new TH2D("Momentum_time", BackScatter_title2.c_str(), 70,0,100, 20, 0, 0.2);


	std::stringstream subdetector_names;

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {
	  subdetector_names << SubDetectors->at(subdetector_iterator)->GetName();

		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
			TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			float actualtime = 0.0;
			float creationtime = 0.0;
			float x_cal = 0.0;
			float y_cal = 0.0;
			float z_cal = 0.0;
			double momentum_x_cal = 0.0;
			double momentum_y_cal = 0.0;
			double momentum_z_cal = 0.0;
			double x_tracker = 0.0;
			double y_tracker = 0.0;
			double z_tracker = 0.0;
			float momentum_x_tracker = 0.0;
			float momentum_y_tracker = 0.0;
			float momentum_z_tracker = 0.0;
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
							tree->SetBranchStatus("HitPosition_x", 1);
							tree->SetBranchStatus("HitPosition_y", 1);
							tree->SetBranchStatus("HitPosition_z", 1);
							tree->SetBranchStatus("HitMotherMomentum_x", 1);
							tree->SetBranchStatus("HitMotherMomentum_y", 1);
							tree->SetBranchStatus("HitMotherMomentum_z", 1);
							tree->SetBranchAddress("HitContrTime", &actualtime);
							tree->SetBranchAddress("HitPosition_x", &x_cal);
							tree->SetBranchAddress("HitPosition_y", &y_cal);
							tree->SetBranchAddress("HitPosition_z", &z_cal);
							tree->SetBranchAddress("HitMotherMomentum_x", &momentum_x_cal);
							tree->SetBranchAddress("HitMotherMomentum_y", &momentum_y_cal);
							tree->SetBranchAddress("HitMotherMomentum_z", &momentum_z_cal);
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
							tree->SetBranchStatus("HitTime", 1);
							tree->SetBranchStatus("HitParticleCreationTime", 1);
							tree->SetBranchStatus("HitPosition_x", 1);
							tree->SetBranchStatus("HitPosition_y", 1);
							tree->SetBranchStatus("HitPosition_z", 1);
							tree->SetBranchStatus("HitMomentum_x", 1);
							tree->SetBranchStatus("HitMomentum_y", 1);
							tree->SetBranchStatus("HitMomentum_z", 1);
							tree->SetBranchAddress("HitTime", &actualtime);
							tree->SetBranchAddress("HitParticleCreationTime", &creationtime);
							tree->SetBranchAddress("HitPosition_x", &x_tracker);
							tree->SetBranchAddress("HitPosition_y", &y_tracker);
							tree->SetBranchAddress("HitPosition_z", &z_tracker);
							tree->SetBranchAddress("HitMomentum_x", &momentum_x_tracker);
							tree->SetBranchAddress("HitMomentum_y", &momentum_y_tracker);
							tree->SetBranchAddress("HitMomentum_z", &momentum_z_tracker);
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}

			//Now we loop through the tree
			//Combine the two Cell ID's into a single new Cell ID
			//See how often the new Cell ID occurs in total, this is the occupancy

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
        ParticleCreationTime->Fill(creationtime);
        if( creationtime > 1){
          BackScatter_Momentum_time->Fill(actualtime, std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));

          if (actualtime < 10.0) BackScatter_Momentum_histos.at(0)->Fill(sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
          else if (actualtime >= 10.0 && actualtime < 20.0) BackScatter_Momentum_histos.at(1)->Fill(std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
          else if (actualtime >= 20.0 && actualtime < 30.0) BackScatter_Momentum_histos.at(2)->Fill(std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
          else if (actualtime >= 30.0 && actualtime < 50.0) BackScatter_Momentum_histos.at(3)->Fill(std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
        }
        //else{
          Momentum_time->Fill(actualtime, std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));

          if (actualtime < 10.0) Momentum_histos.at(0)->Fill(sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
          else if (actualtime >= 10.0 && actualtime < 20.0) Momentum_histos.at(1)->Fill(std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
          else if (actualtime >= 20.0 && actualtime < 30.0) Momentum_histos.at(2)->Fill(std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
          else if (actualtime >= 30.0 && actualtime < 50.0) Momentum_histos.at(3)->Fill(std::sqrt(std::pow(momentum_x_tracker,2)+std::pow(momentum_y_tracker,2)+std::pow(momentum_z_tracker,2)));
       // }
			}
			file->Close();
		}
	}
	gStyle->SetOptStat(1);

	//Plot the histogram and save it
	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);
	TCanvas *canvas2 = new TCanvas("canvas2", "canvas", 800, 600);
	TCanvas *canvas3 = new TCanvas("canvas3", "canvas", 800, 600);

	canvas1->cd();
  canvas1->SetLogx(0);
  canvas1->SetLogy(0);
	Momentum_time->Draw("colz");
	canvas1->Update();
	TPaveStats *st1 = (TPaveStats*)Momentum_time->GetListOfFunctions()->FindObject("stats");
	st1->SetX1NDC(0.65); //new x start position
	st1->SetX2NDC(0.85); //new x end position
	st1->SetY1NDC(0.7); //new x start position
	st1->SetY2NDC(0.9); //new x end position

	canvas1->Print(("output/momentum_time_"+subdetector_names.str()+".pdf").c_str());
	canvas1->Print(("output/momentum_time_"+subdetector_names.str()+".cxx").c_str());

  BackScatter_Momentum_time->Draw("colz");
	canvas1->Update();
	TPaveStats *BackScatter_st1 = (TPaveStats*)BackScatter_Momentum_time->GetListOfFunctions()->FindObject("stats");
	BackScatter_st1->SetX1NDC(0.65); //new x start position
	BackScatter_st1->SetX2NDC(0.85); //new x end position
	BackScatter_st1->SetY1NDC(0.7); //new x start position
	BackScatter_st1->SetY2NDC(0.9); //new x end position

	canvas1->Print(("output/backscatter_momentum_time_"+subdetector_names.str()+".pdf").c_str());
	canvas1->Print(("output/backscatter_momentum_time_"+subdetector_names.str()+".cxx").c_str());

  std::cout << canvas2 << std::endl;
	canvas2->cd();
  canvas2->SetLogx(0);
  canvas2->SetLogy(1);
  double momentummax=GetMinMaxForMultipleOverlappingHistograms(Momentum_histos,true).second;
  for(int time_iterator = 0; time_iterator < Momentum_histos.size(); ++time_iterator){
    Momentum_histos.at(time_iterator)->SetMaximum(momentummax);
  }
	Momentum_histos.at(0)->Draw();
  canvas2->Update();
	TPaveStats *st2 = (TPaveStats*)Momentum_histos.at(0)->GetListOfFunctions()->FindObject("stats");
	st2->SetX1NDC(0.65); //new x start position
	st2->SetX2NDC(0.85); //new x end position
	st2->SetY1NDC(0.73); //new x start position
	st2->SetY2NDC(0.9); //new x end position
	std::vector<TPaveStats*> st_vec;
	st_vec.push_back(st2);
	float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

  for(int time_iterator=1; time_iterator <=3; ++time_iterator){
	  Momentum_histos.at(time_iterator)->SetLineColor(time_iterator+1);
	  Momentum_histos.at(time_iterator)->Draw("SAMES");
		canvas2->Update();
		st_vec.push_back(new TPaveStats());
	  st_vec.at(time_iterator)= (TPaveStats*)Momentum_histos.at(time_iterator)->GetListOfFunctions()->FindObject("stats");
		st_vec.at(time_iterator)->SetLineColor(time_iterator+1);
		st_vec.at(time_iterator)->SetX1NDC(0.65); //new x start position
		st_vec.at(time_iterator)->SetX2NDC(0.85); //new x end position
		st_vec.at(time_iterator)->SetY2NDC(st_vec.at(time_iterator-1)->GetY1NDC()); //new x end position
		st_vec.at(time_iterator)->SetY1NDC(st_vec.at(time_iterator)->GetY2NDC()-boxsize); //new x start position
  }

	canvas2->Print(("output/momentum_histo_"+subdetector_names.str()+".pdf").c_str());
	canvas2->Print(("output/momentum_histo_"+subdetector_names.str()+".cxx").c_str());

  //****Backscatters:
  double backscattermax=GetMinMaxForMultipleOverlappingHistograms(BackScatter_Momentum_histos,true).second;
  for(int time_iterator=0; time_iterator < BackScatter_Momentum_histos.size(); ++time_iterator){
    BackScatter_Momentum_histos.at(time_iterator)->SetMaximum(backscattermax);
  }
  BackScatter_Momentum_histos.at(0)->Draw();
  canvas2->Update();
	TPaveStats *BackScatter_st2 = (TPaveStats*)BackScatter_Momentum_histos.at(0)->GetListOfFunctions()->FindObject("stats");
	BackScatter_st2->SetX1NDC(0.65); //new x start position
	BackScatter_st2->SetX2NDC(0.85); //new x end position
	BackScatter_st2->SetY1NDC(0.73); //new x start position
	BackScatter_st2->SetY2NDC(0.9); //new x end position
	std::vector<TPaveStats*> BackScatter_st_vec;
	BackScatter_st_vec.push_back(BackScatter_st2);
	float BackScatter_boxsize = BackScatter_st_vec.at(0)->GetY2NDC()-BackScatter_st_vec.at(0)->GetY1NDC();

  for(int time_iterator=1; time_iterator <=3; ++time_iterator){
	  BackScatter_Momentum_histos.at(time_iterator)->SetLineColor(time_iterator+1);
	  BackScatter_Momentum_histos.at(time_iterator)->Draw("SAMES");
		canvas2->Update();
		BackScatter_st_vec.push_back(new TPaveStats());
		BackScatter_st_vec.at(time_iterator)= (TPaveStats*)BackScatter_Momentum_histos.at(time_iterator)->GetListOfFunctions()->FindObject("stats");
		BackScatter_st_vec.at(time_iterator)->SetLineColor(time_iterator+1);
		BackScatter_st_vec.at(time_iterator)->SetX1NDC(0.65); //new x start position
		BackScatter_st_vec.at(time_iterator)->SetX2NDC(0.85); //new x end position
		BackScatter_st_vec.at(time_iterator)->SetY2NDC(BackScatter_st_vec.at(time_iterator-1)->GetY1NDC()); //new x end position
		BackScatter_st_vec.at(time_iterator)->SetY1NDC(BackScatter_st_vec.at(time_iterator)->GetY2NDC()-BackScatter_boxsize); //new x start position
  }

	canvas2->Print(("output/backscatter_momentum_histo_"+subdetector_names.str()+".pdf").c_str());
	canvas2->Print(("output/backscatter_momentum_histo_"+subdetector_names.str()+".cxx").c_str());

  canvas3->cd();
  canvas3->SetLogx(1);
  canvas3->SetLogy(1);
  ParticleCreationTime->GetXaxis()->SetRangeUser(0.1,2000);
  ParticleCreationTime->GetXaxis()->SetLimits(0.1,2000);
  ParticleCreationTime->GetYaxis()->SetRangeUser(4,5000000);//for both, SiVertexBarrel and SiVertexEndcap
  //ParticleCreationTime->GetYaxis()->SetRangeUser(4,2000000);
  ParticleCreationTime->Draw();
	canvas3->Update();
	TPaveStats *st3 = (TPaveStats*)ParticleCreationTime->GetListOfFunctions()->FindObject("stats");
	st3->SetX1NDC(0.65); //new x start position
	st3->SetX2NDC(0.85); //new x end position
	st3->SetY1NDC(0.7); //new x start position
	st3->SetY2NDC(0.9); //new x end position

	canvas3->Print(("output/creationtime_histo_"+subdetector_names.str()+".pdf").c_str());
	canvas3->Print(("output/creationtime_histo_"+subdetector_names.str()+".cxx").c_str());

	return 0;
}
