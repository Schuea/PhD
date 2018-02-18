#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TLegend.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include "Subdetector_class_new.h"
#include "UsefulFunctions.h"
#include "GeneralFunctions_SiDBkgSim_new.h"
#include "Style.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
	UsePhDStyle();
  
  std::vector< std::string > *inputfilenames = new std::vector< std::string >();
  double weight = 0.0;
  std::vector< std::string > *argument_subdetectors = new std::vector< std::string >();

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
  std::string outputfile_name = oss.str();

  bool weights_set = false;
  bool inputfile_set = false;
  bool subdetector_set = false;

  for (int i = 1; i < argc; i++) {
	  if (argv[i] == std::string("-i")) {
		  if (argv[i + 1] != NULL && 
				  argv[i + 1] != std::string("-s") && 
				  argv[i + 1] != std::string("-w") && 
				  argv[i + 1] != std::string("-o")) {
			  int j = 1;
			  while (argv[i + j] != NULL && 
					  argv[i + j] != std::string("-w") && 
					  argv[i + j] != std::string("-i") && 
					  argv[i + j] != std::string("-o") && 
					  argv[i + j] != std::string("-s")) {
				  if( access( argv[i + j], F_OK ) != -1 ){
					  inputfilenames->push_back( argv[i + j] );
					  j++;
				  }
				  else{
					  std::cerr
						  << "The inputfiles " << argv[i + j] << " does not exist!"
						  << std::endl;
					  exit(1);
				  }
			  }
			  inputfile_set = true;
		  } else {
			  std::cerr << "You didn't give arguments for the inputfile(s)!" << std::endl;
		  }
	  }
    else if (argv[i] == std::string("-w")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-s") && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-o")) {
          weight = atof(argv[i + 1]);
          weights_set = true;
      }
      else{
        std::cerr << "You didn't give an argument for the weight of the inputfile(s)!" << std::endl;
      }
    }  
    else if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-w") && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-o")) {
        int j = 1;
        while (argv[i + j] != NULL && 
            argv[i + j] != std::string("-w") && 
            argv[i + j] != std::string("-i") && 
            argv[i + j] != std::string("-o") && 
            argv[i + j] != std::string("-s")) {
          argument_subdetectors->push_back( argv[i + j] );
          j++;
        }
        subdetector_set = true;
      } 
      else {
        std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
      }
    }
	  else if (argv[i] == std::string("-o")) {
		  if (argv[i + 1] != NULL && 
				  argv[i + 1] != std::string("-w") && 
				  argv[i + 1] != std::string("-i") && 
				  argv[i + 1] != std::string("-s")) {
			  outputfile_name = argv[i + 1];
		  } else {
			  std::cerr << "You didn't give an argument for the outputfile name!" << std::endl;
		  }
	  }
  }
  if (!inputfile_set || !subdetector_set || !weights_set) {
	  std::cerr
		  << "You didn't give the name for the subdector, the inputfiles or the weights. Please try again!"
		  << std::endl;
	  exit(1);
  }

	//Make histogram for storing the information
	//std::vector< TH2D* > Hits_Time_rtime_2D_;
	float timemax = 100.0; //ns (one bunch spacing is 554 ns)
	float timemin = 0.;
	int timerange = (int)(timemax - timemin);

	std::vector<TH1D*> Hits_Time_;
	std::string const histo_name_All = "HitTime_All";
	std::string const histo_title_All = "Hit time for all subdetectors;Hit time [ns];Number of hits";
  TH1D *Hits_Time_All = new TH1D(histo_name_All.c_str(), histo_title_All.c_str(), timerange, timemin, timemax);

  for (size_t subdetector_iterator = 0; subdetector_iterator < argument_subdetectors->size(); ++subdetector_iterator) {
    Subdetector det( argument_subdetectors->at(subdetector_iterator) );

    bool endcap = false;
    bool barrel = false;
    if (det.getShape().find("endcap") != std::string::npos){//If Endcap, calculate CellID with x and y
      std::cout <<"Found 'endcap'!" << std::endl;
      endcap = true;
    }
    else if (det.getShape().find("barrel") != std::string::npos){
      std::cout <<"Found 'barrel'!" << std::endl;
      barrel = true;
    }
    else{
      std::cerr << "The given subdetector shape was not recognized!" << std::endl; 
      exit(-1);
    }
    bool Silicon = false;
    bool Calo = false;
    if (det.getName().find("Si") != std::string::npos){
      std::cout << "Silicon detector found!" << std::endl;
      Silicon = true;
    }
    else if (det.getName().find("cal",0) != std::string::npos
        || det.getName().find("Cal",0) != std::string::npos
        || det.getName().find("Muon",0) != std::string::npos){
      std::cout << "Calorimeter found!" << std::endl;
      Calo = true;
    }
    else{
      std::cerr << "The given subdetector name was not recognized!" << std::endl; 
      exit(-1);
    }

	  std::string const histo_name = det.getName();
	  std::string const histo_title = "Hit time;Hit time [ns];Number of hits";
    TH1D *temp = new TH1D(histo_name.c_str(), histo_title.c_str(), timerange, timemin, timemax);
    Hits_Time_.push_back( temp );

    for (size_t file_iterator = 0; file_iterator < inputfilenames->size(); ++file_iterator) {
      TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
      TTree *tree = Get_TTree(file, det.getName());

      //Set the branches
      float actualtime = 0.0;
      tree->SetBranchStatus("*", 0);

      if (Calo) {
							tree->SetBranchStatus("HitContrTime", 1);
							tree->SetBranchAddress("HitContrTime", &actualtime);
			}
			else if (Silicon){
							tree->SetBranchStatus("HitTime", 1);
							tree->SetBranchAddress("HitTime", &actualtime);
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}


			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				Hits_Time_.at(subdetector_iterator)->Fill(actualtime-14, weight);
				Hits_Time_All->Fill(actualtime-14, weight);
			}
			file->Close();
		}
	}

	//Plot the histogram and save it
	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);
  //std::vector< TPaveStats* > *st_vec = new std::vector< TPaveStats* >();
  //float boxsize = 0.0;

	canvas1->cd();
	canvas1->SetLogy();
	gStyle->SetOptStat(0);
  TLegend* leg = new TLegend(0.6,0.5,0.9,0.9);
  leg->SetMargin(0.1);
  for (size_t subdetector_iterator = 0; subdetector_iterator < argument_subdetectors->size(); ++subdetector_iterator) {
    Hits_Time_.at(subdetector_iterator)->SetLineColor(subdetector_iterator+1);
    Hits_Time_.at(subdetector_iterator)->Sumw2(1);
    if(subdetector_iterator == 0){
      Hits_Time_.at(0)->Draw("e,hist");
      //canvas1->Update();
      //st_vec->push_back(new TPaveStats());
      //st_vec->at(0) = (TPaveStats*)Hits_Time_.at(0)->GetListOfFunctions()->FindObject("stats");
      //st_vec->at(0)->SetLineColor(0+1);
      //st_vec->at(0)->SetX1NDC(0.65); //new x start position
      //st_vec->at(0)->SetX2NDC(0.85); //new x end position
      //st_vec->at(0)->SetY1NDC(0.83); //new x start position
      //st_vec->at(0)->SetY2NDC(0.9); //new x end position
      //boxsize = st_vec->at(0)->GetY2NDC()-st_vec->at(0)->GetY1NDC();
    }
    else{
      Hits_Time_.at(subdetector_iterator)->Draw("e,hist,SAMES");
      //canvas1->Update();
      //st_vec->push_back(new TPaveStats());
      //st_vec->at(subdetector_iterator)= (TPaveStats*)Hits_Time_.at(subdetector_iterator)->GetListOfFunctions()->FindObject("stats");
      //st_vec->at(subdetector_iterator)->SetLineColor(subdetector_iterator+1);
      //st_vec->at(subdetector_iterator)->SetX1NDC(0.65); //new x start position
      //st_vec->at(subdetector_iterator)->SetX2NDC(0.85); //new x end position
      //st_vec->at(subdetector_iterator)->SetY2NDC(st_vec->at(subdetector_iterator-1)->GetY1NDC()); //new x end position
      //st_vec->at(subdetector_iterator)->SetY1NDC(st_vec->at(subdetector_iterator)->GetY2NDC()-boxsize); //new x start position
    }
      leg->AddEntry(Hits_Time_.at(subdetector_iterator),Hits_Time_.at(subdetector_iterator)->GetName(),"l");
  }
  leg->Draw();
	canvas1->Print(("output/hittime_"+outputfile_name+"_superimposed.pdf").c_str());
	canvas1->Print(("output/hittime_"+outputfile_name+"_superimposed.cxx").c_str());

  TLegend* leg2 = new TLegend(0.7,0.8,0.9,0.9);
  leg2->SetMargin(0.1);
	gStyle->SetOptStat(0);
  Hits_Time_All->Sumw2(1);
  Hits_Time_All->SetLineColor(kPink-1);
  Hits_Time_All->Draw("e,hist");
  leg2->SetHeader("Hit time","C"); // option "C" allows to center the header
  std::ostringstream entries;
  entries << "    Entries ";
  entries << (int)(Hits_Time_All->GetEntries()*weight);
  leg2->AddEntry(Hits_Time_All, entries.str().c_str(),"");
  leg2->Draw();
	canvas1->Print(("output/hittime_"+outputfile_name+"_All.pdf").c_str());
	canvas1->Print(("output/hittime_"+outputfile_name+"_All.cxx").c_str());
	return 0;
}

