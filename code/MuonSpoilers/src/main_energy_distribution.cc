#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"

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

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
  std::string outputfile_name = oss.str();

  bool weights_set = false;
  bool inputfile_set = false;

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
  if (!inputfile_set || !weights_set) {
    std::cerr
      << "You didn't give the name for the subdector, the inputfiles or the weights. Please try again!"
      << std::endl;
    exit(1);
  }

  //Make histogram for storing the information
  TFile* Outputfile = new TFile(("output/Energy_"+outputfile_name+".root").c_str(),"RECREATE");
  
  float energymax = 500./2. + 10.0; 
  //float energymax = atof( inputfilenames->at(0).substr(3,3).c_str() )/2.0 + 10.0; 
  float energymin = 0.;
  int energyrange = (int)((energymax - energymin)/3.0);

  std::vector<TH1D*> Hits_Time_;
  std::string const histo_name_All = "Muon_Energy";
  std::string const histo_title_All = "Muon Energy;Energy [GeV];Number of muons";
  TH1D *Energy_hist = new TH1D(histo_name_All.c_str(), histo_title_All.c_str(), energyrange, energymin, energymax);

  for (size_t file_iterator = 0; file_iterator < inputfilenames->size(); ++file_iterator) {
    TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, "MCP");

    //Set the branches
    double energy = 0.0;
    int PDG = 0;
    tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("Energy", 1);
    tree->SetBranchAddress("Energy", &energy);
    tree->SetBranchStatus("Particle_PDG", 1);
    tree->SetBranchAddress("Particle_PDG", &PDG);


    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      if (PDG == 13 || PDG == -13) Energy_hist->Fill(energy, weight);
    }
    file->Close();
  }

  //Plot the histogram and save it

  TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);

  canvas1->SetLogy();
  gStyle->SetOptStat(0);

  TLegend* leg2 = new TLegend(0.6,0.75,0.9,0.9);
  leg2->SetMargin(0.1);
  gStyle->SetOptStat(0);
  Energy_hist->Sumw2(1);
  Energy_hist->SetLineColor(kPink-1);
  Energy_hist->Draw("hist,e");
  leg2->SetHeader("Energy","C"); // option "C" allows to center the header
  std::ostringstream entries;
  entries << " Entries ";
  entries << (int)(Energy_hist->GetEntries()*weight);
  leg2->AddEntry(Energy_hist, entries.str().c_str(),"");
  leg2->SetTextSize(0.05);
  leg2->Draw();
  canvas1->Print(("output/energy_"+outputfile_name+"_All.pdf").c_str());
  canvas1->Print(("output/energy_"+outputfile_name+"_All.cxx").c_str());
  
  Outputfile->Write();
  return 0;
}

