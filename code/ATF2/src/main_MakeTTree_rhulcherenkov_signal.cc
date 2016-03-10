#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"

#include <iostream>
#include <fstream>
#include <sstream>

int main(int const argc, char const * const * const argv){
 
  std::string inputfilename;
  std::string outputfilename;
  
  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-o")){
        inputfilename = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the inputfilename!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-o")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-i")) {
        outputfilename = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the outputfilename!"
          << std::endl;
      }
    }
  }

  float Coll_aperture = 0.0;
  int voltage1 = 0;
  int signal1 = 0;

	TFile* ROOTFile = new TFile(outputfilename.c_str(),"CREATE","RHUL_Cherenkov_detector_signal");
  TTree* Detector1 = new TTree("Tree_Detector1","TTree for detector 1");
  
  Detector1->Branch("CollAperture",&Coll_aperture,"CollAperture/F");
  Detector1->Branch("Voltage",&voltage1,"Voltage/I");
  Detector1->Branch("Signal",&signal1,"Signal/I");
  
  std::ifstream inputfile(inputfilename);
  std::string line;

  while (!inputfile.eof()){
    std::getline(inputfile, line);
    if (line == 1 || line == 2){
      continue;
    }
    std::istringstream in(line);
    std::string col1, col2, col3, col4, col5;
    in >> col1 >> col2 >> col3 >> col4 >> col5;

    Coll_aperture = std::atoi(col1.c_str());
    voltage1 = std::atoi(col2.c_str());
    signal1 = std::atoi(col4.c_str());

    Detector1->Fill();
  }
  inputfile.close();
  ROOTFile->Write();
  ROOTFile->Close();

  return 0;
}
