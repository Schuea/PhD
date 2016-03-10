#include "TFile.h"
#include "TChain.h"
#include "TObject.h"

#include <iostream>
#include <fstream>
#include <sstream>

int main(int const argc, char const * const * const argv){

  std::vector<std::string> inputfilenames;
  std::string outputfilename;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-o")){
        int j = 1;
        while(argv[i + j] != NULL
            && argv[i + j] != std::string("-o")){
          inputfilenames.push_back(argv[i + j]);
          ++j;
        }
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

  float Coll_Aperture = 0.0;
  int voltage1 = 0;
  int signal1 = 0;

  TChain* ROOTChain = new TChain("Tree_Detector1");

  ROOTChain->SetBranchAddress("CollAperture",&Coll_Aperture);
  ROOTChain->SetBranchAddress("Voltage",&voltage1);
  ROOTChain->SetBranchAddress("Signal",&signal1);

  for(size_t inputfile_it = 0; inputfile_it < inputfilenames.size(); ++inputfile_it){
    ROOTChain->Add(inputfilenames.at(inputfile_it).c_str());
  }  
  ROOTChain->Merge(outputfilename.c_str());

  return 0;
}
