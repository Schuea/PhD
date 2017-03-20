#include "TFile.h"
#include "TChain.h"
#include "TObject.h"
#include "TH1.h"


#include <iostream>
#include <fstream>
#include <sstream>

int main(int const argc, char const * const * const argv){

  std::string inputfilename;
  bool file_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")) {
      if (argv[i + 1] != NULL) 
      {
        inputfilename = argv[i + 1];
        file_set = true;
      } else {
        std::cerr << "You didn't give an argument for the inputfilename!"
          << std::endl;
        std::exit(1);
      }
    }
      }
  if (!file_set){
    std::cerr << "You didn't give the required arguments: -i inputfile.root!"
      << std::endl;
    std::exit(1);
  }


  TFile* inputfile = TFile::Open(inputfilename.c_str(),"UPDATE");
  TTree* Detector = nullptr;
  inputfile->GetObject("Tree_Detector1",Detector);

  int ismc = 0;
  TBranch* IsMC = Detector->Branch("IsMC",&ismc, "IsMC/I");
  
  long long int entries = Detector->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    ismc = 0;
    IsMC->Fill();
  }
  Detector->Write("",TObject::kOverwrite);
  inputfile->Write();
  inputfile->Close();

  return 0;
}
