#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1.h"

#include <iomanip>
#include <iostream>
#include <vector>
#include <map>

#include "UsefulFunctions.h"
#include "Style.h"

//enum EEntryStatus {
//  kEntryValid = 0, // data read okay
//  kEntryNotLoaded, // no entry has been loaded yet
//  kEntryNoTree, // the tree does not exist
//  kEntryNotFound, // the tree entry number does not exist
//  kEntryChainSetupError, // problem in accessing a chain element, e.g. file without the tree
//  kEntryChainFileError, // problem in opening a chain's file
//  kEntryDictionaryError, // problem reading dictionary info from tree
//};

bool Check_Set_InputArguments(int const argc, char const * const * const argv, std::string &in, std::string &out, float &ap, float &up, float &lo);
bool Check_TTreeReader_EntryStatus(TTreeReader const & redr);
bool CheckValue(ROOT::Internal::TTreeReaderValueBase* value); 

int main(int const argc, char const * const * const argv){
  UsePhDStyle();

  std::string inputfilename;
  std::string outputfilename;
  float aperture = 0.0;
  float upperjaw = 0.0;
  float lowerjaw = 0.0;
  
  if (Check_Set_InputArguments(argc, argv, inputfilename, outputfilename, aperture, upperjaw, lowerjaw) == false) return -1;

  TFile* inputfile = new TFile(inputfilename.c_str());

  // Create a TTreeReader named "MyTree" from the given TDirectory.
  // The TTreeReader gives access to the TTree to the TTreeReaderValue and
  // TTreeReaderArray objects. It knows the current entry number and knows
  // how to iterate through the TTree.
  TTreeReader reader("Event", inputfile);

  if(reader.IsZombie()){
    std::cerr << "Tree or filename is invalid" << std::endl;
    exit(-3);
  }

  // Read a single float value in each tree entries:
  TTreeReaderValue< std::vector< float > > x (reader, "RHUL_detector2.x");
  TTreeReaderValue< std::vector< float > > y (reader, "RHUL_detector2.y");

  // Variables for ROOT output file:
  int BkgLevel = 0;
  float Coll_Aperture = aperture;
  float Coll_UpperJaw_Position = upperjaw;
  float Coll_LowerJaw_Position = lowerjaw;


  // Now iterate through the TTree entries and fill a histogram.
  while (reader.Next()) {
    if (Check_TTreeReader_EntryStatus(reader) == false) return -1; 

    for(auto i = x->begin(); i != x->end(); ++i){
      //std::cout << *i << std::endl;
      //std::cout << std::fixed << std::setprecision(3) << (*x)[*i] << std::endl;
      for(auto j = y->begin(); j != y->end(); ++j){
        if(*i >= -0.03  && *i <= 0.0 
            && *j >= -0.015 && *j <= 0.015 ){
          BkgLevel++;
        }
      }
    }
  } // TTree entry / event loop
  inputfile->Close();
  
  TFile* OutputROOTFile = new TFile(outputfilename.c_str(),"CREATE","RHUL_Cherenkov_detector_signal_simulation");
  TTree* Detector1 = new TTree("Tree_Detector1","TTree for detector 1");
  
  Detector1->Branch("CollAperture",&Coll_Aperture,"CollAperture/F");
  Detector1->Branch("CollUpperJawPosition",&Coll_UpperJaw_Position,"CollUpperJawPosition/F");
  Detector1->Branch("CollLowerJawPosition",&Coll_LowerJaw_Position,"CollLowerJawPosition/F");
  Detector1->Branch("Signal",&BkgLevel,"Signal/I");
  
  Detector1->Fill();
  OutputROOTFile->Write();
  OutputROOTFile->Close();

  return 0;
}

bool Check_TTreeReader_EntryStatus(TTreeReader const & redr){
  if (redr.GetEntryStatus() == TTreeReader::kEntryValid) {
    //std::cout << "Loaded entry " << redr.GetCurrentEntry() << '\n';
    return true;
  } else { 
    if(redr.GetEntryStatus() == TTreeReader::kEntryNotLoaded){
      std::cerr << "Error: TTreeReader has not loaded any data yet!\n";
    } else if(redr.GetEntryStatus() == TTreeReader::kEntryNoTree){
      std::cerr << "Error: TTreeReader cannot find a tree names \"MyTree\"!\n";
    } else if(redr.GetEntryStatus() == TTreeReader::kEntryNotFound){
      // Can't really happen as Next() TTreeReader::knows when to stop.
      std::cerr << "Error: The entry number doe not exist\n";
    } else if(redr.GetEntryStatus() == TTreeReader::kEntryChainSetupError){
      std::cerr << "Error: TTreeReader cannot access a chain element, e.g. file without the tree\n";
    } else if(redr.GetEntryStatus() == TTreeReader::kEntryChainFileError){
      std::cerr << "Error: TTreeReader cannot open a chain element, e.g. missing file\n";
    } else if(redr.GetEntryStatus() == TTreeReader::kEntryDictionaryError){
      std::cerr << "Error: TTreeReader cannot find the dictionary for some data\n";
    } else{
      std::cerr << "Unknown Error in switch case" << std::endl;
    }
    return false;
  }
}
bool CheckValue(ROOT::Internal::TTreeReaderValueBase* value) {
  if (value->GetSetupStatus() < 0) {
    std::cerr << "Error " << value->GetSetupStatus()
      << "setting up reader for " << value->GetBranchName() << '\n';
    return false;
  }
  return true;
}
bool Check_Set_InputArguments(int const argc, char const * const * const argv, std::string &in, std::string &out, float &ap, float &up, float &lo)
{
  bool inputfilename_set = false;
  bool outputfilename_set = false;
  bool aperture_set = false;
  bool upperjaw_set = false;
  bool lowerjaw_set = false;
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-o")
					&& argv[i + 1] != std::string("-u")
					&& argv[i + 1] != std::string("-l")
					&& argv[i + 1] != std::string("-a")){
				in = argv[i + 1];
				inputfilename_set = true;
			} else {
				std::cerr << "You didn't give an argument for the inputfilename!"
					<< std::endl;
			}
		}
		if (argv[i] == std::string("-o")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-a")
					&& argv[i + 1] != std::string("-u")
					&& argv[i + 1] != std::string("-l")
					&& argv[i + 1] != std::string("-i")) {
				out = argv[i + 1];
				outputfilename_set = true;
			} else {
				std::cerr << "You didn't give an argument for the outputfilename!"
					<< std::endl;
			}
		}
		if (argv[i] == std::string("-a")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-o")
					&& argv[i + 1] != std::string("-u")
					&& argv[i + 1] != std::string("-l")
					&& argv[i + 1] != std::string("-i")) {
				ap = std::stof(argv[i + 1]);
				aperture_set = true;
			} else {
				std::cerr << "You didn't give an argument for the collimator aperture!"
					<< std::endl;
			}
		}
		if (argv[i] == std::string("-u")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-o")
					&& argv[i + 1] != std::string("-a")
					&& argv[i + 1] != std::string("-l")
					&& argv[i + 1] != std::string("-i")) {
				up = std::stof(argv[i + 1]);
				upperjaw_set = true;
			} else {
				std::cerr << "You didn't give an argument for the upper jaw position!"
					<< std::endl;
			}
		}
		if (argv[i] == std::string("-l")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-o")
					&& argv[i + 1] != std::string("-a")
					&& argv[i + 1] != std::string("-u")
					&& argv[i + 1] != std::string("-i")) {
				lo = std::stof(argv[i + 1]);
				lowerjaw_set = true;
			} else {
				std::cerr << "You didn't give an argument for the lower jaw position!"
					<< std::endl;
			}
		}

	}
	if(!inputfilename_set || !outputfilename_set || !aperture_set || !lowerjaw_set || !upperjaw_set){
		std::cerr << "You didn't set the input arguments correctly!" << std::endl;
		std::cerr << "Try:\n";
		std::cerr << "./CollimatorBkgLevelAna -i inputfilename -o outputfilename -a 24 -u 12 -l 12\n";
		std::cerr << "-a: Collimator Apertur\n";
		std::cerr << "-u: Upper jaw position\n";
		std::cerr << "-l: Lower jaw position" << std::endl;
		return false;
	}
  else return true;
}
