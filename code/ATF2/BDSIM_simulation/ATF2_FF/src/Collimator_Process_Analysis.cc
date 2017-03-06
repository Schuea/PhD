#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1.h"

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

bool Check_TTreeReader_EntryStatus(TTreeReader const & redr);
bool CheckValue(ROOT::Internal::TTreeReaderValueBase* value); 

int main(){
  UsePhDStyle();
  std::string inputfilename = "test1_02032017_event.root";
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
  TTreeReaderValue< std::map< int, int > > trackIndex_modelIndex(reader, "Trajectory.trackIndex_modelIndex");

  // Make histogramms:

  // Now iterate through the TTree entries and fill a histogram.
  while (reader.Next()) {
   if (Check_TTreeReader_EntryStatus(reader) == false) return -1; 

    for(auto i = trackIndex_modelIndex->begin(); i != trackIndex_modelIndex->end(); ++i){
      std::cout << i->second << std::endl;
    }
  } // TTree entry / event loop
  return 0;
}

bool Check_TTreeReader_EntryStatus(TTreeReader const & redr){
  if (redr.GetEntryStatus() == TTreeReader::kEntryValid) {
    std::cout << "Loaded entry " << redr.GetCurrentEntry() << '\n';
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
