#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1.h"

#include <iostream>
#include <vector>
#include <map>

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
//bool CheckValue(ROOT::TTreeReaderValueBase* value); 
bool CheckValue(ROOT::Internal::TTreeReaderValueBase* value); 

int main(int const argc, char const * const * const argv){
  std::string outputfilename("");
  std::string inputfilename("");
  for(int i = 1; i < argc; ++i){
    if(std::string(argv[i]) == "-i") inputfilename = argv[i+1];
    if(std::string(argv[i]) == "-o") outputfilename = argv[i+1];
  }
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
  TFile *outputfile = new TFile(outputfilename.c_str(), "CREATE", "Track and Model Index");
  TTree *outputtree = new TTree("TrackModelIndex","Contains the tracking and model index for each particle");

  // Read a single float value in each tree entries:
  TTreeReaderValue< std::map< int, int > > trackIndex_modelIndex(reader, "Trajectory.trackIndex_modelIndex");
  //TTreeReaderValue< std::map< int, std::vector< int > > > modelIndex_trackIndex(reader, "Trajectory.modelIndex_trackIndex");
  // Make histogramms:
  int trackIndex(0), modelIndex(0);
  outputtree->Branch("trackIndex",&trackIndex,"trackIndex/I");
  outputtree->Branch("modelIndex",&modelIndex,"modelIndex/I");

  // Now iterate through the TTree entries and fill a histogram.
  while (reader.Next()) {
   if (Check_TTreeReader_EntryStatus(reader) == false) return -1; 

    for(auto i = trackIndex_modelIndex->begin(); i != trackIndex_modelIndex->end(); ++i){
      trackIndex = i->first;
      modelIndex = i->second;
      outputtree->Fill();
    }
  } // TTree entry / event loop
  outputtree->Write();
  outputfile->Close();
  inputfile->Close();
  return 0;
}

bool Check_TTreeReader_EntryStatus(TTreeReader const & redr){
  if (redr.GetEntryStatus() == TTreeReader::kEntryValid) {
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
//bool CheckValue(ROOT::TTreeReaderValueBase* value) {
  if (value->GetSetupStatus() < 0) {
    std::cerr << "Error " << value->GetSetupStatus()
      << "setting up reader for " << value->GetBranchName() << '\n';
    return false;
  }
  return true;
}
