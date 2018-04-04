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

  float x = 0.0;
  float y = 0.0;
  float z = 0.0;
  float cx = 0.0;
  float cy = 0.0;
  float cz = 0.0;
  float ekin = 0.0;
  float time = 0.0;
  float weight = 0.0;

	TFile* ROOTFile = new TFile(outputfilename.c_str(),"CREATE","EXT_neutrons");
  TTree* Tree = new TTree("Tree_Neutrons","TTree for EXT neutrons");
  
  Tree->Branch("X",&x,"X/F");
  Tree->Branch("Y",&y,"Y/F");
  Tree->Branch("Z",&z,"Z/F");
  Tree->Branch("CX",&cx,"CX/F");
  Tree->Branch("CY",&cy,"CY/F");
  Tree->Branch("CZ",&cz,"CZ/F");
  Tree->Branch("Ekin",&ekin,"Ekin/F");
  Tree->Branch("Time",&time,"Time/F");
  Tree->Branch("Weight",&weight,"Weight/F");
  
  std::ifstream inputfile(inputfilename);
  std::string line;
  //Now start reading in:
  while (!inputfile.eof()){
    std::getline(inputfile, line);
    std::istringstream in(line);
    std::string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11;

    //The positions are given in cm, time in s!
    //Therefore transform units to nicer ones whilst storing in TTree:
    x = std::atoi(col1.c_str())*10.0;//unit: mm
    y = std::atoi(col2.c_str())*10.0;//unit: mm
    z = std::atoi(col3.c_str())*10.0;//unit: mm
    cx = std::atoi(col4.c_str());
    cy = std::atoi(col5.c_str());
    cz = std::atoi(col6.c_str());
    ekin = std::atoi(col7.c_str());//unit: GeV
    time = std::atoi(col8.c_str())*std::pow(1.0,9.0);//unit: ns
    weight = std::atoi(col9.c_str());

    Tree->Fill();
  }
  inputfile.close();
  ROOTFile->Write();
  ROOTFile->Close();

  return 0;
}
