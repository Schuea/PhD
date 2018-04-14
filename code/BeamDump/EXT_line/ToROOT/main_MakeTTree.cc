#include "TFile.h"
#include "TH1F.h"
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

	TFile* ROOTFile = new TFile(outputfilename.c_str(),"RECREATE","EXT_neutrons");
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
  
  TH1F* X_histo  = new TH1F("X_histo","X_histo",100,-2000,6000);
  TH1F* Y_histo  = new TH1F("Y_histo","Y_histo",100,-2000,2000);
  TH1F* CX_histo  = new TH1F("CX_histo","CX_histo",100,-1,1);
  TH1F* CXlow_histo  = new TH1F("CXlow_histo","CXlow_histo",100,-0.05,0.05);
  TH1F* CX_histo_ex  = new TH1F("CX_histo_ex","CX_histo_ex",100,-1,1);
  TH1F* CY_histo  = new TH1F("CY_histo","CY_histo",100,-1,1);
  TH1F* CY_histo_ex  = new TH1F("CY_histo_ex","CY_histo_ex",100,-1,1);
  TH1F* CYlow_histo  = new TH1F("CYlow_histo","CYlow_histo",100,-0.05,0.05);
  TH1F* CZ_histo  = new TH1F("CZ_histo","CZ_histo",100,-1,1);
  TH1F* CZlow_histo  = new TH1F("CZlow_histo","CZlow_histo",100,-1,-0.9995);
  TH1F* E_histo  = new TH1F("E_histo","E_histo",100,0,0.06);
  TH1F* Elow_histo  = new TH1F("Elow_histo","Elow_histo",100,0,0.005);
  TH1F* Elowlow_histo  = new TH1F("Elowlow_histo","Elowlow_histo",100,0,0.00000005);
  TH1F* T_histo  = new TH1F("T_histo","T_histo",100,0,0.012);
  TH1F* Tlow_histo  = new TH1F("Tlow_histo","Tlow_histo",100,0,0.00003);

  std::ifstream inputfile(inputfilename);
  std::string line;
  int count = 1;
  //Now start reading in:
  while (!inputfile.eof()){
    std::getline(inputfile, line);
    std::istringstream in(line);
    std::string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11;

    if (count == 1){
      std::cout << col1 << ", " << col2 << ", " << col3 << ", " << col4 << ", " << col5 << ", " << col6 << ", " << col7 << ", " << col8 << ", " << col9 << std::endl;
    }
    //The positions are given in cm, time in s!
    //Therefore transform units to nicer ones whilst storing in TTree:
    x = std::atof(col1.c_str())*10.0;//unit: mm
    y = std::atof(col2.c_str())*10.0;//unit: mm
    z = std::atof(col3.c_str())*10.0;//unit: mm
    cx = std::atof(col4.c_str());
    cy = std::atof(col5.c_str());
    cz = std::atof(col6.c_str());
    ekin = std::atof(col7.c_str());//unit: GeV
    time = std::atof(col8.c_str())*std::pow(1.0,9.0);//unit: ns
    weight = std::atof(col9.c_str());

    Tree->Fill();
    ++count;

    X_histo->Fill(x, weight);
    Y_histo->Fill(y, weight);
    CX_histo->Fill(cx, weight);
    if(cx < -0.01 || cx > 0.022) CX_histo_ex->Fill(cx, weight);
    CXlow_histo->Fill(cx, weight);
    CY_histo->Fill(cy, weight);
    if(cy > 0.006 || cy < -0.006) CY_histo_ex->Fill(cy, weight);
    CYlow_histo->Fill(cy, weight);
    CZ_histo->Fill(cz, weight);
    CZlow_histo->Fill(cz, weight);
    E_histo->Fill(ekin, weight);
    Elow_histo->Fill(ekin, weight);
    Elowlow_histo->Fill(ekin, weight);
    T_histo->Fill(time, weight);
    Tlow_histo->Fill(time, weight);
  }

  inputfile.close();
  ROOTFile->Write();
  ROOTFile->Close();

  return 0;
}
