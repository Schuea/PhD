#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"


#include <bitset>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "UsefulFunctions.h"
#include "GeneralFunctions_SiDBkgSim_new.h"
#include "CellHits_class_new.h"
#include "Subdetector_class_new.h"
#include "Style.h"

using namespace std;

double CalculatePhi(double x, double y);
uint64_t CalculateLayer(uint64_t const id, Subdetector const & SubDetector); 
uint32_t MakeNewCellID(double const x, double const y, Subdetector const & component);


int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  std::vector< std::string > *inputfilenames = new std::vector< std::string >();
  std::string argument_subdetectors;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
  std::string outputfile_name = oss.str();

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
	  else if (argv[i] == std::string("-s")) {
		  if (argv[i + 1] != NULL && 
				  argv[i + 1] != std::string("-i") && 
				  argv[i + 1] != std::string("-w") && 
				  argv[i + 1] != std::string("-o")) {
			  argument_subdetectors = argv[i + 1];
			  subdetector_set = true;
		  } else {
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
  if (!inputfile_set || !subdetector_set ) {
	  std::cerr
		  << "You didn't give the name for the subdector, the inputfiles or the weights. Please try again!"
		  << std::endl;
	  exit(1);
  }


  Subdetector det( argument_subdetectors );

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

  std::string subdetectorname = det.getName();
  TFile* Outputfile = new TFile(("output/Occupancy_"+subdetectorname+"_"+outputfile_name+"_tree.root").c_str(),"RECREATE");
  TTree* Outputtree = new TTree("CellHits","CellHits");

  int layer = 0;
  uint32_t HitCellID1 = 0;
  uint64_t combined_cell_id = 0;
  int HitCellID0(0);
  double HitPosition_x;
  double HitPosition_y;
  double HitPosition_z;
  int vector_element = 0;

  Outputtree->Branch("CellID0",&HitCellID0,"CellID0/I");
  Outputtree->Branch("CellID1",&HitCellID1,"CellID1/i");
  Outputtree->Branch("CellID",&combined_cell_id,"CellID/l");
  Outputtree->Branch("CellID_small",&vector_element,"CellID_small/I");
  Outputtree->Branch("Layer",&layer,"Layer/I");
  Outputtree->Branch("Pos_x",&HitPosition_x,"Pos_x/D");
  Outputtree->Branch("Pos_y",&HitPosition_y,"Pos_y/D");
  Outputtree->Branch("Pos_z",&HitPosition_z,"Pos_z/D");

  CellHits cellhits( &det );
  for (size_t file_iterator = 0; file_iterator < inputfilenames->size(); ++file_iterator) {
    //std::cout << inputfilenames->at(file_iterator) << std::endl;
    TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, det.getName());

    //Set the branches
    int pdg(0);
    double D_HitPosition_x(0.0), D_HitPosition_y(0.0), D_HitPosition_z(0.0); //double for Silicon detectors
    float F_HitPosition_x(0.0), F_HitPosition_y(0.0), F_HitPosition_z(0.0);

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("HitCellID0", 1);
    tree->SetBranchAddress("HitCellID0", &HitCellID0);

    tree->SetBranchStatus("HitPosition_x", 1);
    tree->SetBranchStatus("HitPosition_y", 1);
    tree->SetBranchStatus("HitPosition_z", 1);

    if (Silicon){
      tree->SetBranchStatus("HitParticle_PDG", 1);
      tree->SetBranchAddress("HitParticle_PDG", &pdg);

      tree->SetBranchAddress("HitPosition_x", &D_HitPosition_x);
      tree->SetBranchAddress("HitPosition_y", &D_HitPosition_y);
      tree->SetBranchAddress("HitPosition_z", &D_HitPosition_z);
    }
    else if (Calo){
      tree->SetBranchStatus("HitMotherParticle_PDG", 1);
      tree->SetBranchAddress("HitMotherParticle_PDG", &pdg);

      tree->SetBranchAddress("HitPosition_x", &F_HitPosition_x);
      tree->SetBranchAddress("HitPosition_y", &F_HitPosition_y);
      tree->SetBranchAddress("HitPosition_z", &F_HitPosition_z);
    }
    else{
      std::cerr << "This subdetector name was not recognized!" << std::endl; 
      exit(-1);
    }
 
    //Now we loop through the tree
    //Combine the two Cell ID's into a single new Cell ID
    //See how often the new Cell ID occurs in total, this is the occupancy

    long long int const entries = tree->GetEntries();
    std::cout << "There are " << entries << " entries" << std::endl;
    int skipped(0);
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      if (Silicon){
        HitPosition_x = D_HitPosition_x;
        HitPosition_y = D_HitPosition_y;
        HitPosition_z = D_HitPosition_z;
      }
      else if (Calo){
        HitPosition_x = F_HitPosition_x;
        HitPosition_y = F_HitPosition_y;
        HitPosition_z = F_HitPosition_z;
      }
      //if (endcap && HitPosition_z < 0){
      //  skipped++;
      //  continue;//Only compute one of the endcaps
      //}
      //Make a combined cell ID
      HitCellID1 = 0;
      combined_cell_id = 0;
      if (endcap){
        HitCellID1 = MakeNewCellID(HitPosition_x,HitPosition_y,det);
      }
      else if (barrel){
        double phi = CalculatePhi(HitPosition_x, HitPosition_y);
        //If Barrel, calculate CellID with z and phi on a certain radius
        if (Silicon) HitCellID1 = MakeNewCellID(HitPosition_z,phi*det.getRMin().at( CalculateLayer(HitCellID0,det)-1 ),det);//-1 for Silicon detectors only, because layer count starts from 1
        else if (Calo) HitCellID1 = MakeNewCellID(HitPosition_z,phi*det.getRMin().at( CalculateLayer(HitCellID0,det) ),det);
      }
      else{
        std::cerr << "This subdetector shape was not recognized!" << std::endl; 
        exit(-1);
      }
      combined_cell_id = (uint64_t(HitCellID1) << 32) | uint32_t(HitCellID0);
      //Use the CellHits class for storing the hit cells and their hitcounts
      vector_element = cellhits.Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y, HitPosition_z, det);

      layer = cellhits.CalculateLayer(combined_cell_id,det);
      Outputtree->Fill();
    }
    file->Close();
    std::cout << skipped << " events skipped" << std::endl;
  }
	Outputfile->Write();
  return 0;
}

uint64_t CalculateLayer(uint64_t const id, Subdetector const & SubDetector) {

  uint64_t LayerID64;

  LayerID64 = id << (64 - SubDetector.getLengthBitLayer() - SubDetector.getStartBitLayer());
  LayerID64 = LayerID64 >> (64 - SubDetector.getLengthBitLayer());

  return LayerID64;
}

double CalculatePhi(double x, double y){
  double phi =0;
  if(x>0) return phi = atan(y/x);
  if(x<0 && y>=0) return phi = atan(y/x) + M_PI;
  if(x<0 && y<0)  return phi = atan(y/x) - M_PI;
  if(x==0 && y>0)  return phi = 0.5*M_PI;
  if(x==0 && y<0)  return phi = -0.5*M_PI;
  else return -101;
}

uint32_t MakeNewCellID(double const x, double const y, Subdetector const & component){
  uint16_t newX = static_cast<uint16_t>(x/component.getCellSizeX());
  uint16_t newY = static_cast<uint16_t>(y/component.getCellSizeY());
  //std::cout << 3000/component.getCellSizeX() << std::endl;
  //std::cout << 3000/component.getCellSizeY() << std::endl;
  //std::cout << 3000/component.getCellSizeX()*3000/component.getCellSizeY() << std::endl;
  if(x >= 0) ++newX;
  if(y >= 0) ++newY;
  return (uint32_t(newX) << 16) | newY;
}
