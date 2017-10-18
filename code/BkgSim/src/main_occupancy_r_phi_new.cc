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
int Find_Cells_per_Bin(int layer, int bin, double binsize, int which_histo, Subdetector const & SubDetector);
void Draw_single_plots ( TH1D* histo, TCanvas* canvas, bool normalize, double normalization_factor, int integral_startbin, bool integral_numhits);
void Draw_multiple_plots (std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::vector< long long int > normalization_factor, int integral_startbin, bool integral_numhits);
void Print_multiple_plots_from_same_vec (std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::vector< long long int > normalization_factor, int integral_startbin, bool integral_numhits, std::string output);


int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  std::vector<std::string> *inputfilenames = new std::vector<std::string>();
  std::string argument_subdetectors;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
  std::string outputfile_name = oss.str();

  int bufferdepth = 4;
  int NUMBER_OF_FILES = 0;
  bool NUMBER_OF_FILES_set = false;
  bool inputfile_set = false;
  bool subdetector_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-o") && 
          argv[i + 1] != std::string("-b") && 
          argv[i + 1] != std::string("-s") && 
          argv[i + 1] != std::string("-i")) {
        NUMBER_OF_FILES = std::stoi(argv[i + 1]);
        std::cout << "Number of input files = " << NUMBER_OF_FILES << std::endl;
        NUMBER_OF_FILES_set = true;
      } else {
        std::cerr << "You didn't give an argument for the number of files!" << std::endl;
      }
    }
  }
  for (int i = 1; i < argc; i++) {
     if (argv[i] == std::string("-i") && 
        argv[i + 1] != std::string("-b") && 
        argv[i + 1] != std::string("-n") && 
        argv[i + 1] != std::string("-o") && 
        argv[i + 1] != std::string("-s")) {
      if (argv[i + 1] != NULL) {
        std::string filelist = argv[i + 1];
        if (access(filelist.c_str(), F_OK) == -1 ) {
          std::cerr << "The text file does not exist!" << std::endl;
          exit(2);
        }
        else {
          std::ifstream inputfilelist( filelist );
          int line = 1;
          std::string file;
          do {
            std::getline(inputfilelist, file);
            inputfilenames->push_back(file);
            ++line;
          } while (line <= NUMBER_OF_FILES);
          inputfile_set = true;
        }
      } else {
        std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
      }
    }
    if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-b") && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-o")) {
        argument_subdetectors = argv[i + 1];
        subdetector_set = true;
      } else {
        std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
      }
    }
    if (argv[i] == std::string("-b")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-s") && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-o")) {
        bufferdepth = std::stoi(argv[i + 1]);
      } else {
        std::cerr << "You didn't give an argument for the buffer depth!" << std::endl;
      }
    }
    if (argv[i] == std::string("-o")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-b") && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-s")) {
        outputfile_name = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the outputfile name!" << std::endl;
      }
    }
  }
  if (!inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set) {
    std::cerr
      << "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
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

  CellHits cellhits( &det );
  for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
    //std::cout << inputfilenames->at(file_iterator) << std::endl;
    TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, det.getName());

    //Set the branches
    int pdg(0);
    int HitCellID0(0);
    double D_HitPosition_x(0.0), D_HitPosition_y(0.0), D_HitPosition_z(0.0); //double for Silicon detectors
    float F_HitPosition_x(0.0), F_HitPosition_y(0.0), F_HitPosition_z(0.0);
    double HitPosition_x;
    double HitPosition_y;
    double HitPosition_z;

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
      if (endcap && HitPosition_z < 0) continue;//Only compute one of the endcaps
      
      //Make a combined cell ID
      uint32_t HitCellID1 = 0;
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
      uint64_t const combined_cell_id = (uint64_t(HitCellID1) << 32) | HitCellID0;
      //Use the CellHits class for storing the hit cells and their hitcounts
      cellhits.Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y, HitPosition_z, det);
    }
    file->Close();
  }

  std::string subdetectorname = det.getName();
  TFile* Outputfile = new TFile(("output/Occupancy_Phi_"+subdetectorname+"_"+outputfile_name+".root").c_str(),"RECREATE");
  //Make histogram vectors for storing the histograms
  std::vector < std::string > all_name;
  all_name.emplace_back( "All_layers_normocc1" );
  all_name.emplace_back( "All_layers_losthits1" );
  all_name.emplace_back( "All_layers_deadcells1" );
  all_name.emplace_back( "All_layers_normocc2" );
  all_name.emplace_back( "All_layers_losthits2" );
  all_name.emplace_back( "All_layers_deadcells2" );
  std::vector < std::string > title;
  title.emplace_back( "Normalized occupancy for " + subdetectorname + ";#phi;Average number of hits per cell" );
  std::stringstream title_temp;
  title_temp << "Number of hits lost for a buffer depth of " << bufferdepth << "for " << subdetectorname << ";#phi;Average number of hits lost";
  title.emplace_back( title_temp.str() );
  std::stringstream title_temp2;
  title_temp2 << "Number of dead cells for a buffer depth of " << bufferdepth << "for " << subdetectorname << ";#phi;Number of dead cells";
  title.emplace_back( title_temp2.str() );
  title.emplace_back( "Normalized occupancy for " + subdetectorname + ";z [mm];Average number of hits per cell" );
  std::stringstream title_temp3;
  title_temp3 << "Number of hits lost for a buffer depth of " << bufferdepth << "for " << subdetectorname << ";z [mm];Average number of hits lost";
  title.emplace_back( title_temp3.str() );
  std::stringstream title_temp4;
  title_temp4 << "Number of dead cells for a buffer depth of " << bufferdepth << "for " << subdetectorname << ";z [mm];Number of dead cells";
  title.emplace_back( title_temp4.str() );

  std::vector< TH1D* > normoccupancy_1;
  std::vector< TH1D* > normlosthits_1;
  std::vector< TH1D* > normdeadcells_1;
  std::vector< TH1D* > normoccupancy_2;
  std::vector< TH1D* > normlosthits_2;
  std::vector< TH1D* > normdeadcells_2;

  //std::cout << __LINE__ << std::endl;

  //Find out the maximum number of hits per cell and the total number of hits overall
  std::vector< long long int > tot_num_hits_per_layer;
  int tot_no_hits;
  int max_no_hits = 0;

  for (size_t vecpos = 0; vecpos < cellhits.Get_HitCount().size(); ++vecpos) {
    if (cellhits.Get_HitCount().at(vecpos) > max_no_hits){
      max_no_hits = cellhits.Get_HitCount().at(vecpos);
    }
    tot_no_hits += cellhits.Get_HitCount().at(vecpos);
  }
  int max_num_layers = det.getNumberOfLayers();

  //std::cout << __LINE__ << std::endl;
  int xbins1, xbins2 = 0;  
  double xmax1, xmin1 = 0.0;  
  double xmax2, xmin2 = 0.0;  
  if(barrel){
	  xmin1 = -1.0*M_PI;
	  xmax1 = M_PI;
	  xbins1 = 100;
	  xmin2 = -1.0*det.getZHalf().at(max_num_layers-1);
	  xmax2 = det.getZHalf().at(max_num_layers-1);
	  xbins2 = (int)det.getZHalf().at(max_num_layers-1);
  }
  else if(endcap){
	  xmin1 = -1.0*det.getRMax().at(max_num_layers-1);
	  xmax1 = det.getRMax().at(max_num_layers-1);
	  xbins1 = (int)det.getRMax().at(max_num_layers-1);
	  xmin2 = xmin1;
	  xmax2 = xmax1;
	  xbins2 = (int)det.getRMax().at(max_num_layers-1);
  }
  else std::cerr << "Detector shape not recognized. Histograms cannot be filled!" << std::endl;
  //std::cout << __LINE__ << std::endl;

  //Make the histograms
  TH1D* All_Layers_normoccupancy_1   = new TH1D(all_name.at(0).c_str(), title.at(0).c_str(), xbins1, xmin1, xmax1);
  TH1D* All_Layers_normlosthits_1    = new TH1D(all_name.at(1).c_str(), title.at(1).c_str(), xbins1, xmin1, xmax1);
  TH1D* All_Layers_normdeadcells_1   = new TH1D(all_name.at(2).c_str(), title.at(2).c_str(), xbins1, xmin1, xmax1);
  TH1D* All_Layers_normoccupancy_2   = new TH1D(all_name.at(3).c_str(), title.at(3).c_str(), xbins2, xmin2, xmax2);
  TH1D* All_Layers_normlosthits_2    = new TH1D(all_name.at(4).c_str(), title.at(4).c_str(), xbins2, xmin2, xmax2);
  TH1D* All_Layers_normdeadcells_2   = new TH1D(all_name.at(5).c_str(), title.at(5).c_str(), xbins2, xmin2, xmax2);

  if(endcap){
	All_Layers_normoccupancy_1 ->GetXaxis()->SetTitle("x [mm]");
  	All_Layers_normlosthits_1  ->GetXaxis()->SetTitle("x [mm]");
  	All_Layers_normdeadcells_1 ->GetXaxis()->SetTitle("x [mm]");
  	All_Layers_normoccupancy_2 ->GetXaxis()->SetTitle("y [mm]");
  	All_Layers_normlosthits_2  ->GetXaxis()->SetTitle("y [mm]");
  	All_Layers_normdeadcells_2 ->GetXaxis()->SetTitle("y [mm]");
  }
  //std::cout << __LINE__ << std::endl;
  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
    std::stringstream layername_occ1, layername_lost1, layername_dead1, layername_occ2, layername_lost2, layername_dead2;
    layername_occ1  << "Layer_" << number_layer << "_occupancy1";
    layername_lost1 << "Layer_" << number_layer << "_losthits1";
    layername_dead1 << "Layer_" << number_layer << "_deadcells1";
    layername_occ2  << "Layer_" << number_layer << "_occupancy2";
    layername_lost2 << "Layer_" << number_layer << "_losthits2";
    layername_dead2 << "Layer_" << number_layer << "_deadcells2";

    TH1D* temp1 = new TH1D(layername_occ1.str().c_str(),  title.at(0).c_str(), xbins1, xmin1, xmax1);
    TH1D* temp2 = new TH1D(layername_lost1.str().c_str(), title.at(1).c_str(), xbins1, xmin1, xmax1);
    TH1D* temp3 = new TH1D(layername_dead1.str().c_str(), title.at(2).c_str(), xbins1, xmin1, xmax1);
    TH1D* temp4 = new TH1D(layername_occ2.str().c_str(),  title.at(3).c_str(), xbins2, xmin2, xmax2);
    TH1D* temp5 = new TH1D(layername_lost2.str().c_str(), title.at(4).c_str(), xbins2, xmin2, xmax2);
    TH1D* temp6 = new TH1D(layername_dead2.str().c_str(), title.at(5).c_str(), xbins2, xmin2, xmax2);
    if(endcap){
       	 temp1->GetXaxis()->SetTitle("x [mm]");
    	 temp2->GetXaxis()->SetTitle("x [mm]");
    	 temp3->GetXaxis()->SetTitle("x [mm]");
    	 temp4->GetXaxis()->SetTitle("y [mm]");
    	 temp5->GetXaxis()->SetTitle("y [mm]");
    	 temp6->GetXaxis()->SetTitle("y [mm]");
    }
  //std::cout << __LINE__ << std::endl;
    normoccupancy_1.push_back(   temp1);
    normlosthits_1.push_back(	 temp2);
    normdeadcells_1.push_back(   temp3);
    normoccupancy_2.push_back(   temp4);
    normlosthits_2.push_back(	 temp5);
    normdeadcells_2.push_back(   temp6);
  }
  std::cout << __LINE__ << std::endl;

  //std::cout << __LINE__ << std::endl;
  std::vector< std::vector< int > > ave_hits_1;
  std::vector< std::vector< int > > ave_hits_2;
  std::vector< std::vector< int > > ave_losthits_1;
  std::vector< std::vector< int > > ave_losthits_2;
  std::vector< std::vector< int > > num_cells_per_bin_1;
  std::vector< std::vector< int > > num_cells_per_bin_2;
  long long int tot_num_cells = 0;
  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
	  tot_num_cells += det.getNumberOfCells().at(number_layer);
	  tot_num_hits_per_layer.push_back(0);

	//Prepare vectors that will contain all hits in certain bins of x
	  std::vector < int > temp1;
	  std::vector < int > temp2;
	  for (size_t bins = 0; bins < xbins1; ++bins) {
		  temp1.push_back( 0 );
	  }
	  for (size_t bins = 0; bins < xbins2; ++bins) {
		  temp2.push_back( 0 );
	  }
	  ave_hits_1.push_back( temp1 );
	  ave_hits_2.push_back( temp2 );
	  ave_losthits_1.push_back( temp1 );
	  ave_losthits_2.push_back( temp2 );
  	num_cells_per_bin_1.push_back( temp1 );
  	num_cells_per_bin_2.push_back( temp2 );
  }
  std::cout << __LINE__ << std::endl;
  double fill_x_1, fill_x_2 = 0.0;
  double binsize_1, binsize_2 = 0.0;
  for (size_t vecpos = 0; vecpos < cellhits.Get_HitCount().size(); ++vecpos) {
	  int current_layer = cellhits.Get_Layer().at(vecpos);
	  if (Silicon){
	          current_layer -= 1;//-1 for Silicon detectors only, because layer count starts from 1
	  }

	  if (barrel){
		  fill_x_1 = CalculatePhi( cellhits.Get_HitPosition('x').at(vecpos),cellhits.Get_HitPosition('y').at(vecpos) );
      binsize_1 = (2.0*M_PI)/(double)xbins1;
		  fill_x_2 = cellhits.Get_HitPosition('z').at(vecpos);
      binsize_2 = (2.0*det.getZHalf().at(max_num_layers - 1))/(double)xbins2;
	  }
	  else if (endcap){
		  fill_x_1 = cellhits.Get_HitPosition('x').at(vecpos);
      binsize_1 = (2.0*det.getRMax().at(max_num_layers - 1))/(double)xbins1;
		  fill_x_2 = cellhits.Get_HitPosition('y').at(vecpos);
      binsize_2 = (2.0*det.getRMax().at(max_num_layers - 1))/(double)xbins2;
	  }
	  else std::cerr << "Detector shape not recognized. Histograms cannot be filled!" << std::endl;
  std::cout << __LINE__ << std::endl;
	  std::cout << "fill_x_1 = " << fill_x_1 << std::endl;
	  std::cout << "fill_x_2 = " << fill_x_2 << std::endl;
	  int phi_bin_1, phi_bin_2 = 0;
    phi_bin_1 = (int)((fill_x_1+xmax1)/binsize_1); //Find in which bin fill_x_1 lies
	  std::cout << "phi_bin_1 = " << phi_bin_1 << std::endl;
    //int Find_Cells_per_Bin(int layer, int bin, double binsize, int which_histo, std::string subdetector_type);
  	num_cells_per_bin_1.at( current_layer ).at(phi_bin_1) = Find_Cells_per_Bin(current_layer, phi_bin_1, binsize_1, 1, det);
	  ave_hits_1.at( current_layer ).at(phi_bin_1) += cellhits.Get_HitCount().at(vecpos);
	  if(cellhits.Get_HitCount().at(vecpos) > bufferdepth) ave_losthits_1.at( current_layer ).at(phi_bin_1) += cellhits.Get_HitCount().at(vecpos) - bufferdepth;

	  phi_bin_2 = (int)((fill_x_2+xmax2)/binsize_2); //Find in which bin fill_x_2 lies
	  std::cout << "phi_bin_2 = " << phi_bin_2 << std::endl;
  	num_cells_per_bin_2.at( current_layer ).at(phi_bin_2) = Find_Cells_per_Bin(current_layer, phi_bin_2, binsize_2, 2, det);
	  ave_hits_2.at( current_layer ).at(phi_bin_2) += cellhits.Get_HitCount().at(vecpos);
	  if(cellhits.Get_HitCount().at(vecpos) > bufferdepth) ave_losthits_2.at( current_layer ).at(phi_bin_2) += cellhits.Get_HitCount().at(vecpos) - bufferdepth;
  }

  std::cout << __LINE__ << std::endl;

  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
	  for (size_t bins = 0; bins < xbins1; ++bins) {
	    	  normoccupancy_1.at(number_layer)->SetBinContent( bins, (double)ave_hits_1.at(number_layer).at(bins)/(double)num_cells_per_bin_1.at(number_layer).at(bins) );
		  All_Layers_normoccupancy_1->SetBinContent( bins, (double)ave_hits_1.at(number_layer).at(bins)/(double)num_cells_per_bin_1.at(number_layer).at(bins) );
	  }
	  for (size_t bins = 0; bins < xbins2; ++bins) {
	    	  normoccupancy_2.at(number_layer)->SetBinContent( bins, (double)ave_hits_2.at(number_layer).at(bins)/(double)num_cells_per_bin_2.at(number_layer).at(bins) );
		  All_Layers_normoccupancy_2->SetBinContent( bins, (double)ave_hits_2.at(number_layer).at(bins)/(double)num_cells_per_bin_2.at(number_layer).at(bins) );
	  }
  }
  for (size_t vecpos = 0; vecpos < cellhits.Get_HitCount().size(); ++vecpos) {
    if(cellhits.Get_HitCount().at(vecpos) > 0){
	    int current_layer = cellhits.Get_Layer().at(vecpos);
	    if (Silicon){
		    current_layer -= 1;//-1 for Silicon detectors only, because layer count starts from 1
	    }
      //std::cout << "cellhits.Get_HitCount().at(vecpos) = " << cellhits.Get_HitCount().at(vecpos) << std::endl;
	    if(cellhits.Get_HitCount().at(vecpos) >= bufferdepth){
		    normdeadcells_1.at(current_layer)->Fill( fill_x_1 );
		    normdeadcells_2.at(current_layer)->Fill( fill_x_2 );

		    All_Layers_normdeadcells_1->Fill( fill_x_1 );
		    All_Layers_normdeadcells_2->Fill( fill_x_2 );
	    }
	    tot_num_hits_per_layer.at(current_layer) += cellhits.Get_HitCount().at(vecpos);
    }
  }
  std::cout << __LINE__ << std::endl;


  //Plot the histogram and save it
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetLogy(1);

  std::stringstream output;
  output << "output/normoccupancy_1_" << subdetectorname << "_" << outputfile_name;
  Print_multiple_plots_from_same_vec (normoccupancy_1, canvas, true, tot_num_hits_per_layer , 1, false,  output.str());

  std::stringstream output2;
  output2 << "output/normoccupancy_2_" << subdetectorname << "_" << outputfile_name;
  Print_multiple_plots_from_same_vec (normoccupancy_2, canvas, true, tot_num_hits_per_layer, 1, false,  output2.str());

  std::stringstream output3;
  output3 << "output/normlosthits_1_" << subdetectorname << "_" << outputfile_name;
  Print_multiple_plots_from_same_vec (normlosthits_1, canvas, true, tot_num_hits_per_layer, 1, false,  output3.str());

  std::stringstream output4;
  output4 << "output/normlosthits_2_" << subdetectorname << "_" << outputfile_name;
  Print_multiple_plots_from_same_vec (normlosthits_2, canvas, true, tot_num_hits_per_layer, 1, false, output4.str());

  std::stringstream output5;
  output5 << "output/normdeadcells_1_" << subdetectorname << "_" << outputfile_name;
  Print_multiple_plots_from_same_vec (normdeadcells_1, canvas, true, det.getNumberOfCells(), 1, false,  output5.str());

  std::stringstream output6;
  output6 << "output/normdeadcells_2_" << subdetectorname << "_" << outputfile_name;
  Print_multiple_plots_from_same_vec (normdeadcells_2, canvas, true, det.getNumberOfCells(), 1, false, output6.str());

  //std::cout << __LINE__ << std::endl;

  Draw_single_plots( All_Layers_normoccupancy_1,canvas, true, tot_no_hits, 1, false); 
  std::stringstream All_output;
  All_output << "output/normoccupancy_1_all_layers_" << subdetectorname << "_" << outputfile_name;
  canvas->Print((All_output.str() + ".pdf").c_str());
  canvas->Print((All_output.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_normoccupancy_2,canvas, true, tot_no_hits, 1, false); 
  std::stringstream All_output2;
  All_output2 << "output/normoccupancy_2_all_layers_" << subdetectorname << "_" << outputfile_name;
  canvas->Print((All_output2.str() + ".pdf").c_str());
  canvas->Print((All_output2.str() + ".cxx").c_str());

  Draw_single_plots( All_Layers_normlosthits_1,canvas, true, tot_no_hits, 1, false); 
  std::stringstream All_output3;
  All_output3 << "output/normlosthits_1_all_layers_" << subdetectorname << "_" << outputfile_name;
  canvas->Print((All_output3.str() + ".pdf").c_str());
  canvas->Print((All_output3.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_normlosthits_2,canvas, true, tot_no_hits, 1, false); 
  std::stringstream All_output4;
  All_output4 << "output/normlosthits_2_all_layers_" << subdetectorname << "_" << outputfile_name;
  canvas->Print((All_output4.str() + ".pdf").c_str());
  canvas->Print((All_output4.str() + ".cxx").c_str());

  Draw_single_plots( All_Layers_normdeadcells_1,canvas, true, tot_num_cells, 1, false); 
  std::stringstream All_output5;
  All_output5 << "output/normdeadcells_1_all_layers_" << subdetectorname << "_" << outputfile_name;
  canvas->Print((All_output5.str() + ".pdf").c_str());
  canvas->Print((All_output5.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_normdeadcells_2,canvas, true, tot_num_cells, 1, false); 
  std::stringstream All_output6;
  All_output6 << "output/normdeadcells_2_all_layers_" << subdetectorname << "_" << outputfile_name;
  canvas->Print((All_output6.str() + ".pdf").c_str());
  canvas->Print((All_output6.str() + ".cxx").c_str());



  
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
  phi = atan2(y, x);
  //if(x>0) phi = atan(y/x);
  //if(x<0 && y>=0) phi = atan(y/x) + M_PI;
  //if(x<0 && y<0)  phi = atan(y/x) - M_PI;
  //if(x==0 && y>0) phi = 0.5*M_PI;
  //if(x==0 && y<0) phi = -0.5*M_PI;
  //else phi = -101;
  //phi += M_PI;
  return phi;
}

uint32_t MakeNewCellID(double const x, double const y, Subdetector const & component){
  uint16_t newX = static_cast<uint16_t>(x/component.getCellSizeX());
  uint16_t newY = static_cast<uint16_t>(y/component.getCellSizeY());
  if(x >= 0) ++newX;
  if(y >= 0) ++newY;
  return (uint32_t(newX) << 16) | newY;
}
void Draw_single_plots ( TH1D* histo, TCanvas* canvas, bool normalize, double normalization_factor, int integral_startbin, bool integral_numhits){
  int tot = 0;
  for(int bin = integral_startbin; bin <= histo->GetNbinsX(); ++bin){
    if (integral_numhits == true) tot += histo->GetBinContent(bin)*histo->GetBinLowEdge(bin);
    else tot += histo->GetBinContent(bin);
  }
  
  histo->SetStats(0);
  histo->Sumw2(1);
  histo->SetMinimum(0.1);
  if(normalize == true){
    if(normalization_factor == 0.0) histo->Scale(1.0/histo->GetBinContent(1));
    else histo->Scale(1.0/normalization_factor);
    histo->SetMinimum( pow(10,-12) );
  }
  histo->SetLineColor(2);
  histo->Draw("hist");
  canvas->Update();

  TPaveText *text1 = new TPaveText(0.75,0.8,0.95,0.9,"brNDC");
  text1->SetTextFont(62);
  text1->SetTextColor(2);
  text1->SetFillColor(0);
  std::stringstream title;
  title << histo->GetName();
  text1->AddText(title.str().c_str());
  text1->AddLine(0,0.5,1,0.5);
  std::stringstream entries;
  entries << "Entries = " << tot; 
  text1->AddText(entries.str().c_str());
  text1->Draw();
}

int Find_Cells_per_Bin(int layer, int bin, double binsize, int which_histo, Subdetector const & det){
  int num_cells = 0;
  if (det.getShape().find("barrel") != std::string::npos) {
    if (which_histo == 1) num_cells = (det.getZHalf().at(layer)*2.0/det.getCellSizeX()) * (binsize*det.getRMax().at(layer)/det.getCellSizeX());//number of cells along barrel * number of cells on part of circumference
    else if (which_histo == 2) num_cells = (binsize/det.getCellSizeX()) * (2.0*M_PI*det.getRMax().at(layer)/det.getCellSizeX());//number of cells in z bin * number of cells along whole circumference
    else{
      std::cerr << "This histo config does not exist!" << std::endl; 
      exit(-2);
    }
  }
  else if (det.getShape().find("endcap") != std::string::npos) {
    if (which_histo == 1){
      if (std::abs(bin) <= det.getRMin().at(layer)){
        num_cells = (2.0*(det.getRMax().at(layer)*sin( acos( std::abs(bin)/det.getRMax().at(layer) ) )) - 2.0*(det.getRMin().at(layer)*sin( acos( std::abs(bin)/det.getRMin().at(layer) ) )))/det.getCellSizeY();
      }
      else num_cells = 2.0*(det.getRMax().at(layer)*sin( acos( std::abs(bin)/det.getRMax().at(layer) ) ))/det.getCellSizeY();
    }
    else if (which_histo == 2){
      if (std::abs(bin) <= det.getRMin().at(layer)){
        num_cells = (2.0*(det.getRMax().at(layer)*sin( acos( std::abs(bin)/det.getRMax().at(layer) ) )) - 2.0*(det.getRMin().at(layer)*sin( acos( std::abs(bin)/det.getRMin().at(layer) ) )))/det.getCellSizeX();
      }
      else num_cells = 2.0*(det.getRMax().at(layer)*sin( acos( std::abs(bin)/det.getRMax().at(layer) ) ))/det.getCellSizeX();
    }
    else{
      std::cerr << "This histo config does not exist!" << std::endl; 
      exit(-2);
    }
  }
  else std::cerr << "Detector shape not recognized. Histograms cannot be filled!" << std::endl;
  return num_cells;
}

void Print_multiple_plots_from_same_vec (std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::vector< long long int > normalization_factor, int integral_startbin, bool integral_numhits, std::string output){
  Draw_multiple_plots(histos, canvas, normalize, normalization_factor, integral_startbin, integral_numhits);
  canvas->Print((output + ".pdf").c_str());
  canvas->Print((output + ".cxx").c_str());
}

void Draw_multiple_plots (std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::vector< long long int > normalization_factor, int integral_startbin, bool integral_numhits){
  int i = 0;
  int color = 2; // Very first histogram will be drawn with the color 2, then counted up
  int marker = 20; // Very first histogram will be drawn with the marker style 20, then counted up
  for (size_t number_histo = 0; number_histo< histos.size(); ++number_histo) {
    histos.at(number_histo)->SetStats(0);
    histos.at(number_histo)->Sumw2(1);
  }
  double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
  for (size_t number_histo = 0; number_histo< histos.size(); ++number_histo) {
    if(number_histo == 0){
      int tot = 0;
      for(int bin = integral_startbin; bin <= histos.at(number_histo)->GetNbinsX(); ++bin){
        if (integral_numhits == true) tot += histos.at(number_histo)->GetBinContent(bin)*histos.at(number_histo)->GetBinLowEdge(bin);
        else tot += histos.at(number_histo)->GetBinContent(bin);
      }

      if(normalize == true){
        if( normalization_factor.empty() ) histos.at(number_histo)->Scale(1.0/histos.at(number_histo)->GetBinContent(1));
        else histos.at(number_histo)->Scale( 1.0/normalization_factor.at(number_histo) );
        histos.at(number_histo)->SetMinimum( pow(10,-12) );
      }
      else{
        histos.at(number_histo)->SetMaximum(max);
        histos.at(number_histo)->SetMinimum(0.1);
      }
      histos.at(number_histo)->SetLineColor(color);
      histos.at(number_histo)->SetMarkerColor(color);
      histos.at(number_histo)->SetMarkerStyle(marker);
      histos.at(number_histo)->Draw("P");
      canvas->Update();
      TPaveText *text1;
      if (histos.size() > 5){
        text1 = new TPaveText(0.6,0.8,0.8,0.9,"brNDC");
      }
      else{
        text1 = new TPaveText(0.75,0.8,0.95,0.9,"brNDC");
      }
      text1->SetTextFont(62);
      text1->SetTextColor(color);
      text1->SetFillColor(0);
      std::stringstream title;
      title << histos.at(number_histo)->GetName();
      text1->AddText(title.str().c_str());
      text1->AddLine(0,0.5,1,0.5);
      std::stringstream entries1;
      entries1 << "Entries = " << tot; 
      text1->AddText(entries1.str().c_str());
      text1->Draw();
      i=1;
    }
    if(number_histo > 0){
      int tot = 0;
      for(int bin = integral_startbin; bin <= histos.at(number_histo)->GetNbinsX(); ++bin){
        if (integral_numhits == true) tot += histos.at(number_histo)->GetBinContent(bin)*histos.at(number_histo)->GetBinLowEdge(bin);
        else tot += histos.at(number_histo)->GetBinContent(bin);
      }

      if(normalize == true){
        if( normalization_factor.empty() ) histos.at(number_histo)->Scale(1.0/histos.at(number_histo)->GetBinContent(1));
        else histos.at(number_histo)->Scale( 1.0/normalization_factor.at(number_histo) );
        histos.at(number_histo)->SetMinimum( pow(10,-12) );
      }
      else{
        histos.at(number_histo)->SetMaximum(max);
        histos.at(number_histo)->SetMinimum(0.1);
      }
      color++;
      marker++;
      if(color == 5 || color == 10) color += 1; // 5 would be yellow, 10 would be very light gray 
      histos.at(number_histo)->SetLineColor(color);
      histos.at(number_histo)->SetMarkerColor(color);
      histos.at(number_histo)->SetMarkerStyle(marker);
      histos.at(number_histo)->Draw("P,SAMES");
      canvas->Update();
      TPaveText *text2;
      if (histos.size() > 5) {
        if(number_histo >= 5){
          if(number_histo == 5) {
            i=0;
          }
          text2 = new TPaveText(0.75,0.8-i*0.1,0.95,0.9-i*0.1,"brNDC");
        }
        else {
          text2 = new TPaveText(0.6,0.8-i*0.1,0.8,0.9-i*0.1,"brNDC");
        }
      } 
      else{
        text2 = new TPaveText(0.75,0.8-i*0.1,0.95,0.9-i*0.1,"brNDC");
      }
      text2->SetTextFont(62);
      text2->SetTextColor(color);
      text2->SetFillColor(0);
      std::stringstream title2;
      title2 << histos.at(number_histo)->GetName();
      text2->AddText(title2.str().c_str());
      text2->AddLine(0,0.5,1,0.5);
      std::stringstream entries2;
      entries2 << "Entries = " << tot;//-1 because of entries in first bin count 1 because of setbincontent
      text2->AddText(entries2.str().c_str());
      text2->Draw();
      i++;
    }
  }
}

