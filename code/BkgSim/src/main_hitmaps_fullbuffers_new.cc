#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TPaveStats.h"


#include <bitset>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
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
int CalculateLayer(long long int const id, Subdetector const & SubDetector); 
long long int MakeNewCellID(double const x, double const y, Subdetector const & component);
std::pair< double, double > Get_historange(int layer, Subdetector detector);
std::pair< double, double > Get_historangebins(int layer, Subdetector detector);
void Draw_single_plots ( TH2D* histo );
void Draw_single_plots ( TH3D* histo );

double weight = 1.0;

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  std::vector<std::string> *inputfilenames = new std::vector<std::string>();
  std::string argument_subdetectors;
  int bufferdepth = 4;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
  std::string outputfile_name = oss.str();

  int NUMBER_OF_FILES = 0;
  bool NUMBER_OF_FILES_set = false;
  bool inputfile_set = false;
  bool subdetector_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-b") && 
          argv[i + 1] != std::string("-o") && 
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
        argv[i + 1] != std::string("-n") && 
        argv[i + 1] != std::string("-o") && 
        argv[i + 1] != std::string("-b") && 
        argv[i + 1] != std::string("-s")) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-n") && 
              argv[i + j] != std::string("-o") &&
              argv[i + j] != std::string("-b") &&
              argv[i + j] != std::string("-s")) {
            inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_FILES);
        inputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
      }
    }
    if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-b") && 
          argv[i + 1] != std::string("-o")) {
        argument_subdetectors = argv[i + 1];
        subdetector_set = true;
      } else {
        std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
      }
    }
    if (argv[i] == std::string("-o")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-b") && 
          argv[i + 1] != std::string("-s")) {
        outputfile_name = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the outputfile name!" << std::endl;
      }
    }
    if (argv[i] == std::string("-b")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-i") && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-o") && 
          argv[i + 1] != std::string("-s")) {
        bufferdepth = std::stoi( argv[i + 1] );
      } else {
        std::cerr << "You didn't give an argument for the bufferdepth!" << std::endl;
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
      //if (pdg != 13 && pdg != -13) continue;
      //if (endcap && *HitPosition_z < 0) continue;
      //std::cout << "HitPosition_x = " << HitPosition_x << std::endl;
      //std::cout << "HitPosition_y = " << HitPosition_y << std::endl;
      //std::cout << "HitPosition_z = " << HitPosition_z << std::endl;
      //Make a combined cell ID
      long long int HitCellID1 = 0;
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
      //std::cout << "HitCellID1 = " << HitCellID1 << std::endl;
      long long int const combined_cell_id = (long long) HitCellID1 << 32 | HitCellID0;
      //Use the CellHits class for storing the hit cells and their hitcounts
      cellhits.Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y, HitPosition_z);
    }
    file->Close();
  }

  std::string subdetectorname = det.getName();
	TFile* Outputfile = new TFile(("output/OccupancyMap_"+subdetectorname+"_"+outputfile_name+".root").c_str(),"RECREATE");
  //Make histogram vectors for storing the histograms
  std::vector < std::string > title;
  if (barrel){
    title.emplace_back( "Occupancy map for " + subdetectorname + ";z [mm];phi [rad];Number of hits per cell" );
    title.emplace_back( "Number of dead cells for" + subdetectorname + ";z [mm];phi [rad];Number of dead cells" );
  }
  else {
    title.emplace_back( "Occupancy map for " + subdetectorname + ";x [mm];y [mm];Number of hits per cell" );
    title.emplace_back( "Number of dead cells for" + subdetectorname + ";x [mm];y [mm];Number of dead cells" );
  }
  title.emplace_back( "Hitpositions for" + subdetectorname + ";x [mm];y [mm];z [mm];Number of hits" );

  std::vector< TH2D* > histos;
  std::vector< TH2D* > histos_deadcells;
  std::vector< TH3D* > histos3D;

  //Find out the maximum number of hits per cell and the total number of hits overall
  int tot_no_hits;
  int max_no_hits = 0;
  for (size_t vecpos = 0; vecpos < cellhits.Get_HitCount().size(); ++vecpos) {
    if (cellhits.Get_HitCount().at(vecpos) > max_no_hits){
      max_no_hits = cellhits.Get_HitCount().at(vecpos);
    }
    tot_no_hits += cellhits.Get_HitCount().at(vecpos);
  }


  //Make the histograms
  int max_num_layers = det.getNumberOfLayers();
  double xrange = 1.0;
  double xrangebins = 1.0;
  double yrange = 1.0;
  double yrangebins = 1.0;

  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
    std::stringstream layername, layername2, layername3;
    layername << "Layer_" << number_layer;
    layername2 << "Layer_" << number_layer << "_deadcells";
    layername3 << "Layer_" << number_layer << "_3D";
    if (Silicon){
      xrange = Get_historange(number_layer, det).first;
      yrange = Get_historange(number_layer, det).second;
      xrangebins = Get_historangebins(number_layer, det).first;
      yrangebins = Get_historangebins(number_layer, det).second;
    }
    else if (Calo) {//+1 for Calorimeters only, because their layer count starts from 1
      xrange = Get_historange(number_layer+1, det).first;
      yrange = Get_historange(number_layer+1, det).second;
      xrangebins = Get_historangebins(number_layer+1, det).first;
      yrangebins = Get_historangebins(number_layer+1, det).second;
    }

    TH2D* temp1 = new TH2D(layername.str().c_str(),  title.at(0).c_str(), 300, -xrange, xrange, 200, -yrange, yrange);
    TH2D* temp2 = new TH2D(layername2.str().c_str(), title.at(1).c_str(), 300, -xrange, xrange, 200, -yrange, yrange);
    //TH2D* temp1 = new TH2D(layername.str().c_str(),  title.at(0).c_str(), xrangebins, -xrange, xrange, yrangebins, -yrange, yrange);
    //TH2D* temp2 = new TH2D(layername2.str().c_str(), title.at(1).c_str(), xrangebins, -xrange, xrange, yrangebins, -yrange, yrange);
    histos.push_back( temp1 );
    histos_deadcells.push_back( temp2 );
    TH3D* temp3 = new TH3D(layername3.str().c_str(), title.at(2).c_str(), 300, -det.getRLayer().at(number_layer)*1.1, det.getRLayer().at(number_layer)*1.1, 200, -det.getRLayer().at(number_layer)*1.1, det.getRLayer().at(number_layer)*1.1, 200, -det.getZHalf().at(number_layer)*1.1, det.getZHalf().at(number_layer)*1.1 );
    histos3D.push_back( temp3 );
  }
  //Filling the primary histograms with the entries from the cellhits
  double xvalue = 0.0;
  double yvalue = 0.0;
  for (size_t vecpos = 0; vecpos < cellhits.Get_HitCount().size(); ++vecpos) {
    if(cellhits.Get_HitCount().at(vecpos) > 0){
      //std::cout << "Layer: " << cellhits.Get_Layer().at(vecpos) << std::endl;
      int current_layer = cellhits.Get_Layer().at(vecpos);
      if (barrel){
        xvalue = cellhits.Get_HitPosition('z').at(vecpos);
        yvalue = CalculatePhi( cellhits.Get_HitPosition('x').at(vecpos), cellhits.Get_HitPosition('y').at(vecpos));
      }
      else {
        xvalue = cellhits.Get_HitPosition('x').at(vecpos);
        yvalue = cellhits.Get_HitPosition('y').at(vecpos);
      }
      if (Silicon){
        histos.at(current_layer -1 ) -> Fill( xvalue, yvalue, cellhits.Get_HitCount().at(vecpos) );//-1 for Silicon detectors only, because layer count starts from 1
        histos3D.at(current_layer -1 ) -> Fill( cellhits.Get_HitPosition('x').at(vecpos), cellhits.Get_HitPosition('y').at(vecpos), cellhits.Get_HitPosition('z').at(vecpos) );//-1 for Silicon detectors only, because layer count starts from 1
      }
      else{
        histos.at(current_layer) -> Fill( xvalue, yvalue, cellhits.Get_HitCount().at(vecpos) );
        histos3D.at(current_layer) -> Fill( cellhits.Get_HitPosition('x').at(vecpos), cellhits.Get_HitPosition('y').at(vecpos), cellhits.Get_HitPosition('z').at(vecpos) );//-1 for Silicon detectors only, because layer count starts from 1
      }
    }
  }

  //Filling bufferdepth plots:
  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
      for (int bin = 1; bin < histos.at(number_layer)->GetNbinsX()*histos.at(number_layer)->GetNbinsY(); ++bin) {//go through the previous histos (showing the hit cells)
        if (histos.at(number_layer)->GetBinContent(bin) > bufferdepth){ 
          histos_deadcells.at(number_layer)->SetBinContent(bin, 1);
        }
        else{
          histos_deadcells.at(number_layer)->SetBinContent(bin, 0);
        }
      }
  }

  std::cout<< "---------------" <<std::endl;
  std::cout<< "Total number of hits counted for subdetector "<< subdetectorname <<std::endl;
  std::cout << tot_no_hits <<std::endl;
  std::cout<< "---------------" <<std::endl;

  //Plot the histogram and save it
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);

  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
    Draw_single_plots( histos.at(number_layer) ); 
    std::stringstream histos_output;
    histos_output << "output/OccupancyMap_" << subdetectorname << "_layer_" << number_layer << "_" << outputfile_name;
    canvas->Update();
    canvas->Print((histos_output.str() + ".pdf").c_str());
    canvas->Print((histos_output.str() + ".cxx").c_str());

    Draw_single_plots ( histos_deadcells.at(number_layer) );
    std::stringstream histos_deadcells_output;
    histos_deadcells_output << "output/OccupancyMap_deadcells_" << subdetectorname << "_layer_" << number_layer << "_" << outputfile_name;
    canvas->Update();
    canvas->Print((histos_deadcells_output.str() + ".pdf").c_str());
    canvas->Print((histos_deadcells_output.str() + ".cxx").c_str());
  }
  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
  for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetPalette(1);
    TPad *th3Pad  = new TPad("box", "box", 0.01, 0.01, 0.99, 0.84);
    th3Pad->SetTopMargin(0);
    th3Pad->SetBottomMargin(0);
    th3Pad->Draw();
    th3Pad->cd();
    //histo->SetContour(Max);
    histos3D.at(number_layer)->SetDirectory(0);
    //histo->SetStats(1);
    histos3D.at(number_layer)->Draw("glcolz");
    histos3D.at(number_layer)->Write();
    //Draw_single_plots( histos3D.at(number_layer) ); 
    std::stringstream histos3D_output;
    histos3D_output << "output/OccupancyMap_3D_" << subdetectorname << "_layer_" << number_layer << "_" << outputfile_name;
    canvas2->Write();
    canvas2->SaveAs((histos3D_output.str() + ".png").c_str());
    canvas2->SaveAs((histos3D_output.str() + ".pdf").c_str());
    canvas2->SaveAs((histos3D_output.str() + ".cxx").c_str());
  }

	Outputfile->Write();
  return 0;
}

int CalculateLayer(long long int const id, Subdetector const & SubDetector) {
  std::bitset<64> cellidbit(id);
  std::string CellID_ = cellidbit.to_string();
  int LayerInt = -1;
  std::stringstream LayerID;

  //This for loop calculates the layer id
  //From a sring of 0's and 1's, e.g. 00001011010010
  //The StartBin is the first bin in the string we are interested in (when reading from right to left)
  //The LengthBin is the length of the string we are interested in
  //We read from left to right, but we specify the start position from right to left
  //There is a magic +1 in there because strings start at element 0.
  for (int i = CellID_.size() - (SubDetector.getStartBitLayer() + SubDetector.getLengthBitLayer()); i <= CellID_.size() - (SubDetector.getStartBitLayer() + 1);
      ++i) {
    LayerID << CellID_.at(i);
  }

  std::bitset<64> LayerIDbit(LayerID.str());
  LayerInt = LayerIDbit.to_ulong();

  return LayerInt;
}

double CalculatePhi(double x, double y){
  double phi = 0.0;
  phi = atan2(y,x);
  //std::cout << "x = " << x << ", y = " << y << " -> phi = " << phi << std::endl;
  return phi;
  //if(x>0) return phi = atan(y/x);
  //if(x<0 && y>=0) return phi = atan(y/x) + M_PI;
  //if(x<0 && y<0)  return phi = atan(y/x) - M_PI;
  //if(x==0 && y>0)  return phi = 0.5*M_PI;
  //if(x==0 && y<0)  return phi = -0.5*M_PI;
  //else return -101;
}

long long int MakeNewCellID(double const x, double const y, Subdetector const & component){
  int newX = static_cast<int>(x/component.getCellSizeX()); //Check if Cell Size Area is the same as Cell Dimension
  int newY = static_cast<int>(y/component.getCellSizeY()); //Check if Cell Size Area is the same as Cell Dimension
  if(x >= 0) ++newX;
  if(y >= 0) ++newY;
  bitset<32> bitY = newY;
  newY = 0;
  for(int i = 0; i < 31; ++i){
    newY += bitY[i]*pow(2,i);
  }
  return (long long int) newX << 32 | newY;
}

std::pair< double, double > Get_historange(int layer, Subdetector detector){
  double xrange = 0.0;
  double yrange = 0.0;
  if (detector.getShape().find("barrel") != std::string::npos){
    xrange = detector.getZHalf().at(layer);
    yrange = M_PI;
  }
  else if (detector.getShape().find("endcap") != std::string::npos) {
    xrange = detector.getRLayer().at(layer);
    yrange = detector.getRLayer().at(layer);
  }
  else {
    std::cerr << "The given subdetector shape was not recognized!" << std::endl; 
    exit(-1);
  }
  std::pair< double, double > pair (xrange, yrange);
  return pair;
}

std::pair< double, double > Get_historangebins(int layer, Subdetector detector){
  double xrangebins = 0.0;
  double yrangebins = 0.0;
  if (detector.getShape().find("barrel") != std::string::npos){
    xrangebins = int( 2*detector.getZHalf().at(layer)/detector.getCellSizeX() );
    yrangebins = int( (2*M_PI*detector.getRLayer().at(layer))/detector.getCellSizeY() );
  }
  else if (detector.getShape().find("endcap") != std::string::npos) {
    xrangebins = int( 2*detector.getRLayer().at(layer)/detector.getCellSizeX() );
    yrangebins = int( 2*detector.getRLayer().at(layer)/detector.getCellSizeY() );
  }
  else {
    std::cerr << "The given subdetector shape was not recognized!" << std::endl; 
    exit(-1);
  }
  std::pair< double, double > pair (xrangebins, yrangebins);
  return pair;
}

void Draw_single_plots ( TH2D* histo){
  gStyle->SetCanvasPreferGL(kFALSE);
  histo->SetStats(1);
  //histo->SetMinimum(0.1);
  histo->Draw("colz");
}
void Draw_single_plots ( TH3D* histo){
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetPalette(1);
  //histo->SetContour(Max);
  histo->SetDirectory(0);
  histo->SetStats(1);
  //histo->SetMinimum(0.1);
  histo->Draw("glcolz");
  histo->Write();
}
