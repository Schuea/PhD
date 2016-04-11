#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"

#include "Style.h"
#include "UsefulFunctions.h"

#include <vector>
#include <map>
#include <array>
#include <iostream>
#include <sstream>

int const n = 10; //Number of apertures that were recorded
float HistoMax = 0.0; //To find the maximum entry to the histogram, so that the y-axis can be scaled appropriately.
float HistoMin = 1000000.0; //To find the minimum entry to the histogram, so that the y-axis can be scaled appropriately.

void GetAverageSignals(float* SignalAverage, bool GetError, const float* beamintensity, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* firstjawposition_branch, const float* secondjawposition_branch, const float* signal_branch, const int* voltage_branch);

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  //The input is a TTree ROOT file(s)
  //The output is .pdf and .C files

  std::string inputfilename;
  std::string outputfilename;
  float recorded_beamIntensity;
  float around_halfaperture = 0.0;

  bool around_halfaperture_set = false;
  bool inputfile_set = false;
  bool outputfile_set = false;
  bool beamintensity_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")){
      if (argv[i + 1] != std::string("-o")
          && argv[i + 1] != std::string("-a")
          && argv[i + 1] != std::string("-b")
          && argv[i + 1] != NULL) {
        inputfilename = argv[i + 1];
        inputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the input filename!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-o")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-a")
          && argv[i + 1] != std::string("-b")
          && argv[i + 1] != std::string("-i")) {
        outputfilename = argv[i + 1];
        outputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the output filename!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-b")) {
      if (argv[i + 1] != NULL
          && argv[i + 1] != std::string("-a")
          && argv[i + 1] != std::string("-i")
          && argv[i + 1] != std::string("-o")) {
        recorded_beamIntensity = std::stof(argv[i + 1]);
        beamintensity_set = true;
      } else {
        std::cerr << "You didn't give an argument for the beamintensity!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-a")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-b")
          && argv[i + 1] != std::string("-o")
          && argv[i + 1] != std::string("-i")) {
        around_halfaperture = std::stof(argv[i + 1]);
        around_halfaperture_set = true;
      } else {
        std::cerr << "You didn't give an argument for the half aperture around which the asymmetric scan was done!"
          << std::endl;
      }
    }
  }
  if (!inputfile_set || !outputfile_set || !beamintensity_set || !around_halfaperture_set) {
    std::cerr
      << "You didn't give the name for the inputfile, the outputfile, or the beam intensity or the half aperture you're interested in. Please try again!"
      << std::endl;
    exit(1);
  }
  std::cout << "Inputfile: " << inputfilename << std::endl;
  std::cout << "Output: " << outputfilename << std::endl;
  std::cout << "Half aperture scan point [mm]: " << around_halfaperture << std::endl;
  std::cout << "Beam intensity [10^8]: " << recorded_beamIntensity << std::endl;
  recorded_beamIntensity /= 100.0;//change back to a intensity unit of 10^10, as the intensity is given like this in the ROOT inputfiles


  TFile* inputfile = TFile::Open(inputfilename.c_str());
  TTree* Detector = nullptr;
  inputfile->GetObject("Tree_Detector1", Detector);

  //Set the branches
  Detector->SetBranchStatus("*", 0);
  Detector->SetBranchStatus("BeamIntensity",1);
  Detector->SetBranchStatus("CollAperture",1);
  Detector->SetBranchStatus("CollUpperJawPosition",1);
  Detector->SetBranchStatus("CollLowerJawPosition",1);
  Detector->SetBranchStatus("NoiseSubtractedSignal",1);
  Detector->SetBranchStatus("Voltage",1);

  float beamintensity = 0.0;
  float collaperture = 0.0;
  float lowerjawposition = 0.0;
  float upperjawposition = 0.0;
  float signal = 0;
  int voltage = 0;

  Detector->SetBranchAddress("BeamIntensity", &beamintensity);
  Detector->SetBranchAddress("CollAperture", &collaperture);
  Detector->SetBranchAddress("CollUpperJawPosition", &upperjawposition);
  Detector->SetBranchAddress("CollLowerJawPosition", &lowerjawposition);
  Detector->SetBranchAddress("NoiseSubtractedSignal", &signal);
  Detector->SetBranchAddress("Voltage", &voltage);
/*
  float UpperJawPosition[n], LowerJawPosition[n];
  if(around_halfaperture == 4){
    UpperJawPosition[n] = {2.6,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5};
    LowerJawPosition[n] = {5.4,5,4.75,4.5,4.25,4,3.75,3.5,3.25,3};
  }
  if(around_halfaperture == 5){
    UpperJawPosition[n] = {3,3.5,4,4.5,5,5.5,6,6.5,7};
    LowerJawPosition[n] = {7,6.5,6,5.5,5,4.5,4,3.5,3};
  }
  float JawPositionError[n] = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};
*/
  TCanvas* PDFcanvas = new TCanvas();
  std::string PDFTitle = outputfilename + "_AllPlots.pdf";
  PDFcanvas->Print(("output/"+PDFTitle+"[").c_str());
  
  TCanvas* canvas = new TCanvas();
  
  std::map< int, std::vector<int> > AverageSignalMap;
  std::string title = "Average signal as a function of the upper and lower jaw position;Upper jaw position [mm];Lower jaw position [mm];Average signal strength [a.u.]";
  TH2D* AsymmetricScanProfile = new TH2D("AsymmetricScan",title.c_str(),40,2.3,7.5,40,2.3,7.5);
  long long int const entries =  Detector->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    Detector->GetEntry(i);
    if(voltage > 0
        && collaperture >= 2*around_halfaperture-0.1 && collaperture <= 2*around_halfaperture+0.1
        && beamintensity >= recorded_beamIntensity-0.03 && beamintensity <= recorded_beamIntensity+0.03){
      //AsymmetricScanProfile->Fill(upperjawposition,lowerjawposition,signal);
      if(signal > 0){
        int global_bin = AsymmetricScanProfile->FindBin(upperjawposition,lowerjawposition);
        AverageSignalMap[global_bin].push_back(signal);
      }
    }
  }
  std::cout << "Map.size = " << AverageSignalMap.size() << std::endl;
  for(const auto& map_iterator : AverageSignalMap){
    float average = 0.0;
    for(int vector_iterator = 0; vector_iterator < map_iterator.second.size(); ++vector_iterator){
      average += map_iterator.second.at(vector_iterator);
    }  
    average /= (float)map_iterator.second.size();
    std::cout << "map_iterator.first = " << map_iterator.first << std::endl;
    std::cout << "map_iterator.second.size = " << map_iterator.second.size() << std::endl;
    std::cout << "average = " << average << std::endl;
    if(average>0) AsymmetricScanProfile->SetBinContent(map_iterator.first, average);
  }
  gStyle->SetOptStat(0);
  canvas->SetGrid();
  
  canvas->SetTheta(30);
  canvas->SetPhi(40);
  if(around_halfaperture == 4){
    AsymmetricScanProfile->GetXaxis()->SetRangeUser(2.3,5.7);
    AsymmetricScanProfile->GetYaxis()->SetRangeUser(2.3,5.7);
  }
  if(recorded_beamIntensity < 1) AsymmetricScanProfile->GetZaxis()->SetRangeUser(9000,15000);//beam intensity 0.84
  if(recorded_beamIntensity > 1) AsymmetricScanProfile->GetZaxis()->SetRangeUser(11500,17000);//beam intensity 1.01
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(6000,8400);//beam intensity 1.07, 4 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(9000,12500);//beam intensity 0.84, 4 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(12500,14700);//beam intensity 1.01, 4 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(6000,15000);//beam intensity 0.84, 5 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(9000,16900);//beam intensity 1.01, 5 mm
  AsymmetricScanProfile->Draw("lego 0");
  PRINT(canvas, outputfilename+"_lego", PDFTitle);
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(7100,7900);//beam intensity 1.07, 4 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(9700,11700);//beam intensity 0.84, 4 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(12800,13900);//beam intensity 1.01, 4 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(7000,15000);//beam intensity 0.84, 5 mm
  //AsymmetricScanProfile->GetZaxis()->SetRangeUser(11000,17000);//beam intensity 1.01, 5 mm
  AsymmetricScanProfile->Draw("colz");
  PRINT(canvas, outputfilename+"_colz", PDFTitle);
  
  PDFcanvas->Print(("output/"+PDFTitle+"]").c_str());
  inputfile->Close();
}  
