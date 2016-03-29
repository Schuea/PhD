#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLegend.h"

#include "Style.h"
#include "UsefulFunctions.h"

#include <vector>
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
  std::vector< float > recorded_beamIntensities;
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
        int j = 1;
        while(argv[i + j] != NULL
            && argv[i + j] != std::string("-a")
            && argv[i + j] != std::string("-i")
            && argv[i + j] != std::string("-o")) {
          recorded_beamIntensities.push_back(std::stoi(argv[i + j]));
          j++;
        }
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
        around_halfaperture = std::stoi(argv[i + 1]);
        around_halfaperture_set = true;
      } else {
        std::cerr << "You didn't give an argument for the half aperture around which the asymmetric scan was done!"
          << std::endl;
      }
    }
  }
  if (!inputfile_set || !outputfile_set || !beamintensity_set) {
    std::cerr
      << "You didn't give the name for the inputfile, the outputfile, or the beam intensity or the half aperture you're interested in. Please try again!"
      << std::endl;
    exit(1);
  }
  std::cout << "Inputfile: " << inputfilename << std::endl;
  std::cout << "Output: " << outputfilename << std::endl;
  std::cout << "Beam intensities [10^8]: " << std::endl;
  for(int intensity_iterator = 0; intensity_iterator < recorded_beamIntensities.size(); ++intensity_iterator){
    std::cout << recorded_beamIntensities.at(intensity_iterator) << std::endl;
    recorded_beamIntensities.at(intensity_iterator) /= 100.0;//change back to a intensity unit of 10^10, as the intensity is given like this in the ROOT inputfiles
  }


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
  TCanvas* canvas = new TCanvas();
  std::string title = "Signal as a function of the upper and lower jaw position;Upper jaw position [mm];Lower jaw position [mm];Weighted signal strength [a.u.]";
  TH2D* AsymmetricScanProfile = new TH2D("AsymmetricScan",title.c_str(),40,2.3,5.7,40,2.3,5.7);
  AsymmetricScanProfile->GetZaxis()->SetRangeUser(6200000,7150000);
  long long int const entries =  Detector->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    Detector->GetEntry(i);
    if(voltage > 0
        && collaperture >= 2*around_halfaperture-0.1 && collaperture <= 2*around_halfaperture+0.1
        && beamintensity >= recorded_beamIntensities.at(0)-0.04 && beamintensity <= recorded_beamIntensities.at(0)+0.04){
      AsymmetricScanProfile->Fill(upperjawposition,lowerjawposition,signal);
    }
  }
  AsymmetricScanProfile->Draw("colz");
  gStyle->SetOptStat(0);
  canvas->SetGrid();
 
  std::string PDFTitle = outputfilename + "_AllPlots.pdf";
  PDFcanvas->Print(("output/"+PDFTitle+"[").c_str());
  PRINT(canvas, outputfilename, PDFTitle);
  PDFcanvas->Print(("output/"+PDFTitle+"]").c_str());
  inputfile->Close();
}  
