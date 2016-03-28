#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"

#include "Style.h"

#include <vector>
#include <array>
#include <iostream>
#include <sstream>

int const n = 14; //Number of apertures that were recorded
float HistoMax = 0.0; //To find the maximum entry to the histogram, so that the y-axis can be scaled appropriately.
float HistoMin = 1000000.0; //To find the minimum entry to the histogram, so that the y-axis can be scaled appropriately.

void GetAverageSignals(float* Average, bool GetError, const float* beamintensity, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* collaperture_branch, const float* signal_branch, const int* voltage_branch);

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  gStyle->SetEndErrorSize(1.5);
  //The input is a TTree ROOT file(s)
  //The output is .pdf and .C files

  std::string inputfilename;
  std::string outputfilename;
  std::vector< float > recorded_beamIntensities;
  //std::vector< float > recorded_beamIntensities(6,0.0);//Initialize the vector with 6 zeroes, as there are up to 6 complete aperture scans (for different beam intensities)

  bool inputfile_set = false;
  bool outputfile_set = false;
  bool beamintensity_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")){
      if (argv[i + 1] != std::string("-o")
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
          && argv[i + 1] != std::string("-i")
          && argv[i + 1] != std::string("-o")) {
        int j = 1;
        while(argv[i + j] != NULL
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
  }
  if (!inputfile_set || !outputfile_set || !beamintensity_set) {
    std::cerr
      << "You didn't give the name for the inputfiles, the outputfile, or the beamintensity. Please try again!"
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
  Detector->SetBranchStatus("NoiseSubtractedSignal",1);
  Detector->SetBranchStatus("Voltage",1);

  float beamintensity = 0.0;
  float collaperture = 0.0;
  float signal = 0;
  int voltage = 0;

  Detector->SetBranchAddress("BeamIntensity", &beamintensity);
  Detector->SetBranchAddress("CollAperture", &collaperture);
  Detector->SetBranchAddress("NoiseSubtractedSignal", &signal);
  Detector->SetBranchAddress("Voltage", &voltage);
 
  float Aperture[n] = {6,7,8,9,10,11,12,13,14,16,18,20,22,24};
  float ApertureError[n] = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};
  
  //Fill the arrays with the average and the RMS/sqrt(N) of the signals from the TTree for the different beam intensities:
  std::vector< TGraphErrors*> All_TGraphErrors;
  for(int intensity_iterator = 0; intensity_iterator < recorded_beamIntensities.size(); ++intensity_iterator){

    float SignalAverage[n];
    GetAverageSignals(SignalAverage, false, &recorded_beamIntensities.at(intensity_iterator), Aperture, n, Detector, &beamintensity, &collaperture, &signal, &voltage);
    float SignalAverageError[n];
    GetAverageSignals(SignalAverageError, true, &recorded_beamIntensities.at(intensity_iterator), Aperture, n, Detector, &beamintensity, &collaperture, &signal, &voltage);

    TGraphErrors* AverageSignal_CollAperture = new TGraphErrors(n,Aperture,SignalAverage,ApertureError,SignalAverageError);
    AverageSignal_CollAperture->SetTitle("Average signal strength for different beam halo collimator apertures;Collimator aperture [mm];Average RHUL cherenkov signal [a.u.]");
    AverageSignal_CollAperture->SetMarkerColorAlpha(4 + intensity_iterator,0.5);//change the color for every new graph
    AverageSignal_CollAperture->SetMarkerSize(0.85);
    AverageSignal_CollAperture->SetMarkerStyle(8);
   
    //for(int point_iterator = 0; point_iterator < n; ++point_iterator){
    //  if(SignalAverage[point_iterator] < 0.1){
    //    AverageSignal_CollAperture->RemovePoint(point_iterator);//Don't draw all points that are not filled (because data was not taken for all apertures)
    //  }
    //}
    All_TGraphErrors.push_back(AverageSignal_CollAperture);
  }  
  
  //Plot the TGraphErrors for the different intensities onto the same canvas:
  TCanvas* canvas = new TCanvas();
  TLegend* legend = new TLegend(0.5,0.8,0.65,0.93);
  std::stringstream legend_text_unit;
  legend_text_unit << "*10^10";
  
  for(int intensity_iterator = 0; intensity_iterator < recorded_beamIntensities.size(); ++intensity_iterator){
    All_TGraphErrors.at(intensity_iterator)->SetMaximum(HistoMax + 0.1*HistoMax);
    All_TGraphErrors.at(intensity_iterator)->SetMinimum(1000);
    //All_TGraphErrors.at(intensity_iterator)->SetMinimum(HistoMin - 0.1*HistoMin);//will always be 0 because there are points with value 0!
  
    std::stringstream legend_text;
    legend_text.str(std::string());
    legend_text << recorded_beamIntensities.at(intensity_iterator) << legend_text_unit.str();
    legend->AddEntry(All_TGraphErrors.at(intensity_iterator),legend_text.str().c_str(),"lep");

    if(intensity_iterator == 0) All_TGraphErrors.at(intensity_iterator)->Draw("AP");
    else All_TGraphErrors.at(intensity_iterator)->Draw("P SAME");
  }
  legend->Draw();

  std::string outputname1 = "output/" + outputfilename + ".pdf";
  std::string outputname2 = "output/" + outputfilename + ".cxx";
  canvas->Print(outputname1.c_str());
  canvas->Print(outputname2.c_str());

  inputfile->Close();
}  

void GetAverageSignals(float* SignalAverage, bool GetError, const float* beamintensity, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* collaperture_branch, const float* signal_branch, const int* voltage_branch){
  std::string title1D = "RHUL cherenkov detector signal;Signal [a.u.];Count";
  std::vector<TH1D*> Signal_CollAperture;
  //Push back as many TH1 histograms as needed in order to have one for each aperture
  for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
    Signal_CollAperture.emplace_back(new TH1D("signal-noise", title1D.c_str(),1000,0,100000));
  }

  long long int const entries =  tree->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    tree->GetEntry(i);
    if(*voltage_branch > 0){
      if(*beamintensity_branch >= *beamintensity-0.03 && *beamintensity_branch <= *beamintensity+0.03){
        
        for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
          //Fill the TH1 in the vector with signals for an aperture, that corresponds to the desired apertures in the aperture vector:
          if(*collaperture_branch > apertures[number_apertures]-0.1 && *collaperture_branch < apertures[number_apertures]+0.1){
            Signal_CollAperture.at(number_apertures)->Fill(*signal_branch);
          }
        }

      }
      else continue;
    }
    else continue;
  }

  for(size_t iterator; iterator < Signal_CollAperture.size();++iterator){
    //If the average signal is desired, get the mean from the signal distributions in the TH1 vector
    if (GetError==false){
      SignalAverage[iterator] = Signal_CollAperture.at(iterator)->GetMean();  
      if (SignalAverage[iterator] > HistoMax) HistoMax = SignalAverage[iterator];
      if (SignalAverage[iterator] < HistoMin) HistoMin = SignalAverage[iterator];
    }
    //If the error is desired, get the RMS from the signal distributions in the TH1 vector devided by the square root of entries->standard deviation of the mean
    else{
      SignalAverage[iterator] = Signal_CollAperture.at(iterator)->GetRMS();  
      //SignalAverage[iterator] = Signal_CollAperture.at(iterator)->GetRMS()/std::sqrt(Signal_CollAperture.at(iterator)->GetEntries());  
    }
    delete Signal_CollAperture.at(iterator);
  }
}
