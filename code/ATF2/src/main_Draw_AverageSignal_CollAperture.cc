#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include "Style.h"

#include <vector>
#include <array>
#include <iostream>

int const n = 10; //Number of apertures that were recorded
void GetAverageSignals(float* Average, bool GetError, const float* beamintensity, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* collaperture_branch, const float* signal_branch, const int* voltage_branch);

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  //The input is a TTree ROOT file(s)
  //The output is .pdf and .C files

  std::string inputfilename;
  std::string outputfilename;

  bool inputfile_set = false;
  bool outputfile_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")){
      if (argv[i + 1] != std::string("-o")
          && argv[i + 1] != NULL) {
        inputfilename = argv[i + 1];
        inputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the inputfile(s)!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-o")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-i")) {
        outputfilename = argv[i + 1];
        outputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the sample!"
          << std::endl;
      }
    }

  }
  if (!inputfile_set || !outputfile_set) {
    std::cerr
      << "You didn't give the name for the inputfiles or the outfile. Please try again!"
      << std::endl;
    exit(1);
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
 
  float Aperture[n] = {3,4,5,6,7,8,9,10,11,12};
  float ApertureError[n] = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};
  
  //Fill the arrays with the average and the RMS of the signals from the TTree for the different beam intensities:
  float recorded_beamIntensity = 0.44;
  float SignalAverage[n];
  GetAverageSignals(SignalAverage, false, &recorded_beamIntensity, Aperture, n, Detector, &beamintensity, &collaperture, &signal, &voltage);
  float SignalAverageError[n];
  GetAverageSignals(SignalAverageError, true, &recorded_beamIntensity, Aperture, n, Detector, &beamintensity, &collaperture, &signal, &voltage);
  TGraphErrors* AverageSignal_CollAperture = new TGraphErrors(n,Aperture,SignalAverage,ApertureError,SignalAverageError);
  AverageSignal_CollAperture->SetTitle("Average signal strength for different beam halo collimator apertures;Collimator aperture [mm];Average RHUL cherenkov signal [a.u.]");
  AverageSignal_CollAperture->SetMarkerColor(4);
  AverageSignal_CollAperture->SetMarkerStyle(8);
  
  float recorded_beamIntensity2 = 0.21;
  float SignalAverage_secondIntensity[n];
  GetAverageSignals(SignalAverage_secondIntensity, false, &recorded_beamIntensity2, Aperture, n, Detector, &beamintensity, &collaperture, &signal, &voltage);
  float SignalAverageError_secondIntensity[n];
  GetAverageSignals(SignalAverageError_secondIntensity, true, &recorded_beamIntensity2, Aperture, n, Detector, &beamintensity, &collaperture, &signal, &voltage);
  TGraphErrors* AverageSignal_secondIntensity_CollAperture = new TGraphErrors(n,Aperture,SignalAverage_secondIntensity,ApertureError,SignalAverageError_secondIntensity);
  AverageSignal_secondIntensity_CollAperture->SetTitle("Average signal strength for different beam halo collimator apertures;Collimator aperture [mm];Average RHUL cherenkov signal [a.u.]");
  AverageSignal_secondIntensity_CollAperture->SetMarkerColor(4);
  AverageSignal_secondIntensity_CollAperture->SetMarkerStyle(8);

  //Plot the TGraphErrors for the different intensities onto the same canvas:
  TCanvas* canvas = new TCanvas();
  AverageSignal_CollAperture->Draw("APE");
  AverageSignal_secondIntensity_CollAperture->Draw("APE,SAME");

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
    Signal_CollAperture.emplace_back(new TH1D("signal-noise", title1D.c_str(),1000,0,10000));
  }

  long long int const entries =  tree->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    tree->GetEntry(i);
    if(*voltage_branch > 0){
      if(*beamintensity_branch > *beamintensity-0.1 && *beamintensity_branch < *beamintensity+0.1){
        
        for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
          //Fill the TH1 in the vector with signals for an aperture, that corresponds to the desired apertures in the aperture vector:
          if(*collaperture_branch > apertures[number_apertures]-0.1 && *collaperture_branch < apertures[number_apertures]+0.1){
            Signal_CollAperture.at(number_apertures)->Fill(*signal_branch);
          }
        }

      }
    }
  }

  for(size_t iterator; iterator<Signal_CollAperture.size();++iterator){
    //If the average signal is desired, get the mean from the signal distributions in the TH1 vector
    if (GetError==false){
      SignalAverage[iterator] = Signal_CollAperture.at(iterator)->GetMean(1);  
    }
    //If the RMS is desired, get the RMS from the signal distributions in the TH1 vector
    else{
      SignalAverage[iterator] = Signal_CollAperture.at(iterator)->GetRMS(1);  
    }
    delete Signal_CollAperture.at(iterator);
  }
}
