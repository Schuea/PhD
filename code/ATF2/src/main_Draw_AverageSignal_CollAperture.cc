#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include "Style.h"

#include <iostream>

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


  std::string title1D = "RHUL cherenkov detector signal;Signal [a.u.];Count";
  TH1D* Signal_CollAperture3mm = new TH1D("signal_3mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture4mm = new TH1D("signal_4mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture5mm = new TH1D("signal_5mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture6mm = new TH1D("signal_6mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture7mm = new TH1D("signal_7mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture8mm = new TH1D("signal_8mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture9mm = new TH1D("signal_9mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture10mm = new TH1D("signal_10mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture11mm = new TH1D("signal_11mm",title1D.c_str(),1000,0,10000);
  TH1D* Signal_CollAperture12mm = new TH1D("signal_12mm",title1D.c_str(),1000,0,10000);


  TFile* inputfile = TFile::Open(inputfilename.c_str());
  TTree* Detector = nullptr;
  inputfile->GetObject("Tree_Detector1", Detector);

  //Set the branches
  Detector->SetBranchStatus("*", 0);
  Detector->SetBranchStatus("BeamIntensity",1);
  Detector->SetBranchStatus("CollAperture",1);
  Detector->SetBranchStatus("Signal",1);
  Detector->SetBranchStatus("Voltage",1);

  float beamintensity = 0.0;
  float collaperture = 0.0;
  int signal = 0;
  int voltage = 0;

  Detector->SetBranchAddress("BeamIntensity", &beamintensity);
  Detector->SetBranchAddress("CollAperture", &collaperture);
  Detector->SetBranchAddress("Signal", &signal);
  Detector->SetBranchAddress("Voltage", &voltage);
  
  long long int const entries =  Detector->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    Detector->GetEntry(i);
    if(voltage > 0 &&
        beamintensity > 0 && 
        beamintensity < 0.48 && beamintensity > 0.40){
       if(collaperture > 2.9 && collaperture < 3.1) Signal_CollAperture3mm->Fill(signal);
       if(collaperture > 3.9 && collaperture < 4.1) Signal_CollAperture4mm->Fill(signal);
       if(collaperture > 4.9 && collaperture < 5.1) Signal_CollAperture5mm->Fill(signal);
       if(collaperture > 5.9 && collaperture < 6.1) Signal_CollAperture6mm->Fill(signal);
       if(collaperture > 6.9 && collaperture < 7.1) Signal_CollAperture7mm->Fill(signal);
       if(collaperture > 7.9 && collaperture < 8.1) Signal_CollAperture8mm->Fill(signal);
       if(collaperture > 8.9 && collaperture < 9.1) Signal_CollAperture9mm->Fill(signal);
       if(collaperture > 9.9 && collaperture < 10.1) Signal_CollAperture10mm->Fill(signal);
       if(collaperture > 10.9 && collaperture < 11.1) Signal_CollAperture11mm->Fill(signal);
       if(collaperture > 11.9 && collaperture < 12.1) Signal_CollAperture12mm->Fill(signal);
    }
  }

  float SignalAverage[10];
  float SignalAverageError[10];
  SignalAverage[0] = Signal_CollAperture3mm->GetMean(1);  
  SignalAverageError[0] = Signal_CollAperture3mm->GetRMS(1);  
  SignalAverage[1] = Signal_CollAperture4mm->GetMean(1);  
  SignalAverageError[1] = Signal_CollAperture4mm->GetRMS(1);  
  SignalAverage[2] = Signal_CollAperture5mm->GetMean(1);  
  SignalAverageError[2] = Signal_CollAperture5mm->GetRMS(1);  
  SignalAverage[3] = Signal_CollAperture6mm->GetMean(1);  
  SignalAverageError[3] = Signal_CollAperture6mm->GetRMS(1);  
  SignalAverage[4] = Signal_CollAperture7mm->GetMean(1);  
  SignalAverageError[4] = Signal_CollAperture7mm->GetRMS(1);  
  SignalAverage[5] = Signal_CollAperture8mm->GetMean(1);  
  SignalAverageError[5] = Signal_CollAperture8mm->GetRMS(1);  
  SignalAverage[6] = Signal_CollAperture9mm->GetMean(1);  
  SignalAverageError[6] = Signal_CollAperture9mm->GetRMS(1);  
  SignalAverage[7] = Signal_CollAperture10mm->GetMean(1);  
  SignalAverageError[7] = Signal_CollAperture10mm->GetRMS(1);  
  SignalAverage[8] = Signal_CollAperture11mm->GetMean(1);  
  SignalAverageError[8] = Signal_CollAperture11mm->GetRMS(1);  
  SignalAverage[9] = Signal_CollAperture12mm->GetMean(1);  
  SignalAverageError[9] = Signal_CollAperture11mm->GetRMS(1);  
  float Aperture[10] = {3,4,5,6,7,8,9,10,11,12};
  float ApertureError[10] = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  TGraphErrors* AverageSignal_CollAperture = new TGraphErrors(10,Aperture,SignalAverage,ApertureError,SignalAverageError);
  AverageSignal_CollAperture->SetTitle("Average signal strength for different beam halo collimator apertures;Collimator aperture [mm];Average RHUL cherenkov signal [a.u.]");
  AverageSignal_CollAperture->SetMarkerColor(4);
  AverageSignal_CollAperture->SetMarkerStyle(8);

  TCanvas* canvas = new TCanvas();
  AverageSignal_CollAperture->Draw("APE");

  std::string outputname1 = "output/" + outputfilename + ".pdf";
  std::string outputname2 = "output/" + outputfilename + ".cxx";
  canvas->Print(outputname1.c_str());
  canvas->Print(outputname2.c_str());

  inputfile->Close();
  delete Signal_CollAperture3mm;
  delete Signal_CollAperture4mm;
  delete Signal_CollAperture5mm;
  delete Signal_CollAperture6mm;
  delete Signal_CollAperture7mm;
  delete Signal_CollAperture8mm;
  delete Signal_CollAperture9mm;
  delete Signal_CollAperture10mm;
  delete Signal_CollAperture11mm;
  delete Signal_CollAperture12mm;
}  
