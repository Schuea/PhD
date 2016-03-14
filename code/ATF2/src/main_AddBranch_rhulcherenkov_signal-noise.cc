#include "TFile.h"
#include "TChain.h"
#include "TObject.h"
#include "TH1.h"


#include <iostream>
#include <fstream>
#include <sstream>

int main(int const argc, char const * const * const argv){

  std::string SIGNALinputfilename;
  std::string NOISEinputfilename;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-n")
          && argv[i + 1] != std::string("-o")){
        SIGNALinputfilename = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the SIGNALinputfilename!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-s")
          && argv[i + 1] != std::string("-o")){
        NOISEinputfilename = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the NOISEinputfilename!"
          << std::endl;
      }
    }
  }

  std::string title1D = "RHUL cherenkov detector noise;Noise [a.u.];Count";
  TH1D* Noise_voltage0V = new TH1D("noise_0V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage500V = new TH1D("noise_500V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage700V = new TH1D("noise_700V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage800V = new TH1D("noise_800V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage900V = new TH1D("noise_900V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage1000V = new TH1D("noise_1000V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage1100V = new TH1D("noise_1100V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage1200V = new TH1D("noise_1200V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage1300V = new TH1D("noise_1300V",title1D.c_str(),1000,0,10000);

  TFile* NOISEinputfile = TFile::Open(NOISEinputfilename.c_str());
  TTree* Detector_Noise = nullptr;
  NOISEinputfile->GetObject("Tree_Detector1",Detector_Noise);

  //Set the branches
  Detector_Noise->SetBranchStatus("*", 0);
  Detector_Noise->SetBranchStatus("Signal",1);
  Detector_Noise->SetBranchStatus("Voltage",1);

  int signal_noise = 0;
  int voltage_noise = 0;

  Detector_Noise->SetBranchAddress("Signal", &signal_noise);
  Detector_Noise->SetBranchAddress("Voltage", &voltage_noise);
  
  long long int const entries_Noise =  Detector_Noise->GetEntries();
  for (long long int i = 0; i < entries_Noise; ++i){
    Detector_Noise->GetEntry(i);
    if(voltage_noise < 1) Noise_voltage0V->Fill(signal_noise);
    if(voltage_noise > 409 && voltage_noise < 501) Noise_voltage500V->Fill(signal_noise);
    if(voltage_noise > 609 && voltage_noise < 701) Noise_voltage700V->Fill(signal_noise);
    if(voltage_noise > 709 && voltage_noise < 801) Noise_voltage800V->Fill(signal_noise);
    if(voltage_noise > 809 && voltage_noise < 901) Noise_voltage900V->Fill(signal_noise);
    if(voltage_noise > 909 && voltage_noise < 1001) Noise_voltage1000V->Fill(signal_noise);
    if(voltage_noise > 1009 && voltage_noise < 1101) Noise_voltage1100V->Fill(signal_noise);
    if(voltage_noise > 1109 && voltage_noise < 1201) Noise_voltage1200V->Fill(signal_noise);
    if(voltage_noise > 1209 && voltage_noise < 1301) Noise_voltage1300V->Fill(signal_noise);
  }
  NOISEinputfile->Close();

  TFile* SIGNALinputfile = TFile::Open(SIGNALinputfilename.c_str(),"UPDATE");
  TTree* Detector = nullptr;
  SIGNALinputfile->GetObject("Tree_Detector1",Detector);
  
  Detector->SetBranchStatus("Signal",1);
  Detector->SetBranchStatus("Voltage",1);
  int signal = 0;
  int voltage = 0;
  Detector->SetBranchAddress("Signal", &signal);
  Detector->SetBranchAddress("Voltage", &voltage);
  
  float NoiseSubtracted_signal1 = 0.0;
  TBranch* NoiseSubtracted_Signal = Detector->Branch("NoiseSubtractedSignal",&NoiseSubtracted_signal1, "NoiseSubtractedSignal/F");
  
  std::cout << "Noise_voltage0V->GetMean(1): " << Noise_voltage0V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage500V->GetMean(1): " << Noise_voltage500V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage700V->GetMean(1): " << Noise_voltage700V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage800V->GetMean(1): " << Noise_voltage800V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage900V->GetMean(1): " << Noise_voltage900V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage1000V->GetMean(1): " << Noise_voltage1000V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage1100V->GetMean(1): " << Noise_voltage1100V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage1200V->GetMean(1): " << Noise_voltage1200V->GetMean(1) << std::endl;
  std::cout << "Noise_voltage1300V->GetMean(1): " << Noise_voltage1300V->GetMean(1) << std::endl;

  long long int entries = Detector->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    Detector->GetEntry(i);
    std::cout << "signal=" << signal << std::endl;
    if(voltage < 1) NoiseSubtracted_signal1 = signal - Noise_voltage0V->GetMean(1);
    if(voltage > 409 && voltage < 501) NoiseSubtracted_signal1 = signal - Noise_voltage500V->GetMean(1);
    if(voltage > 609 && voltage < 701) NoiseSubtracted_signal1 = signal - Noise_voltage700V->GetMean(1);
    if(voltage > 709 && voltage < 801) NoiseSubtracted_signal1 = signal - Noise_voltage800V->GetMean(1);
    if(voltage > 809 && voltage < 901) NoiseSubtracted_signal1 = signal - Noise_voltage900V->GetMean(1);
    if(voltage > 909 && voltage < 1001) NoiseSubtracted_signal1 = signal - Noise_voltage1000V->GetMean(1);
    if(voltage > 1009 && voltage < 1101) NoiseSubtracted_signal1 = signal - Noise_voltage1100V->GetMean(1);
    if(voltage > 1109 && voltage < 1201) NoiseSubtracted_signal1 = signal - Noise_voltage1200V->GetMean(1);
    if(voltage > 1209 && voltage < 1301) NoiseSubtracted_signal1 = signal - Noise_voltage1300V->GetMean(1);
    std::cout << "NoiseSubtracted_signal1 = " << NoiseSubtracted_signal1 << std::endl;
    NoiseSubtracted_Signal->Fill();
  }
  Detector->Write("",TObject::kOverwrite);
  SIGNALinputfile->Write();
  SIGNALinputfile->Close();

  return 0;
}
