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
  bool signalfile_set = false;
  bool noisefile_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-n")){
        SIGNALinputfilename = argv[i + 1];
        signalfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the SIGNALinputfilename!"
          << std::endl;
        std::exit(1);
      }
    }
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-s")){
        NOISEinputfilename = argv[i + 1];
        noisefile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the NOISEinputfilename!"
          << std::endl;
        std::exit(1);
      }
    }
  }
  if (argc < 4 || (!noisefile_set || !signalfile_set)){
    std::cerr << "You didn't give the required arguments: -s signalfile.root -n noisefile.root !"
      << std::endl;
    std::exit(1);
  }


  std::string title1D = "RHUL cherenkov detector noise;Noise [a.u.];Count";
  TH1D* Noise_voltage0V = new TH1D("noise_0V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage500V = new TH1D("noise_500V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage600V = new TH1D("noise_600V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage650V = new TH1D("noise_650V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage700V = new TH1D("noise_700V",title1D.c_str(),1000,0,10000);
  TH1D* Noise_voltage750V = new TH1D("noise_750V",title1D.c_str(),1000,0,10000);
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
    if(voltage_noise > 0 && voltage_noise < 3) Noise_voltage0V->Fill(signal_noise);
    if(voltage_noise > 497 && voltage_noise < 503) Noise_voltage500V->Fill(signal_noise);
    if(voltage_noise > 597 && voltage_noise < 603) Noise_voltage600V->Fill(signal_noise);
    if(voltage_noise > 647 && voltage_noise < 653) Noise_voltage650V->Fill(signal_noise);
    if(voltage_noise > 697 && voltage_noise < 703) Noise_voltage700V->Fill(signal_noise);
    if(voltage_noise > 747 && voltage_noise < 753) Noise_voltage750V->Fill(signal_noise);
    if(voltage_noise > 797 && voltage_noise < 803) Noise_voltage800V->Fill(signal_noise);
    if(voltage_noise > 897 && voltage_noise < 903) Noise_voltage900V->Fill(signal_noise);
    if(voltage_noise > 997 && voltage_noise < 1003) Noise_voltage1000V->Fill(signal_noise);
    if(voltage_noise > 1097 && voltage_noise < 1103) Noise_voltage1100V->Fill(signal_noise);
    if(voltage_noise > 1197 && voltage_noise < 1203) Noise_voltage1200V->Fill(signal_noise);
    if(voltage_noise > 1297 && voltage_noise < 1303) Noise_voltage1300V->Fill(signal_noise);
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
  
  long long int entries = Detector->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    Detector->GetEntry(i);
    if(voltage > 0 && voltage < 3) NoiseSubtracted_signal1 = signal - Noise_voltage0V->GetMean(1);
    if(voltage > 497 && voltage < 503) NoiseSubtracted_signal1 = signal - Noise_voltage500V->GetMean(1);
    if(voltage > 597 && voltage < 603) NoiseSubtracted_signal1 = signal - Noise_voltage600V->GetMean(1);
    if(voltage > 647 && voltage < 653) NoiseSubtracted_signal1 = signal - Noise_voltage650V->GetMean(1);
    if(voltage > 697 && voltage < 703) NoiseSubtracted_signal1 = signal - Noise_voltage700V->GetMean(1);
    if(voltage > 747 && voltage < 753) NoiseSubtracted_signal1 = signal - Noise_voltage750V->GetMean(1);
    if(voltage > 797 && voltage < 803) NoiseSubtracted_signal1 = signal - Noise_voltage800V->GetMean(1);
    if(voltage > 897 && voltage < 903) NoiseSubtracted_signal1 = signal - Noise_voltage900V->GetMean(1);
    if(voltage > 997 && voltage < 1003) NoiseSubtracted_signal1 = signal - Noise_voltage1000V->GetMean(1);
    if(voltage > 1097 && voltage < 1103) NoiseSubtracted_signal1 = signal - Noise_voltage1100V->GetMean(1);
    if(voltage > 1197 && voltage < 1203) NoiseSubtracted_signal1 = signal - Noise_voltage1200V->GetMean(1);
    if(voltage > 1297 && voltage < 1303) NoiseSubtracted_signal1 = signal - Noise_voltage1300V->GetMean(1);
    NoiseSubtracted_Signal->Fill();
  }
  Detector->Write("",TObject::kOverwrite);
  SIGNALinputfile->Write();
  SIGNALinputfile->Close();

  return 0;
}
