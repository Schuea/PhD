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

int const n = 16; //Number of apertures that were recorded
float HistoMax = 0.0; //To find the maximum entry to the histogram, so that the y-axis can be scaled appropriately.
float HistoMin = 10000000.0; //To find the minimum entry to the histogram, so that the y-axis can be scaled appropriately.

void GetAverageSignals(float* SignalAverage, bool GetError, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* firstjawposition_branch, const float* secondjawposition_branch, const float* signal_branch, const int* voltage_branch);
void GetAverageSignals(bool normalizing_ON, float* SignalAverage, bool GetError, const float* beamintensity, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* firstjawposition_branch, const float* secondjawposition_branch, const float* signal_branch, const int* voltage_branch);

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  //The input is a TTree ROOT file(s)
  //The output is .pdf and .C files

  std::string inputfilename;
  std::string outputfilename;
  int normalizing_ON = 0;
        std::cout <<__LINE__<< " normalizing_ON = " << normalizing_ON << std::endl; 
  std::vector< float > recorded_beamIntensities;

  bool inputfile_set = false;
  bool outputfile_set = false;
  bool beamintensity_set = false;
  bool normalizing_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i")){
      if (argv[i + 1] != std::string("-o")
          && argv[i + 1] != std::string("-b")
          && argv[i + 1] != std::string("-n")
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
          && argv[i + 1] != std::string("-n")
          && argv[i + 1] != std::string("-i")) {
        outputfilename = argv[i + 1];
        outputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the output filename!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-b")
          && argv[i + 1] != std::string("-o")
          && argv[i + 1] != std::string("-i")) {
        normalizing_ON = std::stof(argv[i + 1]);
        std::cout <<__LINE__<< " normalizing_ON = " << normalizing_ON << std::endl; 
        normalizing_set = true;
      } else {
        std::cerr << "You didn't say if the data shall be normalized by charge!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-b")) {
      if (argv[i + 1] != NULL
          && argv[i + 1] != std::string("-n")
          && argv[i + 1] != std::string("-i")
          && argv[i + 1] != std::string("-o")) {
        int j = 1;
        while(argv[i + j] != NULL
            && argv[i + j] != std::string("-n")
            && argv[i + j] != std::string("-i")
            && argv[i + j] != std::string("-o")) {
          std::string str(argv[i + j]);
          recorded_beamIntensities.push_back(std::atof(str.c_str()));
          j++;
        }
        beamintensity_set = true;
      } else {
        std::cerr << "You didn't give an argument for the beamintensity!"
          << std::endl;
      }
    }

  }
  if (!inputfile_set || !outputfile_set || !normalizing_set) {
    std::cerr
      << "You didn't give the name for the inputfile, the outputfile, or if the data shall be normalized. Please try again!"
      << std::endl;
    exit(1);
  }
  if (normalizing_ON == 0 && !beamintensity_set){
    std::cerr << "Please specify which beam intensities you are interested in: e.g. -b 0.8" << std::endl;
    exit(1);
  }
  std::cout << "Inputfile: " << inputfilename << std::endl;
  std::cout << "Output: " << outputfilename << std::endl;
        std::cout <<__LINE__<< " normalizing_ON = " << normalizing_ON << std::endl; 
  if (normalizing_ON == 0){
    std::cout << "Beam intensities [10^10]: " << std::endl;
    for(int intensity_iterator = 0; intensity_iterator < recorded_beamIntensities.size(); ++intensity_iterator){
      std::cout << recorded_beamIntensities.at(intensity_iterator) << std::endl;
      recorded_beamIntensities.at(intensity_iterator);
      //recorded_beamIntensities.at(intensity_iterator) /= 100.0;//change back to a intensity unit of 10^10, as the intensity is given like this in the ROOT inputfiles
    }
  }
  else{
    std::cout << "----------------" << std::endl;
    std::cout << "The data will be normalized by charge!" << std::endl;
    std::cout << "----------------" << std::endl;
  }

  TFile* inputfile = TFile::Open(inputfilename.c_str());
  TTree* Detector = nullptr;
  inputfile->GetObject("Tree_Detector1", Detector);

  //Set the branches
  Detector->SetBranchStatus("*", 0);
  Detector->SetBranchStatus("BeamIntensity",1);
  Detector->SetBranchStatus("CollUpperJawPosition",1);
  Detector->SetBranchStatus("CollLowerJawPosition",1);
  Detector->SetBranchStatus("NoiseSubtractedSignal",1);
  Detector->SetBranchStatus("Voltage",1);

  float beamintensity = 0.0;
  float lowerjawposition = 0.0;
  float upperjawposition = 0.0;
  float signal = 0;
  int voltage = 0;

  Detector->SetBranchAddress("BeamIntensity", &beamintensity);
  Detector->SetBranchAddress("CollUpperJawPosition", &upperjawposition);
  Detector->SetBranchAddress("CollLowerJawPosition", &lowerjawposition);
  Detector->SetBranchAddress("NoiseSubtractedSignal", &signal);
  Detector->SetBranchAddress("Voltage", &voltage);

  float JawPosition[n] = {2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12};
  float JawPositionError[n] = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04};

  //Fill the arrays with the average and the RMS/sqrt(N) of the signals from the TTree for the different beam intensities:
  std::vector< TGraphErrors*> All_TGraphErrors;
 
  int intensity_iterator = 0; 
  do{
    //UPPER JAW:
    float SignalAverage[n];
    float SignalAverageError[n];
    if(recorded_beamIntensities.size()>0){
      GetAverageSignals(normalizing_ON, SignalAverage, false, &recorded_beamIntensities.at(intensity_iterator), JawPosition, n, Detector, &beamintensity, &upperjawposition, &lowerjawposition, &signal, &voltage);
      GetAverageSignals(normalizing_ON, SignalAverageError, true, &recorded_beamIntensities.at(intensity_iterator), JawPosition, n, Detector, &beamintensity, &upperjawposition, &lowerjawposition, &signal, &voltage);
    }
    else{
      GetAverageSignals(SignalAverage, false, JawPosition, n, Detector, &beamintensity, &upperjawposition, &lowerjawposition, &signal, &voltage);
      GetAverageSignals(SignalAverageError, true, JawPosition, n, Detector, &beamintensity, &upperjawposition, &lowerjawposition, &signal, &voltage);
    }
  std::cout << __LINE__ << std::endl;

    TGraphErrors* AverageSignal_CollAperture = new TGraphErrors(n,JawPosition,SignalAverage,JawPositionError,SignalAverageError);
    AverageSignal_CollAperture->SetTitle("Average signal strength for different beam halo collimator apertures;Jaw position [mm];Average RHUL cherenkov signal [a.u.]");
    AverageSignal_CollAperture->SetMarkerColorAlpha(4 + intensity_iterator*2,0.5);//change the color for every new graph
    AverageSignal_CollAperture->SetMarkerSize(0.85);
    AverageSignal_CollAperture->SetMarkerStyle(8);

    All_TGraphErrors.push_back(AverageSignal_CollAperture);

    //LOWER_JAW:
    float SignalAverage2[n];
    float SignalAverageError2[n];
    if(recorded_beamIntensities.size()>0){
    GetAverageSignals(normalizing_ON, SignalAverage2, false, &recorded_beamIntensities.at(intensity_iterator), JawPosition, n, Detector, &beamintensity, &lowerjawposition, &upperjawposition, &signal, &voltage);
    GetAverageSignals(normalizing_ON, SignalAverageError2, true, &recorded_beamIntensities.at(intensity_iterator), JawPosition, n, Detector, &beamintensity, &lowerjawposition, &upperjawposition, &signal, &voltage);
    }
    else{
    GetAverageSignals(SignalAverage2, false, JawPosition, n, Detector, &beamintensity, &lowerjawposition, &upperjawposition, &signal, &voltage);
    GetAverageSignals(SignalAverageError2, true, JawPosition, n, Detector, &beamintensity, &lowerjawposition, &upperjawposition, &signal, &voltage);
    }
  std::cout << __LINE__ << std::endl;

    TGraphErrors* AverageSignal_CollAperture2 = new TGraphErrors(n,JawPosition,SignalAverage2,JawPositionError,SignalAverageError2);
    AverageSignal_CollAperture2->SetTitle("Average signal strength for different beam halo collimator apertures;Collimator aperture [mm];Average RHUL cherenkov signal [a.u.]");
    AverageSignal_CollAperture2->SetMarkerColorAlpha(4 + intensity_iterator*2 + 1,0.5);//change the color for every new graph
    AverageSignal_CollAperture2->SetMarkerSize(0.85);
    AverageSignal_CollAperture2->SetMarkerStyle(8);

    All_TGraphErrors.push_back(AverageSignal_CollAperture2);
    ++intensity_iterator;
  }
  while(intensity_iterator < recorded_beamIntensities.size());

  std::cout << __LINE__ << std::endl;

  //Plot the TGraphErrors for the different intensities onto the same canvas:
  TCanvas* canvas = new TCanvas();
  TLegend* legend = new TLegend(0.43,0.13,0.8,0.36);
  std::stringstream legend_text_unit;
  legend_text_unit << "*10^10";

  for(int graph_iterator = 0; graph_iterator < All_TGraphErrors.size(); ++graph_iterator){
    All_TGraphErrors.at(graph_iterator)->SetMaximum(HistoMax + 0.1*HistoMax);
    All_TGraphErrors.at(graph_iterator)->SetMinimum(HistoMin - 0.3*HistoMin);
  std::cout << __LINE__ << std::endl;

    std::stringstream legend_text1, legend_text2;
    int intensity_iterator = std::floor(((float)graph_iterator)/2.0);//There are two graphs per intensity -> there are twice as many entries in the GraphVector as there are in the IntensityVector -> only count the intensity_iterator up for every number%2==0 basically
    if (normalizing_ON == false && recorded_beamIntensities.size()>0){
      legend_text1 << recorded_beamIntensities.at(intensity_iterator) << legend_text_unit.str() << ", upper jaw moving";
      legend_text2 << recorded_beamIntensities.at(intensity_iterator) << legend_text_unit.str() << ", lower jaw moving";
    }
    else{
      legend_text1 << "Run "<< intensity_iterator+1 <<", upper jaw moving";
      legend_text2 << "Run "<< intensity_iterator+1 <<", lower jaw moving";
    }
  std::cout << __LINE__ << std::endl;

    //All even graphs are for the upper jaw with legend_text1, all odd graphs for the lower jaw with legend_text2:
    if(graph_iterator == 0){//Specify the very first graph separately because here mustn't be written 'SAME' in drawing options
      All_TGraphErrors.at(graph_iterator)->Draw("AP");
      legend->AddEntry(All_TGraphErrors.at(graph_iterator),legend_text1.str().c_str(),"lep");
    }
    else if(graph_iterator > 0 && graph_iterator%2 == 0){
      All_TGraphErrors.at(graph_iterator)->Draw("P SAME");
      legend->AddEntry(All_TGraphErrors.at(graph_iterator),legend_text1.str().c_str(),"lep");
    }
    else{
      All_TGraphErrors.at(graph_iterator)->Draw("P SAME");
      legend->AddEntry(All_TGraphErrors.at(graph_iterator),legend_text2.str().c_str(),"lep");
    }
  }
  legend->Draw();
  std::cout << __LINE__ << std::endl;

  std::string outputname1 = "output/" + outputfilename + ".pdf";
  std::string outputname2 = "output/" + outputfilename + ".cxx";
  canvas->Print(outputname1.c_str());
  canvas->Print(outputname2.c_str());

  inputfile->Close();
}  

void GetAverageSignals(bool normalizing_ON, float* SignalAverage, bool GetError, const float* beamintensity, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* firstjawposition_branch, const float* secondjawposition_branch, const float* signal_branch, const int* voltage_branch){

  std::string title1D = "RHUL cherenkov detector signal;Signal [a.u.];Count";
  std::vector<TH1D*> Signal_CollAperture;
  //Push back as many TH1 histograms as needed in order to have one for each aperture
  for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
    Signal_CollAperture.emplace_back(new TH1D("signal-noise", title1D.c_str(),500,0,1000000));
  }

  long long int const entries =  tree->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    tree->GetEntry(i);
    if(*voltage_branch > 0
        && *secondjawposition_branch >= 11.5 && *secondjawposition_branch <= 12.5) //one jaw was held on the open position, while the other one was moved
    {
      if(normalizing_ON){
        for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
          //Fill the TH1 in the vector with signals for an aperture, that corresponds to the desired apertures in the aperture vector:
          if(*firstjawposition_branch > apertures[number_apertures]-0.1 && *firstjawposition_branch < apertures[number_apertures]+0.1){
            //Normalizing by charge
            Signal_CollAperture.at(number_apertures)->Fill(*signal_branch/ *beamintensity_branch);
          }//end if
        }//end for
      }//end if(normalizing_ON)  
      else{    
        if(*beamintensity_branch >= *beamintensity-0.15 && *beamintensity_branch <= *beamintensity+0.15){
          for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
            //Fill the TH1 in the vector with signals for an aperture, that corresponds to the desired apertures in the aperture vector:
            if(*firstjawposition_branch > apertures[number_apertures]-0.1 && *firstjawposition_branch < apertures[number_apertures]+0.1){
              Signal_CollAperture.at(number_apertures)->Fill(*signal_branch);
            }//end if
          }//end for
        }//end if
      }//end else
    }//end if (*voltage_branch > 0 &&..
    else continue;
  }//end for

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
  }//end for
}//end function

void GetAverageSignals(float* SignalAverage, bool GetError, const float apertures[], const int num_apertures, TTree* tree, const float* beamintensity_branch, const float* firstjawposition_branch, const float* secondjawposition_branch, const float* signal_branch, const int* voltage_branch){

  std::string title1D = "RHUL cherenkov detector signal;Signal [a.u.];Count";
  std::vector<TH1D*> Signal_CollAperture;
  //Push back as many TH1 histograms as needed in order to have one for each aperture
  for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
    Signal_CollAperture.emplace_back(new TH1D("signal-noise", title1D.c_str(),500,0,1000000));
  }

  long long int const entries =  tree->GetEntries();
  for (long long int i = 0; i < entries; ++i){
    tree->GetEntry(i);
    if(*voltage_branch > 0
        && *secondjawposition_branch >= 11.5 && *secondjawposition_branch <= 12.5) //one jaw was held on the open position, while the other one was moved
    {
      for(int number_apertures = 0; number_apertures < num_apertures; ++number_apertures){
        //Fill the TH1 in the vector with signals for an aperture, that corresponds to the desired apertures in the aperture vector:
        if(*firstjawposition_branch > apertures[number_apertures]-0.1 && *firstjawposition_branch < apertures[number_apertures]+0.1){
          //Normalizing by charge
          Signal_CollAperture.at(number_apertures)->Fill(*signal_branch/ *beamintensity_branch);
        }//end if
      }//end for
    }//end if (*voltage_branch > 0 &&..
    else continue;
  }//end for

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
  }//end for
}
