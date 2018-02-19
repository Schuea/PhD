#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "Style.h"

#include "UsefulFunctions.h"
#include "Style.h"

void Plot_everything(std::vector< std::string > inputfilenames, std::string histo_name, std::string outputfilename);

int main(int const argc, char const * const * const argv){
  UsePhDStyle();

	std::vector< std::string > inputfilenames;
	std::string outputfilename;

	if (argc < 3) {
		std::cerr << "Please provide the names of the input and output files: " << std::endl;
		std::cerr << "./CompareEnergy -i file1.root file2.root -o output" << std::endl;
		exit(2);
	}
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i")){
			if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-o")) {
				int j = 1;
				do {
					inputfilenames.push_back(argv[i + j]);
					++j;
				} while (argv[i + j] != NULL 
					  && argv[i + j] != std::string("-s")
            && argv[i + j] != std::string("-o") );
			} else {
				std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
				exit(2);
			}
		}
    if (argv[i] == std::string("-o")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-i")) {
				outputfilename = argv[i + 1];
			} else {
				std::cerr << "You didn't give an argument for the outputfilename!"
					<< std::endl;
				exit(2);
			}
		}
	}

  Plot_everything( inputfilenames, "Muon_Energy", outputfilename);

	return 0;
}

void Plot_everything(std::vector< std::string > inputfilenames, std::string histo_name, std::string outputfilename) {
  std::string x_title, y_title;
  x_title = "Energy [GeV]";
  y_title = "Number of muons";

  TCanvas* canvas = new TCanvas();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  gStyle->SetOptStat(0);

  TLegend* legend;
  legend = new TLegend(0.5,0.7,0.9,0.9);
  legend->SetFillColor(kWhite);
  legend->SetMargin(0.2);
  
  std::vector< TH1D* > histos;
	for(size_t no_files = 0; no_files < inputfilenames.size(); ++ no_files){
    TFile* inputfile = new TFile( inputfilenames.at(no_files).c_str() );
    std::ifstream ifile( inputfilenames.at(no_files).c_str() );
    if( (bool)ifile == 1){
      TH1D* temp = new TH1D();
      temp = (TH1D*) inputfile->Get( histo_name.c_str() );
      histos.push_back(temp);
    }
		else{
			std::cout<<"Error! File "<< inputfilenames.at(no_files) << " not found!";
			exit(1);
		}
  }
  double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
  std::string name;
  for(size_t no_histo = 0; no_histo < histos.size(); ++ no_histo){
    name = " ";
    histos.at(no_histo)->SetMaximum( max );
    histos.at(no_histo)->SetLineWidth(2);
    histos.at(no_histo)->SetMarkerSize(0.7);
    if(no_histo == 0){
      name = "ILC500, 5 spoilers";
      histos.at(no_histo)->SetLineColor(1);
      histos.at(no_histo)->SetMarkerColor(1);
      histos.at(no_histo)->SetMarkerStyle(4);
      histos.at(no_histo)->GetXaxis()->SetTitle( x_title.c_str() );
      histos.at(no_histo)->GetYaxis()->SetTitle( y_title.c_str() );
      //histos.at(no_histo)->Draw("hist,p");
      histos.at(no_histo)->Draw("e");
    }
    else{
      if(no_histo == 1){
        histos.at(no_histo)->SetLineColor(2);
        histos.at(no_histo)->SetMarkerColor(2);
        histos.at(no_histo)->SetMarkerStyle(8);
        name = "ILC500, 5 spoilers + wall";
      }
      if(no_histo == 2){
        histos.at(no_histo)->SetLineColor(3);
        histos.at(no_histo)->SetMarkerColor(3);
        histos.at(no_histo)->SetMarkerStyle(25);
        name = "ILC250, 5 spoilers";
      }
      if(no_histo == 3){
        histos.at(no_histo)->SetLineColor(4);
        histos.at(no_histo)->SetMarkerColor(4);
        histos.at(no_histo)->SetMarkerStyle(21);
        name = "ILC250, 5 spoilers + wall";
      }
      //histos.at(no_histo)->Draw("hist,p,SAMES");
      histos.at(no_histo)->Draw("e,SAMES");
    }
    legend->AddEntry(histos.at(no_histo),name.c_str(),"lep");
  }
  legend->Draw();
  canvas->Update();
  canvas->Print( ("output/Energy_Comparison_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/Energy_Comparison_"+outputfilename+".cxx").c_str() );
}
