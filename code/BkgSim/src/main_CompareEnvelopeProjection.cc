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

int main(int const argc, char const * const * const argv){
  UsePhDStyle();

	std::vector< std::string > inputfilenames;
	std::string outputfilename;

	if (argc < 5) {
		std::cerr << "Please provide the names of the input and output files: " << std::endl;
		std::cerr << "./CompareEnvelopeProjections -i file1.root file2.root -o output" << std::endl;
		exit(2);
	}
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i")){
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-o")) {
				int j = 1;
				do {
					inputfilenames.push_back(argv[i + j]);
					++j;
				} while (argv[i + j] != NULL && argv[i + j] != std::string("-o") );
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


  std::stringstream new_outputfilename;
  new_outputfilename << "output/HelixEnvelope_Projection_Comparison_" << outputfilename << ".root";
	TFile* ROOTFile = new TFile(new_outputfilename.str().c_str(),"RECREATE","Comparison Pair Envelope Prrojections");

  TCanvas* canvas = new TCanvas();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  gStyle->SetOptStat(0);

  TLegend* legend = new TLegend(0.58,0.65,0.98,0.85);
  legend->SetFillColor(kWhite);
  std::string name;
  
	for(size_t no_files = 0; no_files < inputfilenames.size(); ++ no_files){
    name = " ";

    TFile* inputfile = new TFile( inputfilenames.at(no_files).c_str() );
    std::ifstream ifile( inputfilenames.at(no_files).c_str() );
    if( (bool)ifile == 1){
      TH1F* temp = (TH1F*) inputfile->Get("Projection_xz_1Kink");
      temp->SetLineWidth(2);
      temp->SetMarkerSize(0.7);
      if(no_files == 0){
        name = "set  2: TDR";
        temp->SetLineColor(1);
        temp->SetMarkerColor(1);
        temp->SetMarkerStyle(4);
        temp->GetXaxis()->SetTitle("x [mm]");
        temp->GetYaxis()->SetTitle("Number of pair background particles");
        temp->Draw();
      }
      else{
        if(no_files == 1){
          temp->SetLineColor(2);
          temp->SetMarkerColor(2);
          temp->SetMarkerStyle(8);
          name = "set  4: TDR + Emittance_x";
        }
        if(no_files == 2){
          temp->SetLineColor(3);
          temp->SetMarkerColor(3);
          temp->SetMarkerStyle(25);
          name = "set 15: TDR + Emittance_x + Beta_x";
        }
        if(no_files == 3){
          temp->SetLineColor(4);
          temp->SetMarkerColor(4);
          temp->SetMarkerStyle(21);
          name = "set 16: TDR + Emittance_x + Beta_x + Beta_y";
        }
        temp->Draw("SAMES");
      }
      legend->AddEntry(temp,name.c_str(),"p");
		}
		else{
			std::cout<<"Error! File "<< inputfilenames.at(no_files) << " not found!";
			exit(1);
		}

	}
  canvas->Update();

  float ymin = -0.001; 
  float ymax = 4.0*std::pow(10.0,4.0);

  TLine* line = new TLine(-12,ymin,-12,ymax);
  TLine* line2 = new TLine(12,ymin,12,ymax);
  line->SetLineColor(kMagenta);
  line->SetLineWidth(2);
  line2->SetLineColor(kMagenta);
  line2->SetLineWidth(2);
  line->Draw();
  line2->Draw();

  legend->Draw();
  
  canvas->Print( ("output/HelixEnvelope_Projection_Comparison_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/HelixEnvelope_Projection_Comparison_"+outputfilename+".cxx").c_str() );

	std::cout << "Output will be created: " << outputfilename << std::endl;
	ROOTFile->Write();
	ROOTFile->Close();
	return 0;
}
