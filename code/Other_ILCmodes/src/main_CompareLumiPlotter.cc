#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <sstream>

int main(int const argc, char const * const * const argv){

	std::vector< std::string > inputfilenames;
	std::string outputfilename;

	if (argc < 5) {
		std::cerr << "Please provide the names of the input and output files: " << std::endl;
		std::cerr << "./CompareLumiPlotter -i file1.root file2.root -o output" << std::endl;
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
  new_outputfilename << "output/" << outputfilename << ".root";
	TFile* ROOTFile = new TFile(new_outputfilename.str().c_str(),"RECREATE","Luminosity histogram");

  TCanvas* canvas = new TCanvas();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  gStyle->SetOptStat(0);

  TLegend* legend = new TLegend(0.13,0.7,0.55,0.88);
  std::string name;

	for(size_t no_files = 0; no_files < inputfilenames.size(); ++ no_files){
    name = " ";

    TFile* inputfile = new TFile( inputfilenames.at(no_files).c_str() );
    std::ifstream ifile( inputfilenames.at(no_files).c_str() );
    if( (bool)ifile == 1){
      TH1F* temp = (TH1F*) inputfile->Get("Lumi");
      temp->SetLineWidth(2);
      if(no_files == 0){
        name = "set  2: TDR";
        temp->SetLineColor(1);
        temp->GetXaxis()->SetTitle("E_CM [GeV]");
        temp->GetYaxis()->SetTitle("Luminosity");
        temp->Draw();
      }
      else{
        if(no_files == 1){
          temp->SetLineColor(2);
          name = "set  4: TDR + Emittance_x";
        }
        if(no_files == 2){
          temp->SetLineColor(3);
          name = "set 15: TDR + Emittance_x + Beta_x";
        }
        if(no_files == 3){
          temp->SetLineColor(4);
          name = "set 16: TDR + Emittance_x + Beta_x + Beta_y";
        }
        temp->Draw("SAMES");
      }
      legend->AddEntry(temp,name.c_str(),"l");
      //inputfile->Close();
		}
		else{
			std::cout<<"Error! File "<< inputfilenames.at(no_files) << " not found!";
			exit(1);
		}

	}
  legend->Draw();

  canvas->Print( ("output/"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/"+outputfilename+".cxx").c_str() );

	std::cout << "Output will be created: " << outputfilename << std::endl;
	ROOTFile->Write();
	ROOTFile->Close();
	return 0;
}
