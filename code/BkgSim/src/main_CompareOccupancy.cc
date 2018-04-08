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
#include <vector>

#include "UsefulFunctions.h"
#include "Subdetector_class_new.h"
#include "Style.h"

void Plot_everything(std::vector< std::string > inputfilenames, std::string layer, std::string outputfilename);
void Plot_Comparison_histos(std::vector< std::string > inputfilenames, std::string histoname, std::string x_title, std::string y_title);

int main(int const argc, char const * const * const argv){
  UsePhDStyle();

	std::vector< std::string > inputfilenames;
	std::string outputfilename;

	if (argc < 5) {
		std::cerr << "Please provide the names of the input and output files: " << std::endl;
		std::cerr << "./CompareOccupancies -i file1.root file2.root -o output" << std::endl;
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
            && argv[i + j] != std::string("-o") );
			} else {
				std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
				exit(2);
			}
		}
		if (argv[i] == std::string("-o")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-s")
					&& argv[i + 1] != std::string("-i")) {
				outputfilename = argv[i + 1];
			} else {
				std::cerr << "You didn't give an argument for the outputfilename!"
					<< std::endl;
				exit(2);
			}
		}

	}

  Plot_everything( inputfilenames, "All_layers", outputfilename);
  //for(int layer = 0; layer < det.getNumberOfLayers(); ++layer){
  for(int layer = 0; layer < 1; ++layer){
    std::ostringstream num;
    num << "Layer_" << layer;
    Plot_everything( inputfilenames, num.str(), outputfilename);
  }
	
  return 0;
}

void Plot_everything(std::vector< std::string > inputfilenames, std::string layer, std::string outputfilename) {
  TCanvas* canvas = new TCanvas();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  gStyle->SetOptStat(0);

 
  Plot_Comparison_histos( inputfilenames, layer, "Number of hits per cell", "Number of cells" );
  canvas->Update();
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_"+outputfilename+".cxx").c_str() );

  Plot_Comparison_histos( inputfilenames, layer+"_losthits", "Assumed buffer depth", "Number of hits lost" );
  canvas->Update();
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_losthits_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_losthits_"+outputfilename+".cxx").c_str() );

  Plot_Comparison_histos( inputfilenames, layer+"_deadcells", "Assumed buffer depth", "Number of dead cells" );
  canvas->Update();
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_deadcells_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_deadcells_"+outputfilename+".cxx").c_str() );
  
  if (layer.compare("All_layers") == 0 ){
  Plot_Comparison_histos( inputfilenames, layer+"_wrt_#cells", "Number of hits per cell", "Number of cells" );
  canvas->Update();
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_wrt_#cells_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_wrt_#cells_"+outputfilename+".cxx").c_str() );
  }
  else {
  Plot_Comparison_histos( inputfilenames, layer+"_numcells", "Number of hits per cell", "Number of cells" );
  canvas->Update();
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_numcells_"+outputfilename+".pdf").c_str() );
  canvas->Print( ("output/Occupancy_Comparison_"+layer+"_numcells_"+outputfilename+".cxx").c_str() );
  }
}

void Plot_Comparison_histos(std::vector< std::string > inputfilenames, std::string histoname, std::string x_title, std::string y_title){
  TLegend* legend;
  if(inputfilenames.size()<3) legend = new TLegend(0.62,0.75,0.8,0.86);
  else                         legend = new TLegend(0.3,0.6267,0.8,0.86);
  legend->SetBorderSize(1);
  legend->SetFillColor(kWhite);
  legend->SetMargin(0.1);
  //legend->SetMargin(0.2);
  
  std::vector< TH1F* > histos;
	for(size_t no_files = 0; no_files < inputfilenames.size(); ++ no_files){
    TFile* inputfile = new TFile( inputfilenames.at(no_files).c_str() );
    std::ifstream ifile( inputfilenames.at(no_files).c_str() );
    if( (bool)ifile == 1){
      TH1F* temp = new TH1F();
      temp = (TH1F*) inputfile->Get( histoname.c_str() );
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
    std::size_t found = histoname.find("cells");
    if(found != std::string::npos){//substring was found
      histos.at(no_histo)->SetMinimum( std::pow(10.0,-9.0) );
      //histos.at(no_histo)->SetMinimum( std::pow(10.0,-12.0) );
    }else{
      histos.at(no_histo)->SetMinimum( std::pow(10.0,-6.0) );
    }
    histos.at(no_histo)->SetMaximum( max );
    histos.at(no_histo)->SetLineWidth(2);
    histos.at(no_histo)->SetMarkerSize(0.7);
    if(no_histo == 0){
      //name = "ILC500";
      name = "set (TDR)";
      //name = "set (A): old L*, w/o antiDiD";
      histos.at(no_histo)->SetLineColor(1);
      histos.at(no_histo)->SetMarkerColor(1);
      histos.at(no_histo)->SetMarkerStyle(4);
      //histos.at(no_histo)->SetLineColor(kCyan+3);
      //histos.at(no_histo)->SetMarkerColor(kCyan+3);
      histos.at(no_histo)->GetXaxis()->SetTitle( x_title.c_str() );
      histos.at(no_histo)->GetYaxis()->SetTitle( y_title.c_str() );
      histos.at(no_histo)->Draw("e");
    }
    else{
      if(no_histo == 1){
        //histos.at(no_histo)->SetLineColor(kPink-5);
        //histos.at(no_histo)->SetMarkerColor(kPink-5);
        histos.at(no_histo)->SetLineColor(2);
        histos.at(no_histo)->SetMarkerColor(2);
        histos.at(no_histo)->SetMarkerStyle(8);
        //histos.at(no_histo)->SetMarkerStyle(26);
        //histos.at(no_histo)->SetMarkerSize(0.8);
        //name = "ILC250";
        //name = "set (A): old L*, w antiDiD";
        name = "set (A): TDR + Emittance_x";
        //name = "ILC500 LumiUp";
      }
      if(no_histo == 2){
        histos.at(no_histo)->SetLineColor(3);
        histos.at(no_histo)->SetMarkerColor(3);
        histos.at(no_histo)->SetMarkerStyle(25);
        //name = "set (A): new L*, w/o antiDiD";
        name = "set (B): TDR + Emittance_x + Beta_x";
      }
      if(no_histo == 3){
        histos.at(no_histo)->SetLineColor(4);
        histos.at(no_histo)->SetMarkerColor(4);
        histos.at(no_histo)->SetMarkerStyle(21);
        //name = "set (A): new L*, w antiDiD";
        name = "set (C): TDR + Emittance_x + Beta_x + Beta_y";
      }
      //histos.at(no_histo)->Draw("hist,p,SAMES");
      histos.at(no_histo)->Draw("e,SAMES");
    }
    legend->AddEntry(histos.at(no_histo),name.c_str(),"p");
    //legend->AddEntry(histos.at(no_histo),name.c_str(),"lep");
  }
  legend->Draw();
}
