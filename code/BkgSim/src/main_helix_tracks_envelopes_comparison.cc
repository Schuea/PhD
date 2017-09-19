#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"

#include <ctime>
#include <iostream>
#include "Style.h"

using namespace std;

TLegend* MakeLegend(vector< TGraph* > const & all_graphs){
  TLegend* leg = new TLegend(0.65,0.55,0.9,0.9);
  leg->AddEntry(all_graphs.at(0),"99.9% set 2","p");
  leg->AddEntry(all_graphs.at(4),"99.9% set 4","p");
  leg->AddEntry(all_graphs.at(8),"99.9% set 15","p");
  leg->AddEntry(all_graphs.at(12),"99.9% set 16","p");
  leg->AddEntry(all_graphs.at(2),"99.99% set 2","p");
  leg->AddEntry(all_graphs.at(6),"99.99% set 4","p");
  leg->AddEntry(all_graphs.at(10),"99.99% set 15","p");
  leg->AddEntry(all_graphs.at(14),"99.99% set 16","p");
  //leg->AddEntry(all_graphs.at(0),"99.9% 500GeV","p");
  //leg->AddEntry(all_graphs.at(4),"99.9% 350GeV","p");
  //leg->AddEntry(all_graphs.at(8),"99.9% 250GeV","p");
  //leg->AddEntry(all_graphs.at(2),"99.99% 500GeV","p");
  //leg->AddEntry(all_graphs.at(6),"99.99% 350GeV","p");
  //leg->AddEntry(all_graphs.at(10),"99.99% 250GeV","p");
	return leg;
}

void SetAllLineStyles(vector< TLine* > & all_lines){
  for(int i = 0; i < all_lines.size(); ++i){
    all_lines.at(i)->SetLineColor(2);
    all_lines.at(i)->SetLineWidth(6);
  }
}

vector< TLine* > GetAllLines(){
  vector< TLine* > all_lines;
	all_lines.push_back(new TLine(0,12,62.5,12));
	all_lines.push_back(new TLine(0,-12,62.5,-12));
	all_lines.push_back(new TLine(62.5,12,115,15));
	all_lines.push_back(new TLine(62.5,-12,115,-15));
	return all_lines;
}

void DrawAllLines(vector< TLine* > const & all_lines){
  for(size_t i = 0; i < all_lines.size(); ++i){
    all_lines.at(i)->Draw();
  }
}

void DrawAllGraphs(vector< TGraph* > const & all_graphs){
  all_graphs.at(0)->Draw("AP");
  for(size_t i = 1; i < all_graphs.size(); ++i){
    all_graphs.at(i)->Draw("PSAME");
  }
}

void SetAllGraphStyles(vector< TGraph* > &all_graphs, std::string direction){
    int n = 1;
    int m = 1;
  for(size_t i = 0; i < all_graphs.size(); ++i){
    //Comment out the following 2 lines if you want to also see the negative y axis:
    all_graphs.at(i)->GetHistogram()->SetMaximum(15);
    all_graphs.at(i)->GetHistogram()->SetMinimum(0);
    all_graphs.at(i)->SetMarkerStyle(7);
    if(i/2 % 2 == 0) {
      if(n < 4) all_graphs.at(i)->SetMarkerColor(kMagenta + 4/n);
      else all_graphs.at(i)->SetMarkerColor(kMagenta - 9);
      if(i % 2 == 0) ++n;
    }
    else{
      if(m < 4) all_graphs.at(i)->SetMarkerColor(kCyan + 4/m);
      else all_graphs.at(i)->SetMarkerColor(kCyan - 9);
      if(i % 2 == 0) ++m;
    }
    all_graphs.at(i)->SetTitle(("Envelopes outlining fractions of helix tracks from pair backgrounds;z [mm];"+direction+" [mm]").c_str());
  }
}

int main(int const argc, char const * const * const argv){
  UsePhDStyle();
  std::vector< std::string > specialnames;
  std::string outputname;

  if(argc <=4){
    cerr << "Input arguments were invalid" << endl;
    cerr << "Give the specifications of the input file, e.g. 1bunch_500GeV_3T, and the outputfile specs" << endl;
    exit(1);
  }
  for(int i = 1; i < argc; ++i){
    if(string(argv[i]) == "-i" ){
      int j = 1;
      while( j <= 4 ){
        if(string(argv[i+j]) != "-i" 
            && argv[i+1] != NULL
            && string(argv[i+j]) != "-o"){
          specialnames.push_back( argv[i+j] );
          j++;
        }else break;
      };
    }
    if(string(argv[i]) == "-o" 
        && string(argv[i+1]) != "-o"
        && string(argv[i+1]) != "-i"){
      outputname =  argv[i+1] ;
    }
  }

  std::cout << outputname << std::endl;
  for( size_t file_iterator = 0; file_iterator < specialnames.size(); ++ file_iterator){
    std::cout << specialnames.at(file_iterator) << std::endl;
  }

  std::vector< TGraph* > gr_x;
  std::vector< TGraph* > gr_y;

  for( size_t file_iterator = 0; file_iterator < specialnames.size(); ++ file_iterator){
    TFile *file1 = TFile::Open( ("output/HelixEnvelopes_xz_"+specialnames.at(file_iterator)+".root").c_str() );
    TFile *file2 = TFile::Open( ("output/HelixEnvelopes_yz_"+specialnames.at(file_iterator)+".root").c_str() );

    gr_x.push_back( (TGraph*)file1->Get( "gr999" ));
    gr_x.push_back( (TGraph*)file1->Get( "gr999neg" ));
    gr_x.push_back( (TGraph*)file1->Get( "gr9999" ));
    gr_x.push_back( (TGraph*)file1->Get( "gr9999neg" ));

    gr_y.push_back( (TGraph*)file2->Get( "gr999" ));
    gr_y.push_back( (TGraph*)file2->Get( "gr999neg" ));
    gr_y.push_back( (TGraph*)file2->Get( "gr9999" ));
    gr_y.push_back( (TGraph*)file2->Get( "gr9999neg" ));

    file1->Close();
    file2->Close();
  }

  SetAllGraphStyles(gr_x,"x");
  SetAllGraphStyles(gr_y,"y");

	gStyle->SetOptStat(0);
  TCanvas c;
  c.SetGridx();
  c.SetGridy();

  DrawAllGraphs(gr_x);

  vector< TLine* > all_lines = GetAllLines();
  SetAllLineStyles(all_lines);
  DrawAllLines(all_lines);

  TLegend* leg1 = MakeLegend(gr_x);
  leg1->Draw();
  c.Print( ("output/HelixEnvelopes_COMPARISON_xz_"+outputname+".pdf").c_str() );
  c.Print( ("output/HelixEnvelopes_COMPARISON_xz_"+outputname+".cxx").c_str() );

  DrawAllGraphs(gr_y);
  DrawAllLines(all_lines);
  TLegend* leg2 = MakeLegend(gr_y);
  leg2->Draw();
  c.Print( ("output/HelixEnvelopes_COMPARISON_yz_"+outputname+".pdf").c_str() );
  c.Print( ("output/HelixEnvelopes_COMPARISON_yz_"+outputname+".cxx").c_str() );
  return 0;
}
