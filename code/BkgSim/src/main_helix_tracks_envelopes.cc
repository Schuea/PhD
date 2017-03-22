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
double * InvertYaxis(int const n, double *y){
  double * result = new double[n];
  for(int i = 0; i < n; ++i){
    result[i] = -y[i];
  }
  return result;
}
TGraph * MakeTheGraph(TH2D *h, float const acceptance){
  int const sizeX = h->GetNbinsX();
  int const sizeY = h->GetNbinsY();
  double const binSizeY = 29*2/200.; //Maybe be wrong
  double const binSizeX = 300.0/800.; //Maybe be wrong
  int const n = sizeX-1;
  double x1[n], y1[n];
  for(int i = 1; i < sizeX; ++i){ //Maybe be wrong
    int total(0);
    for(int j = 1; j < sizeY; ++j){ //Maybe be wrong
      total += h->GetBinContent(i,j);
    }
    x1[i-1] = 0;
    y1[i-1] = 0;
    int toTenPercent(0);
    for(int j = 1; j < sizeY; ++j){ //Maybe be wrong
      toTenPercent += h->GetBinContent(i,j) + h->GetBinContent(i, sizeY-j);
      if(toTenPercent > (1.0-acceptance)*total){
        x1[i-1] = i * binSizeX;
        y1[i-1] = (sizeY/2.0 - j) * binSizeY;
        break;
      }
    }
  }
  TGraph * graph = new TGraph(n,x1,y1);
  graph->SetMinimum(-15);
  graph->SetMaximum(15);
  graph->SetMarkerStyle(7);
  return graph;
}

TLegend* MakeLegend(vector< TGraph* > const & all_graphs){
  TLegend* leg = new TLegend(0.78,0.55,0.9,0.9);
  leg->AddEntry(all_graphs.at(0),"68%","p");
  leg->AddEntry(all_graphs.at(2),"90%","p");
  leg->AddEntry(all_graphs.at(4),"95%","p");
  leg->AddEntry(all_graphs.at(6),"99%","p");
  leg->AddEntry(all_graphs.at(8),"99.7%","p");
  leg->AddEntry(all_graphs.at(10),"99.9%","p");
  leg->AddEntry(all_graphs.at(12),"99.99%","p");
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
  for(size_t i = 0; i < all_graphs.size(); ++i){
    all_graphs.at(i)->SetMarkerStyle(7);
    all_graphs.at(i)->SetMarkerColor(i/2 + 1);
    all_graphs.at(i)->SetTitle(("Envelopes outlining fractions of helix tracks from pair backgrounds;z [mm];"+direction+" [mm]").c_str());
  }
}

vector< TGraph* > GetAllGraphs(TH2D * hx){
  int const sizeX = hx->GetNbinsX();
  int const n = sizeX-1;
  TGraph *graph_68 = MakeTheGraph(hx, 0.68);
  TGraph *graph_90 = MakeTheGraph(hx, 0.90);
  TGraph *graph_95 = MakeTheGraph(hx, 0.95);
  TGraph *graph_99 = MakeTheGraph(hx, 0.99);
  TGraph *graph_997 = MakeTheGraph(hx, 0.997);
  TGraph *graph_999 = MakeTheGraph(hx, 0.999);
  TGraph *graph_9999 = MakeTheGraph(hx, 0.9999);

  TGraph *graph_68_neg = new TGraph(n, graph_68->GetX(), InvertYaxis(n, graph_68->GetY()));
  TGraph *graph_90_neg = new TGraph(n, graph_90->GetX(), InvertYaxis(n, graph_90->GetY()));
  TGraph *graph_95_neg = new TGraph(n, graph_95->GetX(), InvertYaxis(n, graph_95->GetY()));
  TGraph *graph_99_neg = new TGraph(n, graph_99->GetX(), InvertYaxis(n, graph_99->GetY()));
  TGraph *graph_997_neg = new TGraph(n, graph_997->GetX(), InvertYaxis(n, graph_997->GetY()));
  TGraph *graph_999_neg = new TGraph(n, graph_999->GetX(), InvertYaxis(n, graph_999->GetY()));
  TGraph *graph_9999_neg = new TGraph(n, graph_9999->GetX(), InvertYaxis(n, graph_9999->GetY()));

  vector< TGraph* > all_graphs;
  all_graphs.push_back(graph_68);
  all_graphs.push_back(graph_68_neg);
  all_graphs.push_back(graph_90);
  all_graphs.push_back(graph_90_neg);
  all_graphs.push_back(graph_95);
  all_graphs.push_back(graph_95_neg);
  all_graphs.push_back(graph_99);
  all_graphs.push_back(graph_99_neg);
  all_graphs.push_back(graph_997);
  all_graphs.push_back(graph_997_neg);
  all_graphs.push_back(graph_999);
  all_graphs.push_back(graph_999_neg);
  all_graphs.push_back(graph_9999);
  all_graphs.push_back(graph_9999_neg);

  return all_graphs;
}


int main(int const argc, char const * const * const argv){
  UsePhDStyle();
  std::string specialname("");

  try{
    if(argc <=2) throw;
    for(int i = 0; i < argc; ++i){
      if(string(argv[i]) == "-i") specialname = argv[i+1];
    }
  } catch(...){
    cerr << "Input arguments were invalid" << endl;
    cerr << "Give the specifications of the input file, e.g. 1bunch_500GeV_3T" << endl;
    exit(1);
  }

  std::string plane1 = "Helix_tracks_xz";
  std::string plane2 = "Helix_tracks_yz";
  TFile *file = TFile::Open( ("output/Helix_in_beampipe_"+specialname+".root").c_str() );
  TH2D *hx = (TH2D*)file->Get( plane1.c_str() );
  TH2D *hy = (TH2D*)file->Get( plane2.c_str() );

  vector< TGraph* > all_graphs1 = GetAllGraphs(hx);
  vector< TGraph* > all_graphs2 = GetAllGraphs(hy);

  SetAllGraphStyles(all_graphs1,"x");
  SetAllGraphStyles(all_graphs2,"y");

  TCanvas c;
  c.SetGridx();
  c.SetGridy();

  DrawAllGraphs(all_graphs1);

  vector< TLine* > all_lines = GetAllLines();
  SetAllLineStyles(all_lines);
  DrawAllLines(all_lines);

  TLegend* leg1 = MakeLegend(all_graphs1);
  //TLegend *leg1 = new TLegend(0.18,0.11,0.38,0.3);
  //leg1->SetBorderSize(0);
  //leg1->SetFillColorAlpha(0,0);
  leg1->Draw();
  //resultx->Draw("hist");
  c.Print( ("output/HelixEnvelopes_"+plane1+"_"+specialname+".pdf").c_str() );
  c.Print( ("output/HelixEnvelopes_"+plane1+"_"+specialname+".cxx").c_str() );

  DrawAllGraphs(all_graphs2);
  DrawAllLines(all_lines);
  TLegend* leg2 = MakeLegend(all_graphs2);
  leg2->Draw();
  c.Print( ("output/HelixEnvelopes_"+plane2+"_"+specialname+".pdf").c_str() );
  c.Print( ("output/HelixEnvelopes_"+plane2+"_"+specialname+".cxx").c_str() );
  return 0;
}
