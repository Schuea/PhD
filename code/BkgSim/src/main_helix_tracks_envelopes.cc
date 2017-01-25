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
  double const binSizeY = 29*2/100.; //Maybe be wrong
  double const binSizeX = 300.0/10000.; //Maybe be wrong
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

int main(){
  UsePhDStyle();
  time_t time1 = time(NULL);
  TFile *file = TFile::Open("output/Helix_in_beampipe_1bunch_wo_momentum_cuts.root");
  TH2D *hx = (TH2D*)file->Get("Helix_tracks_xz");
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

  graph_68_neg->SetMarkerStyle(7);
  graph_90_neg->SetMarkerStyle(7);
  graph_95_neg->SetMarkerStyle(7);
  graph_99_neg->SetMarkerStyle(7);
  graph_997_neg->SetMarkerStyle(7);
  graph_999_neg->SetMarkerStyle(7);
  graph_9999_neg->SetMarkerStyle(7);

  graph_68->SetTitle("Envelopes outlining fractions of helix tracks from pair backgrounds;z [mm];x [mm]");

  graph_68->SetMarkerColor(1);
  graph_68_neg->SetMarkerColor(1);
  graph_90->SetMarkerColor(2);
  graph_90_neg->SetMarkerColor(2);
  graph_95->SetMarkerColor(3);
  graph_95_neg->SetMarkerColor(3);
  graph_99->SetMarkerColor(4);
  graph_99_neg->SetMarkerColor(4);
  graph_997->SetMarkerColor(6);
  graph_997_neg->SetMarkerColor(6);
  graph_999->SetMarkerColor(7);
  graph_999_neg->SetMarkerColor(7);
  graph_9999->SetMarkerColor(8);
  graph_9999_neg->SetMarkerColor(8);

	TLine *line = new TLine(0,12,62.5,12);
	TLine *nline = new TLine(0,-12,62.5,-12);
	TLine *line2 = new TLine(62.5,12,115,15);
	TLine *nline2 = new TLine(62.5,-12,115,-15);
	line->SetLineColor(2);
	nline->SetLineColor(2);
	line2->SetLineColor(2);
	nline2->SetLineColor(2);

  TCanvas c;
  c.SetGridx();
  c.SetGridy();

  graph_68->Draw("AP");
  graph_68_neg->Draw("PSAME");
  graph_90->Draw("PSAME");
  graph_90_neg->Draw("PSAME");
  graph_95->Draw("PSAME");
  graph_95_neg->Draw("PSAME");
  graph_99->Draw("PSAME");
  graph_99_neg->Draw("PSAME");
  graph_997->Draw("PSAME");
  graph_997_neg->Draw("PSAME");
  graph_999->Draw("PSAME");
  graph_999_neg->Draw("PSAME");
  graph_9999->Draw("PSAME");
  graph_9999_neg->Draw("PSAME");

	line->Draw();
	nline->Draw();
	line2->Draw();
	nline2->Draw();

  //TLegend *leg = new TLegend(0.18,0.11,0.38,0.3);
	TLegend* leg = new TLegend(0.75,0.7,0.95,0.9);
  //leg->SetBorderSize(0);
  //leg->SetFillColorAlpha(0,0);
  leg->AddEntry(graph_68,"68%","p");
  leg->AddEntry(graph_90,"90%","p");
  leg->AddEntry(graph_95,"95%","p");
  leg->AddEntry(graph_99,"99%","p");
  leg->AddEntry(graph_997,"99.7%","p");
  leg->AddEntry(graph_999,"99.9%","p");
  leg->AddEntry(graph_9999,"99.99%","p");
  leg->Draw();
  //resultx->Draw("hist");
  c.Print("output/HelixEnvelopes_xz.pdf");
  c.Print("output/HelixEnvelopes_xz.cxx");

  return 0;
}
