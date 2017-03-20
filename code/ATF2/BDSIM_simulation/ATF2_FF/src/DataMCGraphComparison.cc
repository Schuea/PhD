#include "Style.h"

#include "TLegend.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

int main(int const argc, char const * const * const argv){
  UsePhDStyle();
  TFile *file1 = TFile::Open(argv[1]);
  TFile *file2 = TFile::Open(argv[2]);
  TGraphErrors *graph1 = (TGraphErrors*)file1->Get("Graph");
  TGraphErrors *graph2 = (TGraphErrors*)file2->Get("Graph");

  //Normalise graph 2 to be like graph 1 for the average of the last three points
  float average1 =
    graph1->GetY()[graph1->GetN()-1] +
    graph1->GetY()[graph1->GetN()-2] +
    graph1->GetY()[graph1->GetN()-3];
  average1 /= 3;

  int N = graph2->GetN();
  float average2 =
    graph2->GetY()[N-1] +
    graph2->GetY()[N-2] +
    graph2->GetY()[N-3];
  average2 /= 3;

  float scale = average1/average2;
  cout << "average1 is " << average1 << endl;
  cout << "average2 is " << average2 << endl;
  cout << "scale is " << scale << endl;
  float x[N], ex[N], y[N], ey[N];
  for(int i = 0; i < N; ++i){
    x[i] = graph2->GetX()[i];
    y[i] = graph2->GetY()[i] * scale;
    ex[i] = graph2->GetEX()[i];
    ey[i] = graph2->GetEY()[i] * scale;
    cout << "x[" << i << "] is " << x[i] << endl;
    cout << "y[" << i << "] is " << y[i] << endl;
  }

  TGraphErrors *scaledgraph = new TGraphErrors(N,x,y,ex,ey);

  TCanvas *c = new TCanvas();
  gStyle->SetOptStat(0);
  graph1->SetMarkerColor(kRed);
  scaledgraph->SetMarkerColor(kBlue);
  graph1->Draw("AP");
  scaledgraph->Draw("PSAME");
  TLegend leg(0.6,0.6,0.77,0.77);
  leg.AddEntry(graph1,"Data","PE");
  leg.AddEntry(scaledgraph,"MC","PE");
  leg.Draw();
  c->Print("output/data_mc_comparison.pdf");
  c->Print("output/data_mc_comparison.cxx");
}
