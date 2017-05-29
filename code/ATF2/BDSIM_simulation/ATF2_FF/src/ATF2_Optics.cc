#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "Style.h"
#include "TLegend.h"
#include "TTree.h"
#include "TH2D.h"
#include "TFile.h"

using namespace std;

int main(int const argc, char const * const * const argv){
  //Use PhD Style
  //Filename is an option (CoreOptics_10032017_event_OpticsOutput.root)
  //TTree name is optics
  //Sigma_x
  //"_y
  //Beta_x
  //"_y
  //All of these plotted against S (x-axis)
  //Plots should have legend
  //Thick points so you can see the distribution
  //Scatter plot

  UsePhDStyle();

  string filename("");

  try{
    if(argc <=2) throw;
    for(int i = 0; i < argc; ++i){
      if(string(argv[i]) == "-i") filename = argv[i+1];
    }
  } catch(...){
    cerr << "Input arguments were invalid" << endl;
    exit(1);
  }

  shared_ptr< TFile > input_file(TFile::Open(filename.c_str()));
  if(input_file->IsZombie()){
    cerr << "Failed to open the TFile" << endl;
    exit(2);
  }

  shared_ptr< TTree > tree(static_cast< TTree* >(input_file->Get("optics")));
  
  shared_ptr< TH2D > histo_sigmax_s(new TH2D("histo_sigmax_s",";The path along the beam line [m];Beam size [mm]",45,0,45,50,0,0.0020));
  int entries_sigmax_s = tree->Draw("Sigma_x:S >>+ histo_sigmax_s","","p");
  cout << "Histogram histo_sigmax_s has " << entries_sigmax_s << " entries" << endl;
  
  shared_ptr< TH2D > histo_sigmay_s(new TH2D("histo_sigmay_s",";The path along the beam line [m];Beam size [mm]",45,0,45,50,0,0.0020));
  int entries_sigmay_s = tree->Draw("Sigma_y:S >>+ histo_sigmay_s","","p");
  cout << "Histogram histo_sigmay_s has " << entries_sigmay_s << " entries" << endl;

  shared_ptr< TH2D > histo_betax_s(new TH2D("histo_betax_s",";The path along the beam line [m];#beta Function [mm^{2}]",45,0,45,120,0,4500));
  int entries_betax_s = tree->Draw("Beta_x:S >>+ histo_betax_s","","p");
  cout << "Histogram histo_betax_s has " << entries_betax_s << " entries" << endl;
  
  shared_ptr< TH2D > histo_betay_s(new TH2D("histo_betay_s",";The path along the beam line [m];#beta Function [mm^{2}]",45,0,45,120,0,4500));
  int entries_betay_s = tree->Draw("Beta_y:S >>+ histo_betay_s","","p");
  cout << "Histogram histo_betay_s has " << entries_betay_s << " entries" << endl;

  histo_sigmax_s->SetMarkerStyle(8);
  histo_sigmax_s->SetMarkerColor(kPink-6);
  histo_sigmay_s->SetMarkerStyle(8);
  histo_sigmay_s->SetMarkerColor(kTeal+3);
  histo_betax_s->SetMarkerStyle(8);
  histo_betax_s->SetMarkerColor(kPink-6);
  histo_betay_s->SetMarkerStyle(8);
  histo_betay_s->SetMarkerColor(kTeal+3);

  TLegend legsigma(0.35,0.5,0.45,0.75), legbeta(0.35,0.5,0.45,0.75);
  legsigma.AddEntry(histo_sigmax_s.get(),"#sigma_{x}","p");
  legsigma.AddEntry(histo_sigmay_s.get(),"#sigma_{y}","p");
  legbeta.AddEntry(histo_betax_s.get(),"#beta_{x}","p");
  legbeta.AddEntry(histo_betay_s.get(),"#beta_{y}","p");
  legbeta.SetTextSize(0.07);

  TCanvas c;
  cout << "histo_sigmax_s has " << histo_sigmax_s->GetEntries() << " entries" << endl;
  //cout << "Values are:" << endl;
  //for(int i = 0; i < histo_sigmax_s->GetNbinsX() + 1; ++i){
  //  cout << histo_sigmax_s->GetBinContent(i) << endl;
  //}

  gStyle->SetOptStat(0);

  histo_sigmax_s->Draw("P");
  histo_sigmay_s->Draw("PSAME");
  legsigma.Draw();
  //c.SetLogy();
  c.Print("output/histo_sigma_s.pdf");
  c.Print("output/histo_sigma_s.cxx");
 
  TCanvas d;
  histo_betax_s->Draw("P");
  histo_betay_s->Draw("PSAME");
  legbeta.Draw();
  d.Print("output/histo_beta_s.pdf");
  d.Print("output/histo_beta_s.cxx");

  return 0;
}
