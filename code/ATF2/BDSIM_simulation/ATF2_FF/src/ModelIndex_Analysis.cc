#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1.h"
#include "TText.h"
#include "TLegend.h"
#include "TGaxis.h"

#include <iostream>

#include "UsefulFunctions.h"
#include "Style.h"

bool Check_Set_InputArguments(int const argc, char const * const * const argv, std::string &in1, std::string &in2, std::string &out);
void DrawTheBarCharts(TH1F* hist, TH1F* hist2, int numbins, std::string outname);

int main(int const argc, char const * const * const argv){
  //UsePhDStyle();

  std::string inputfilename_open;
  std::string inputfilename_closed;
  std::string outputfilename;
  //float aperture = 0.0;
  //float upperjaw = 0.0;
  //float lowerjaw = 0.0;
  
  if (Check_Set_InputArguments(argc, argv, inputfilename_open, inputfilename_closed, outputfilename) == false) return -1;

  TFile* inputfile_open = TFile::Open(inputfilename_open.c_str());
  TTree* Tree_open = nullptr;
  inputfile_open->GetObject("TrackModelIndex", Tree_open);
  TFile* inputfile_closed = TFile::Open(inputfilename_closed.c_str());
  TTree* Tree_closed = nullptr;
  inputfile_closed->GetObject("TrackModelIndex", Tree_closed);

  //Set the branches
  int modelindex_open = 0;
  Tree_open->SetBranchAddress("modelIndex", &modelindex_open);
  int modelindex_closed = 0;
  Tree_closed->SetBranchAddress("modelIndex", &modelindex_closed);

    //Make the histogram
  TH1F* TracksPerModel_openCollimator1 = new TH1F("TracksPerModel_openColl1","Number of Particles created in the ATF2 beam line components",56,0,56);
  TH1F* TracksPerModel_openCollimator2 = new TH1F("TracksPerModel_openColl2","Number of Particles created in the ATF2 beam line components",56,56,112);
  TH1F* TracksPerModel_openCollimator3 = new TH1F("TracksPerModel_openColl3","Number of Particles created in the ATF2 beam line components",56,112,168);
  TH1F* TracksPerModel_openCollimator4 = new TH1F("TracksPerModel_openColl4","Number of Particles created in the ATF2 beam line components",55,168,223);
  std::vector< TH1F* > TracksPerModel_openCollimator;
  TracksPerModel_openCollimator.push_back(TracksPerModel_openCollimator1);
  TracksPerModel_openCollimator.push_back(TracksPerModel_openCollimator2);
  TracksPerModel_openCollimator.push_back(TracksPerModel_openCollimator3);
  TracksPerModel_openCollimator.push_back(TracksPerModel_openCollimator4);
  for(int t = 0; t < TracksPerModel_openCollimator.size(); ++t){
    TracksPerModel_openCollimator.at(t)->SetFillColor(kTeal-6);
    TracksPerModel_openCollimator.at(t)->SetBarWidth(0.4);
    TracksPerModel_openCollimator.at(t)->SetBarOffset(0.1);
    TracksPerModel_openCollimator.at(t)->SetStats(0);
    TracksPerModel_openCollimator.at(t)->SetMinimum(0.01);
    if(t = 0) TracksPerModel_openCollimator.at(t)->SetMaximum(18*std::pow(10,6));
    if(t = 1) TracksPerModel_openCollimator.at(t)->SetMaximum(1*std::pow(10,6));
    if(t = 2) TracksPerModel_openCollimator.at(t)->SetMaximum(5*std::pow(10,4));
    if(t = 3) TracksPerModel_openCollimator.at(t)->SetMaximum(18*std::pow(10,6));
    TracksPerModel_openCollimator.at(t)->GetXaxis()->SetLabelOffset(99);
  }
  TH1F* TracksPerModel_closedCollimator1 = new TH1F("TracksPerModel_closedColl1","Number of Particles created in the ATF2 beam line components",56,0,56);
  TH1F* TracksPerModel_closedCollimator2 = new TH1F("TracksPerModel_closedColl2","Number of Particles created in the ATF2 beam line components",56,56,112);
  TH1F* TracksPerModel_closedCollimator3 = new TH1F("TracksPerModel_closedColl3","Number of Particles created in the ATF2 beam line components",56,112,168);
  TH1F* TracksPerModel_closedCollimator4 = new TH1F("TracksPerModel_closedColl4","Number of Particles created in the ATF2 beam line components",55,168,223);
  std::vector< TH1F* > TracksPerModel_closedCollimator;
  TracksPerModel_closedCollimator.push_back(TracksPerModel_closedCollimator1);
  TracksPerModel_closedCollimator.push_back(TracksPerModel_closedCollimator2);
  TracksPerModel_closedCollimator.push_back(TracksPerModel_closedCollimator3);
  TracksPerModel_closedCollimator.push_back(TracksPerModel_closedCollimator4);
  for(int t = 0; t < TracksPerModel_closedCollimator.size(); ++t){
    TracksPerModel_closedCollimator.at(t)->SetFillColor(kPink+4);
    TracksPerModel_closedCollimator.at(t)->SetBarWidth(0.4);
    TracksPerModel_closedCollimator.at(t)->SetBarOffset(0.5);
    TracksPerModel_closedCollimator.at(t)->SetStats(0);
    TracksPerModel_closedCollimator.at(t)->GetXaxis()->SetLabelOffset(99);
  }

  int entries_open = Tree_open->GetEntries();
  int entries_closed = Tree_closed->GetEntries();

  // Now iterate through the TTree entries and fill a histogram.
  for(int i = 0; i < entries_open; ++i) {
    Tree_open->GetEntry(i);
    for(int t = 0; t < TracksPerModel_openCollimator.size(); ++t){
      TracksPerModel_openCollimator.at(t)->Fill( modelindex_open );
    }
  } // TTree entry / event loop
  for(int i = 0; i < entries_closed; ++i) {
    Tree_closed->GetEntry(i);
    for(int t = 0; t < TracksPerModel_closedCollimator.size(); ++t){
      TracksPerModel_closedCollimator.at(t)->Fill( modelindex_closed, entries_open/entries_closed );
    }
  } // TTree entry / event loop


  DrawTheBarCharts(TracksPerModel_openCollimator.at(0),TracksPerModel_closedCollimator.at(0),0,"output/"+outputfilename+"_firstPart");
  DrawTheBarCharts(TracksPerModel_openCollimator.at(1),TracksPerModel_closedCollimator.at(1),1,"output/"+outputfilename+"_secondPart");
  DrawTheBarCharts(TracksPerModel_openCollimator.at(2),TracksPerModel_closedCollimator.at(2),2,"output/"+outputfilename+"_thirdPart");
  DrawTheBarCharts(TracksPerModel_openCollimator.at(3),TracksPerModel_closedCollimator.at(3),3,"output/"+outputfilename+"_forthPart");

  inputfile_open->Close();
  inputfile_closed->Close();
  return 0;
}

void DrawTheBarCharts(TH1F* hist, TH1F* hist2, int numbins, std::string outname){
  TCanvas* canvas = new TCanvas();
  gStyle->SetOptStat(0);
  canvas->SetGrid();
  canvas->SetBottomMargin(0.15);
  TGaxis::SetMaxDigits(4);
  hist->SetFillColor(kTeal-6);
  hist->SetBarWidth(0.4);
  hist->SetBarOffset(0.1);
  hist->SetStats(0);
  hist->GetXaxis()->SetLabelOffset(99);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->SetMinimum(0.01);
  if(numbins == 0){
    hist->SetMaximum(18*std::pow(10,6));
  }
  if(numbins == 1){
    hist->SetMaximum(7*std::pow(10,5));
  }
  if(numbins == 2){
    hist->SetMaximum(3*std::pow(10,4));
  }
  if(numbins == 3){
    hist->SetMaximum(18*std::pow(10,6));
  }
  hist2->SetFillColor(kPink+4);
  hist2->SetBarWidth(0.4);
  hist2->SetBarOffset(0.5);
  hist2->SetStats(0);
  hist2->GetXaxis()->SetLabelOffset(99);
  hist2->GetYaxis()->SetNoExponent(kFALSE);
  hist->Draw("b");
  hist2->Draw("b,same");


  TLegend* leg = new TLegend(0.25,0.6,0.55,0.8);
  leg->AddEntry(hist,"Open Vertical Collimator","f");
  leg->AddEntry(hist2,"Closed Vertical Collimator","f");
  leg->Draw();

  int nx = 56;
  const Int_t nxall = 223;
  char *models[nxall] = {
      "L200","QM16FF","QM16FF_1","L201A","ZH1FF", "L201B","ZV1FF","L201C","QM15FF","QM15FF_1","L202A","MQM15FF","L202B","QM14FF","QM14FF_1","L203A","MQM14FF",
       "L203B","L203C","L203D","L203E","MFB2FF","L203F","QM13FF","QM13FF_1","L204A","MQM13FF","L204B","QM12FF","QM12FF_1","L205A","MQM12FF","L205B","L205C",
      "QM11FF","QM11FF_1","L206A","MQM11FF","L206BA","COLLBY","L206BB","QD10BFF","QD10BFF_1","L207A","MQD10BFF","COLL","L207B","CREF3","L207C","QD10AFF",
   "QD10AFF_1","L208A","MQD10AFF","L208B","L208C","QF9BFF","QF9BFF_1","L209A","MQF9BFF","L209B","SF6FF","SF6FF_1","L210A","MSF6FF","L210B","QF9AFF","QF9AFF_1",
       "L211A","MQF9AFF","L211B","SK4FF","SK4FF_1","L211C0","L211C","QD8FF","QD8FF_1","L212A","MQD8FF","L212B","QF7FF","QF7FF_1","L213A","MQF7FF","L213B",
   "B5FFA_eve","B5FFA_eve","B5FFA_eve","B5FFB_eve","B5FFB_eve","B5FFB_eve","L214","QD6FF","QD6FF_1","L215A","MQD6FF","L215B","SK3FF","SK3FF_1","L215C","QF5BFF",
   "QF5BFF_1","L216A","MQF5BFF","L216B","SF5FF","SF5FF_1","L217A","MSF5FF","L217B","QF5AFF","QF5AFF_1","L218A","MQF5AFF","L218B","L218C","QD4BFF","QD4BFF_1",
       "L219A","MQD4BFF","L219B","SD4FF","SD4FF_1","L220A","MSD4FF","L220B","QD4AFF","QD4AFF_1","L221A","MQD4AFF","L221B","SK2FF","SK2FF_1","L221C","B2FFA_eve",
   "B2FFA_eve","B2FFA_eve","B2FFB_eve","B2FFB_eve","B2FFB_eve","L222","QF3FF","QF3FF_1","L223A","MQF3FF","L223B","B1FFA_eve","B1FFA_eve","B1FFA_eve","B1FFB_eve",
   "B1FFB_eve","B1FFB_eve","L224A","L224B","QD2BFF","QD2BFF_1","L225A","MQD2BFF","L225B","L225C","L225D","L225E","QD2AFF","QD2AFF_1","L226A","MQD2AFF","L226B",
       "L226C","SK1FF","SK1FF_1","L226D","L226E","SF1FF","SF1FF_1","MSF1FF","L227","QF1FF","QF1FF_1","L228A","L228B","L228C","SD0FF","SD0FF_1","MSD0FF","L229",
       "QD0FF","QD0FF_1","L230A","MPREIP","L230B","IPKICK","L230C","L230D","L230E","L230F","L231","L301A","L301B","L301C","L301D","BDUMPA_ev","BDUMPA_ev",
   "BDUMPA_ev","BDUMPA_ev","BDUMPA_ev","BDUMPA_ev","BDUMPA_ev","BDUMPB_ev","BDUMPB_ev","BDUMPB_ev","BDUMPB_ev","BDUMPB_ev","BDUMPB_ev","BDUMPB_ev","BDUMPC_ev",
   "BDUMPC_ev","BDUMPC_ev","BDUMPC_ev","BDUMPC_ev","BDUMPC_ev","BDUMPC_ev","L302A","L302B","L302C"
  };

  //h->GetXaxis()->SetLabelOffset(99);
  // draw labels along X
  Float_t x, y;
  y = gPad->GetUymin() - 0.2*hist->GetYaxis()->GetBinWidth(1);
  TText t;
  t.SetTextAngle(90);
  t.SetTextSize(0.03);
  t.SetTextAlign(33);
  int start = nx*numbins;
  int end = nx*(numbins+1);
  if (numbins == 3) end = nx*(numbins+1)-1;
  int j = 0;
  for (int i = start; i < end; i++) {
    x = hist->GetXaxis()->GetBinCenter(j+1);
    t.DrawText(x,y,models[i]);
    ++j;
  }

  canvas->Print( (outname+".pdf").c_str() );
  canvas->Print( (outname+".cxx").c_str() );
}
bool Check_Set_InputArguments(int const argc, char const * const * const argv, std::string &in1, std::string &in2, std::string &out)
{
  bool inputfilename1_set = false;
  bool inputfilename2_set = false;
  bool outputfilename_set = false;
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i1")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-o")
					&& argv[i + 1] != std::string("-i2")){
				in1 = argv[i + 1];
				inputfilename1_set = true;
			} else {
				std::cerr << "You didn't give an argument for the open collimator inputfilename!"
					<< std::endl;
			}
		}
    if (argv[i] == std::string("-i2")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-o")
					&& argv[i + 1] != std::string("-i1")){
				in2 = argv[i + 1];
				inputfilename2_set = true;
			} else {
				std::cerr << "You didn't give an argument for the closed collimator inputfilename!"
					<< std::endl;
			}
		}
		if (argv[i] == std::string("-o")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-i1")
					&& argv[i + 1] != std::string("-i2")) {
				out = argv[i + 1];
				outputfilename_set = true;
			} else {
				std::cerr << "You didn't give an argument for the outputfilename!"
					<< std::endl;
			}
		}

	}
	if(!inputfilename1_set || !outputfilename_set || !inputfilename2_set){
		std::cerr << "You didn't set the input arguments correctly!" << std::endl;
		std::cerr << "Try:\n";
		std::cerr << "./CollimatorBkgLevelAna -i1 inputfilename_openColl -i2 inputfilename_closedColl -o outputfilename";
		std::cerr << std::endl;
		return false;
	}
  else return true;
}
