#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"


#include <bitset>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

#include "UsefulFunctions.h"
#include "GeneralFunctions_SiDBkgSim_new.h"
#include "CellHits_class_new.h"
#include "Subdetector_class_new.h"

using namespace std;

void Draw_All_Layer_plots_together ( std::vector< TH1D* > histo, TCanvas* canvas, bool normalize);
void Draw_multiple_plots (int num_layers, int start_layer, std::vector< TH1D* > histos, TCanvas* canvas, bool normalize);
void Print_multiple_plots_from_same_vec (int num_layers, std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::string output);

double weight = 0.08846; //The weight is the same for both scenarios, for both, the electron and the positron line, e.g.: 5sp+wall, elec: (4321/10155 * 898.34)/4321

int main(int const argc, char const * const * const argv) {
  std::vector<std::string> *inputfilenames = new std::vector<std::string>();
  std::string argument_subdetectors;

  int NUMBER_OF_FILES = 0;
  bool NUMBER_OF_FILES_set = false;
  bool inputfile_set = false;
  bool subdetector_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-s") && 
          argv[i + 1] != std::string("-i")) {
        NUMBER_OF_FILES = std::stoi(argv[i + 1]);
        std::cout << "Number of input files = " << NUMBER_OF_FILES << std::endl;
        NUMBER_OF_FILES_set = true;
      } else {
        std::cerr << "You didn't give an argument for the number of files!" << std::endl;
      }
    }
  }
  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")) {
            inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_FILES*2);//2 files per scenario
        inputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
      }
    }
    if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL && 
          argv[i + 1] != std::string("-n") && 
          argv[i + 1] != std::string("-s")) {
        argument_subdetectors = argv[i + 1];
        subdetector_set = true;
      } else {
        std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
      }
    }
  }
  if (!inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set) {
    std::cerr
      << "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
      << std::endl;
    exit(1);
  }


  Subdetector det( argument_subdetectors );
  std::string subdetectorname = det.getName();
  int max_num_layers = det.getNumberOfLayers();

  bool endcap = false;
  bool barrel = false;
  if (det.getShape().find("endcap") != std::string::npos){//If Endcap, calculate CellID with x and y
    std::cout <<"Found 'endcap'!" << std::endl;
    endcap = true;
  }
  else if (det.getShape().find("barrel") != std::string::npos){
    std::cout <<"Found 'barrel'!" << std::endl;
    barrel = true;
  }
  else{
    std::cerr << "This subdetector shape was not recognized!" << std::endl; 
    exit(-1);
  }
  bool Silicon = false;
  bool Calo = false;
  if (det.getName().find("Si") != std::string::npos){
    std::cout << "Silicon detector found!" << std::endl;
    Silicon = true;
  }
  else if (det.getName().find("cal",0) != std::string::npos
      || det.getName().find("Cal",0) != std::string::npos
      || det.getName().find("Muon",0) != std::string::npos){
    std::cout << "Calorimeter found!" << std::endl;
    Calo = true;
  }
  else{
    std::cerr << "This subdetector name was not recognized!" << std::endl; 
    exit(-1);
  }

  //Make histogram vectors for storing the histograms
  std::vector < std::string > all_name;
  all_name.emplace_back( "All_layers_5sp" );
  all_name.emplace_back( "All_layers_5sp_wall" );
  std::vector < std::string > title;
  title.emplace_back( "Hit time for " + subdetectorname + ";Hit time [ns];Number of hits" );
  title.emplace_back( "Hit time for " + subdetectorname +";Hit time [ns];Number of hits" );
  std::vector< TH1D* > All_Layers_histo;
  All_Layers_histo
  std::vector< TH1D* > histos;

  for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES*2; ++file_iterator) {
    TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, det.getName());

    //Set the branches
    int pdg(0);
    float hittime(0.0);

    tree->SetBranchStatus("*", 0);

    if (Silicon){
      tree->SetBranchStatus("HitParticle_PDG", 1);
      tree->SetBranchAddress("HitParticle_PDG", &pdg);

      tree->SetBranchStatus("HitTime", 1);
      tree->SetBranchAddress("HitTime", &hittime);
    }
    else if (Calo){
      tree->SetBranchStatus("HitMotherParticle_PDG", 1);
      tree->SetBranchAddress("HitMotherParticle_PDG", &pdg);

      tree->SetBranchStatus("HitContrTime", 1);
      tree->SetBranchAddress("HitContrTime", &hittime);
    }
    else{
      std::cerr << "This subdetector name was not recognized!" << std::endl; 
      exit(-1);
    }


    //Now we loop through the tree
    //Combine the two Cell ID's into a single new Cell ID
    //See how often the new Cell ID occurs in total, this is the occupancy

    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      //if (pdg != 13 && pdg != -13) continue;
      //if (endcap && *HitPosition_z < 0) continue;

      All_Layers_histo.at(file_iterator/NUMBER_OF_FILES)->Fill(hittime,weight);  
    }
    file->Close();
  }


  //Plot the histogram and save it
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetLogy(1);

  //std::stringstream output;
  //output << "output/muon_occupancy_" << subdetectorname;
  //Print_multiple_plots_from_same_vec (max_num_layers, histos, canvas, false, output.str());

  Draw_All_Layer_plots_together( All_Layers_histo,canvas, false); 
  std::stringstream All_output;
  All_output << "output/muon_occupancy_all_layers_" << subdetectorname;
  canvas->Print((All_output.str() + ".pdf").c_str());
  canvas->Print((All_output.str() + ".cxx").c_str());

  return 0;
}


void Draw_All_Layer_plots_together ( std::vector< TH1D* > histo, TCanvas* canvas, bool normalize){
  int tot = 0;
  for(int bin = 1; bin <= histo.at(0)->GetNbinsX(); ++bin){
    tot += histo.at(0)->GetBinContent(bin);
  }
  int tot2 = 0;
  for(int bin = 1; bin <= histo.at(1)->GetNbinsX(); ++bin){
    tot2 += histo.at(1)->GetBinContent(bin);
  }
  for (size_t vec_entry = 0; vec_entry < histo.size(); ++vec_entry){
    histo.at(vec_entry)->SetStats(0);
    histo.at(vec_entry)->Sumw2(1);
    //histo.at(vec_entry)->Scale(weight);
  }
  double max=GetMinMaxForMultipleOverlappingHistograms(histo,true).second;
  for (size_t vec_entry = 0; vec_entry < histo.size(); ++vec_entry){
    histo.at(vec_entry)->SetMinimum(max);
    histo.at(vec_entry)->SetMinimum(0.1);
    if(normalize == true){
      histo.at(vec_entry)->Scale(1.0/histo.at(vec_entry)->GetBinContent(1));
      histo.at(vec_entry)->SetMinimum( pow(10,-12) );
    }
    if(vec_entry == 0){
      histo.at(vec_entry)->SetLineColor(2);
      histo.at(vec_entry)->Draw("hist");
      canvas->Update();
    }
    else{
      histo.at(vec_entry)->SetLineColor(4);
      histo.at(vec_entry)->Draw("hist,SAMES");
      canvas->Update();
    }
  }
  TPaveText *text1 = new TPaveText(0.75,0.8,0.95,0.9,"brNDC");
  text1->SetTextFont(62);
  text1->SetTextColor(2);
  text1->SetFillColor(0);
  text1->AddText("5 spoilers");
  text1->AddLine(0,0.5,1,0.5);
  std::stringstream entries_sp;
  entries_sp << "Entries = " << tot;
  text1->AddText(entries_sp.str().c_str());
  TPaveText *text2 = new TPaveText(0.75,0.7,0.95,0.8,"brNDC");
  text2->SetTextFont(62);
  text2->SetTextColor(4);
  text2->SetFillColor(0);
  //text->AddLine(0,0.5,1,0.5);
  text2->AddText("5 spoilers + wall");
  text2->AddLine(0,0.5,1,0.5);
  std::stringstream entries_spwall;
  entries_spwall << "Entries = " << tot2;
  text2->AddText(entries_spwall.str().c_str());
  text1->Draw();
  text2->Draw();
}

void Print_multiple_plots_from_same_vec (int num_layers, std::vector< TH1D* > histos, TCanvas* canvas, bool normalize,  std::string output){
  int start_layer = 0;
  Draw_multiple_plots(num_layers, start_layer, histos, canvas, normalize);
  std::stringstream output1;
  output1 << output << "_5spoilers";
  canvas->Print((output1.str() + ".pdf").c_str());
  canvas->Print((output1.str() + ".cxx").c_str());


  start_layer = num_layers;
  Draw_multiple_plots(num_layers, start_layer, histos, canvas, normalize);
  std::stringstream output2;
  output2 << output << "_5spoilers_wall";
  canvas->Print((output2.str() + ".pdf").c_str());
  canvas->Print((output2.str() + ".cxx").c_str());
}
void Draw_multiple_plots (int num_layers, int start_layer, std::vector< TH1D* > histos, TCanvas* canvas, bool normalize ){
  int i = 0;
  int color = 2; // Very first histogram will be drawn with the color 2, then counted up
  int marker = 20; // Very first histogram will be drawn with the marker style 20, then counted up
  for (int number_histo = start_layer; number_histo< start_layer+num_layers; ++number_histo) {
    histos.at(number_histo)->SetStats(0);
    histos.at(number_histo)->Sumw2(1);
    //histos.at(number_histo)->Scale(weight);
  }
  double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
  for (int number_histo = start_layer; number_histo< start_layer+num_layers; ++number_histo) {
    if(number_histo == start_layer){
      int tot = 0;
      for(int bin = 1; bin <= histos.at(number_histo)->GetNbinsX(); ++bin){
        tot += histos.at(number_histo)->GetBinContent(bin);
      }
      if(normalize == true){
        histos.at(number_histo)->Scale(1.0/histos.at(number_histo)->GetBinContent(1));
        histos.at(number_histo)->SetMinimum( pow(10,-12) );
      }
      else{
        histos.at(number_histo)->SetMaximum(max);
        histos.at(number_histo)->SetMinimum(0.1);
      }
      histos.at(number_histo)->SetLineColor(color);
      histos.at(number_histo)->SetMarkerColor(color);
      histos.at(number_histo)->SetMarkerStyle(marker);
      histos.at(number_histo)->Draw("P");
      canvas->Update();
      TPaveText *text1;
      if (num_layers-start_layer > 5){
        text1 = new TPaveText(0.6,0.8,0.8,0.9,"brNDC");
      }
      else{
        text1 = new TPaveText(0.75,0.8,0.95,0.9,"brNDC");
      }
      text1->SetTextFont(62);
      text1->SetTextColor(color);
      text1->SetFillColor(0);
      std::stringstream title;
      title << histos.at(number_histo)->GetName();
      text1->AddText(title.str().c_str());
      text1->AddLine(0,0.5,1,0.5);
      std::stringstream entries1;
      entries1 << "Entries = " << tot; 
      text1->AddText(entries1.str().c_str());
      text1->Draw();
      i=1;
    }
    if(number_histo > start_layer){
      int tot = 0;
      for(int bin = 1; bin <= histos.at(number_histo)->GetNbinsX(); ++bin){
       tot += histos.at(number_histo)->GetBinContent(bin);
      }
      if(normalize == true){
        histos.at(number_histo)->Scale(1.0/histos.at(number_histo)->GetBinContent(1));
        histos.at(number_histo)->SetMinimum( pow(10,-12) );
      }
      else{
        histos.at(number_histo)->SetMaximum(max);
        histos.at(number_histo)->SetMinimum(0.1);
      }
      color++;
      marker++;
      if(color == 5 || color == 10) color += 1; // 5 would be yellow, 10 would be very light gray 
      histos.at(number_histo)->SetLineColor(color);
      histos.at(number_histo)->SetMarkerColor(color);
      histos.at(number_histo)->SetMarkerStyle(marker);
      histos.at(number_histo)->Draw("P,SAMES");
      canvas->Update();
      TPaveText *text2;
      if (num_layers > 5) {
        if(number_histo >= 5+start_layer){
          if(number_histo == 5+start_layer) {
            i=0;
          }
          text2 = new TPaveText(0.75,0.8-i*0.1,0.95,0.9-i*0.1,"brNDC");
        }
        else {
          text2 = new TPaveText(0.6,0.8-i*0.1,0.8,0.9-i*0.1,"brNDC");
        }
      } 
      else{
        text2 = new TPaveText(0.75,0.8-i*0.1,0.95,0.9-i*0.1,"brNDC");
      }
      text2->SetTextFont(62);
      text2->SetTextColor(color);
      text2->SetFillColor(0);
      //text->AddLine(0,0.5,1,0.5);
      std::stringstream title2;
      title2 << histos.at(number_histo)->GetName();
      text2->AddText(title2.str().c_str());
      text2->AddLine(0,0.5,1,0.5);
      std::stringstream entries2;
      entries2 << "Entries = " << tot;//-1 because of entries in first bin count 1 because of setbincontent
      text2->AddText(entries2.str().c_str());
      text2->Draw();
      i++;
    }
  }
}
