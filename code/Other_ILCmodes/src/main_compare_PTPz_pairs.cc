#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include <bitset>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "UsefulFunctions.h"

#include "Style.h"

using namespace std;

void Draw_histos(std::vector< TH1D* > histos, TCanvas* canvas);
bool Check_Set_InputArguments(int const argc, char const * const * const argv, 
    std::vector<std::string> *pair500_inputfilenames, std::vector<std::string> *pair350_inputfilenames, std::vector<std::string> *pair250_inputfilenames,
    int &NUMBER_OF_pair500FILES, int &NUMBER_OF_pair350FILES, int &NUMBER_OF_pair250FILES,
    float &weight_bunches_pair500, float &weight_bunches_pair350, float &weight_bunches_pair250);

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  TH1::SetDefaultSumw2();

  std::vector<std::string> *pair500_inputfilenames = new std::vector<std::string>();
  std::vector<std::string> *pair350_inputfilenames = new std::vector<std::string>();
  std::vector<std::string> *pair250_inputfilenames = new std::vector<std::string>();

  int NUMBER_OF_pair500FILES = 0;
  int NUMBER_OF_pair350FILES = 0;
  int NUMBER_OF_pair250FILES = 0;

  float weight_bunches_pair500 = 0;
  float weight_bunches_pair350 = 0;
  float weight_bunches_pair250 = 0;

  if (Check_Set_InputArguments(argc, argv, 
        pair500_inputfilenames, pair350_inputfilenames, pair250_inputfilenames,
        NUMBER_OF_pair500FILES, NUMBER_OF_pair350FILES, NUMBER_OF_pair250FILES,
        weight_bunches_pair500, weight_bunches_pair350, weight_bunches_pair250) == false){
    return -1;
  }

  std::vector<std::string> *inputfiles = new std::vector<std::string>();
  inputfiles->insert(inputfiles->begin(), pair500_inputfilenames->begin(), pair500_inputfilenames->end());
  inputfiles->insert(inputfiles->end(), pair350_inputfilenames->begin(), pair350_inputfilenames->end());
  inputfiles->insert(inputfiles->end(), pair250_inputfilenames->begin(), pair250_inputfilenames->end());
  std::cout << inputfiles->size() << " files will be processed." << std::endl;

  //Make histogram for storing the information
  std::string const histo_title_PzPT_500 = "P_T vs. P_z of pairs from ILC 500GeV;P_z [GeV];P_T [GeV]";
  std::string const histo_title_PzPT_350 = "P_T vs. P_z of pairs from ILC 350GeV;P_z [GeV];P_T [GeV]";
  std::string const histo_title_PzPT_250 = "P_T vs. P_z of pairs from ILC 250GeV;P_z [GeV];P_T [GeV]";
  TH2D* Pz_PT_pairs500 = new TH2D("Pz_PT_pairs500",histo_title_PzPT_500.c_str(),100,-100,100,50,0,1);
  TH2D* Pz_PT_pairs350 = new TH2D("Pz_PT_pairs350",histo_title_PzPT_350.c_str(),100,-100,100,50,0,1);
  TH2D* Pz_PT_pairs250 = new TH2D("Pz_PT_pairs250",histo_title_PzPT_250.c_str(),100,-100,100,50,0,1);
  std::string const histo_title_Pz = "P_z of pairs from ILC 500, 350, and 250GeV;P_z [GeV];# particles";
  TH1D* Pz_pairs500 = new TH1D("Pz_pairs500",histo_title_Pz.c_str(),100,-100,100);
  TH1D* Pz_pairs350 = new TH1D("Pz_pairs350",histo_title_Pz.c_str(),100,-100,100);
  TH1D* Pz_pairs250 = new TH1D("Pz_pairs250",histo_title_Pz.c_str(),100,-100,100);
  std::string const histo_title_PT = "P_T of pairs from ILC 500, 350, and 250GeV;P_T [GeV];# particles";
  TH1D* PT_pairs500 = new TH1D("PT_pairs500",histo_title_PT.c_str(),30,0,1.7);
  TH1D* PT_pairs350 = new TH1D("PT_pairs350",histo_title_PT.c_str(),30,0,1.7);
  TH1D* PT_pairs250 = new TH1D("PT_pairs250",histo_title_PT.c_str(),30,0,1.7);

  for (size_t file_iterator = 0; file_iterator < inputfiles->size(); ++file_iterator) {

    TFile *file = TFile::Open(inputfiles->at(file_iterator).c_str());
    TTree *tree = Get_TTree(file, "MCP");

    //Set the branches
    double mom_x = 0.0;
    double mom_y = 0.0;
    double mom_z = 0.0;
    bool CreatedInSimulation_Status = 0;
    tree->SetBranchStatus("*", 0);

    if (tree->GetName() == std::string("Tree_MCP")){
      tree->SetBranchStatus("CreatedInSimulation_Status", kTRUE);
      tree->SetBranchAddress("CreatedInSimulation_Status", &CreatedInSimulation_Status);
      tree->SetBranchStatus("Momentumx", 1);
      tree->SetBranchAddress("Momentumx", &mom_x);
      tree->SetBranchStatus("Momentumy", 1);
      tree->SetBranchAddress("Momentumy", &mom_y);
      tree->SetBranchStatus("Momentumz", 1);
      tree->SetBranchAddress("Momentumz", &mom_z);
    }

    //Now we loop through the tree
    long long int const entries = tree->GetEntries();
    float weight(1.0);
    if( file_iterator < NUMBER_OF_pair500FILES){
      weight = 1.0/weight_bunches_pair500;
    }
    else if( file_iterator >= NUMBER_OF_pair500FILES && file_iterator < NUMBER_OF_pair500FILES + NUMBER_OF_pair350FILES){
      weight = 1.0/weight_bunches_pair350;
    }
    else if( file_iterator >= NUMBER_OF_pair500FILES + NUMBER_OF_pair350FILES && file_iterator < NUMBER_OF_pair500FILES + NUMBER_OF_pair350FILES + NUMBER_OF_pair250FILES){
      weight = 1.0/weight_bunches_pair250;
    }
    else{
      std::cerr << "The file iterator went outside the range of inputfiles." << std::endl;
      exit(-1);
    }

    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      //if(CreatedInSimulation_Status == 1) continue;
      if( file_iterator < NUMBER_OF_pair500FILES){
        Pz_pairs500->Fill(mom_z,weight);
        PT_pairs500->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight);
        Pz_PT_pairs500->Fill(mom_z,sqrt(mom_x*mom_x+mom_y*mom_y),weight);
      }
      if( file_iterator >= NUMBER_OF_pair500FILES && file_iterator < NUMBER_OF_pair500FILES + NUMBER_OF_pair350FILES){
        Pz_pairs350->Fill(mom_z,weight);
        PT_pairs350->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight);
        Pz_PT_pairs350->Fill(mom_z,sqrt(mom_x*mom_x+mom_y*mom_y),weight);
      }
      if( file_iterator >= NUMBER_OF_pair500FILES + NUMBER_OF_pair350FILES && file_iterator < NUMBER_OF_pair500FILES + NUMBER_OF_pair350FILES + NUMBER_OF_pair250FILES){
        Pz_pairs250->Fill(mom_z,weight);
        PT_pairs250->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight);
        Pz_PT_pairs250->Fill(mom_z,sqrt(mom_x*mom_x+mom_y*mom_y),weight);
      }
    }
    file->Close();
  }
  gStyle->SetOptStat(0);

  //Plot the histogram and save it
  std::vector< TH1D* > PT_histos;
  std::vector< TH1D* > Pz_histos;
  PT_histos.push_back(PT_pairs500);
  PT_histos.push_back(PT_pairs350);
  PT_histos.push_back(PT_pairs250);
  Pz_histos.push_back(Pz_pairs500);
  Pz_histos.push_back(Pz_pairs350);
  Pz_histos.push_back(Pz_pairs250);
  double max_PT=GetMinMaxForMultipleOverlappingHistograms(PT_histos,true).second;
  double max_Pz=GetMinMaxForMultipleOverlappingHistograms(Pz_histos,true).second;
  for(size_t iterator = 0; iterator < PT_histos.size(); ++iterator){
    PT_histos.at(iterator)->SetMinimum(0.5);
    PT_histos.at(iterator)->SetMaximum(max_PT);
    Pz_histos.at(iterator)->SetMinimum(0.5);
    Pz_histos.at(iterator)->SetMaximum(max_Pz);
  }


  TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);

  canvas1->cd();
  canvas1->SetLogy(1);
  canvas1->SetLogz(0);

  Draw_histos(PT_histos,canvas1);
  canvas1->Print("output/pairs_comparison_PT.pdf");
  canvas1->Print("output/pairs_comparison_PT.cxx");
  Draw_histos(Pz_histos,canvas1);
  canvas1->Print("output/pairs_comparison_Pz.pdf");
  canvas1->Print("output/pairs_comparison_Pz.cxx");

  canvas1->SetLogy(0);
  canvas1->SetLogz(1);
  Pz_PT_pairs500->Draw("colz");
  canvas1->Print("output/pairs500_PTvsPz.pdf");
  canvas1->Print("output/pairs500_PTvsPz.cxx");
  Pz_PT_pairs350->Draw("colz");
  canvas1->Print("output/pairs350_PTvsPz.pdf");
  canvas1->Print("output/pairs350_PTvsPz.cxx");
  Pz_PT_pairs250->Draw("colz");
  canvas1->Print("output/pairs250_PTvsPz.pdf");
  canvas1->Print("output/pairs250_PTvsPz.cxx");

  return 0;
}

void Draw_histos(std::vector< TH1D* > histos, TCanvas* canvas){
  std::vector<TPaveStats*> st_vec;
  float boxsize = 0.0;
  std::vector<TPaveText*> text_vec;
  int color = 1;

  for(size_t hist_num =0; hist_num < histos.size(); ++hist_num){
    if (hist_num+color == 3 || hist_num+color == 5) color++;

    histos.at(hist_num)->SetMarkerColor(hist_num+color);
    if(hist_num == 0) histos.at(hist_num)->Draw();
    else histos.at(hist_num)->Draw("SAMES");
    //canvas->Update();
    //st_vec.push_back(new TPaveStats());
    //st_vec.at(hist_num) = (TPaveStats*)histos.at(hist_num)->GetListOfFunctions()->FindObject("stats");
    //st_vec.at(hist_num)->SetLineColor(hist_num+1);
    //st_vec.at(hist_num)->SetX1NDC(0.65); //new x start position
    //st_vec.at(hist_num)->SetX2NDC(0.85); //new x end position
    //if(hist_num == 0){
    //  st_vec.at(hist_num)->SetY1NDC(0.83); //new x start position
    //  st_vec.at(hist_num)->SetY2NDC(0.9); //new x end position
    //  boxsize = st_vec.at(hist_num)->GetY2NDC()-st_vec.at(hist_num)->GetY1NDC();
    //}
    //else{
    //  st_vec.at(hist_num)->SetY2NDC(st_vec.at(hist_num-1)->GetY1NDC()); //new x end position
    //  st_vec.at(hist_num)->SetY1NDC(st_vec.at(hist_num)->GetY2NDC()-boxsize); //new x start position
    //}

    text_vec.push_back( new TPaveText(0.65,0.83-hist_num*0.07,0.85,0.9-hist_num*0.07,"brNDC") );
    text_vec.at(hist_num)->SetTextFont(62);
    text_vec.at(hist_num)->SetTextColor(hist_num+color);
    text_vec.at(hist_num)->SetFillColor(0);
    text_vec.at(hist_num)->AddText(histos.at(hist_num)->GetName());
    text_vec.at(hist_num)->AddLine(0,0.5,1,0.5);
    std::stringstream entriespairs;
    entriespairs << "Entries = " << histos.at(hist_num)->Integral();
    text_vec.at(hist_num)->AddText(entriespairs.str().c_str());
    text_vec.at(hist_num)->Draw();

  }
}

bool Check_Set_InputArguments(int const argc, char const * const * const argv, 
    std::vector<std::string> *pair500_inputfilenames, std::vector<std::string> *pair350_inputfilenames, std::vector<std::string> *pair250_inputfilenames,
    int &NUMBER_OF_pair500FILES, int &NUMBER_OF_pair350FILES, int &NUMBER_OF_pair250FILES,
    float &weight_bunches_pair500, float &weight_bunches_pair350, float &weight_bunches_pair250){
  bool pair500_inputfiles_set = false;
  bool pair350_inputfiles_set = false;
  bool pair250_inputfiles_set = false;
  bool NUMBER_OF_pair500FILES_set = false;
  bool NUMBER_OF_pair350FILES_set = false;
  bool NUMBER_OF_pair250FILES_set = false;
  bool weight_bunches_pair500_set = false;
  bool weight_bunches_pair350_set = false;
  bool weight_bunches_pair250_set = false;
  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-n500")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-n350") 
          && argv[i + 1] != std::string("-n250") 
          && argv[i + 1] != std::string("-w500")
          && argv[i + 1] != std::string("-w350")
          && argv[i + 1] != std::string("-w250")
          && argv[i + 1] != std::string("-i500")
          && argv[i + 1] != std::string("-i350")
          && argv[i + 1] != std::string("-i250")) {
        NUMBER_OF_pair500FILES = std::stoi(argv[i + 1]);
        std::cout << "Number of 500GeV pair input files = " << NUMBER_OF_pair500FILES << std::endl;
        NUMBER_OF_pair500FILES_set = true;
      } else {
        std::cerr << "You didn't give an argument for the number of 500GeV pair files!" << std::endl;
      }
    }
    if (argv[i] == std::string("-n350")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-n500") 
          && argv[i + 1] != std::string("-n250") 
          && argv[i + 1] != std::string("-w500")
          && argv[i + 1] != std::string("-w350")
          && argv[i + 1] != std::string("-w250")
          && argv[i + 1] != std::string("-i500")
          && argv[i + 1] != std::string("-i350")
          && argv[i + 1] != std::string("-i250")) {
        NUMBER_OF_pair350FILES = std::stoi(argv[i + 1]);
        std::cout << "Number of 350GeV pair input files = " << NUMBER_OF_pair350FILES << std::endl;
        NUMBER_OF_pair350FILES_set = true;
      } else {
        std::cerr << "You didn't give an argument for the number of 350GeV pair files!" << std::endl;
      }
    }
    if (argv[i] == std::string("-n250")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-n350") 
          && argv[i + 1] != std::string("-n200") 
          && argv[i + 1] != std::string("-w500")
          && argv[i + 1] != std::string("-w350")
          && argv[i + 1] != std::string("-w250")
          && argv[i + 1] != std::string("-i500")
          && argv[i + 1] != std::string("-i350")
          && argv[i + 1] != std::string("-i250")) {
        NUMBER_OF_pair250FILES = std::stoi(argv[i + 1]);
        std::cout << "Number of 250GeV pair input files = " << NUMBER_OF_pair250FILES << std::endl;
        NUMBER_OF_pair250FILES_set = true;
      } else {
        std::cerr << "You didn't give an argument for the number of 250GeV pair files!" << std::endl;
      }
    }
  }
  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i500") 
        && argv[i + 1] != std::string("-i350") 
        && argv[i + 1] != std::string("-i250") 
        && argv[i + 1] != std::string("-n500") 
        && argv[i + 1] != std::string("-n350") 
        && argv[i + 1] != std::string("-n250") 
        && argv[i + 1] != std::string("-w500") 
        && argv[i + 1] != std::string("-w350") 
        && argv[i + 1] != std::string("-w250")) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-i350") 
              && argv[i + j] != std::string("-i250")
              && argv[i + j] != std::string("-n500") 
              && argv[i + j] != std::string("-n350") 
              && argv[i + j] != std::string("-n250")
              && argv[i + j] != std::string("-w500") 
              && argv[i + j] != std::string("-w350") 
              && argv[i + j] != std::string("-w250")) {
            pair500_inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_pair500FILES);
        pair500_inputfiles_set = true;
      } else {
        std::cerr << "You didn't give an argument for the 500GeV pair inputfile(s)!" << std::endl;
      }
    }
    if (argv[i] == std::string("-i350") 
        && argv[i + 1] != std::string("-i500") 
        && argv[i + 1] != std::string("-i250") 
        && argv[i + 1] != std::string("-n500") 
        && argv[i + 1] != std::string("-n350") 
        && argv[i + 1] != std::string("-n250") 
        && argv[i + 1] != std::string("-w500") 
        && argv[i + 1] != std::string("-w350") 
        && argv[i + 1] != std::string("-w250")) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-i500") 
              && argv[i + j] != std::string("-i250")
              && argv[i + j] != std::string("-n500") 
              && argv[i + j] != std::string("-n350") 
              && argv[i + j] != std::string("-n250")
              && argv[i + j] != std::string("-w500") 
              && argv[i + j] != std::string("-w350") 
              && argv[i + j] != std::string("-w250")) {
            pair350_inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_pair350FILES);
        pair350_inputfiles_set = true;
      } else {
        std::cerr << "You didn't give an argument for the 350GeV pair inputfile(s)!" << std::endl;
      }
    }
    if (argv[i] == std::string("-i250") 
        && argv[i + 1] != std::string("-i500") 
        && argv[i + 1] != std::string("-i350") 
        && argv[i + 1] != std::string("-n500") 
        && argv[i + 1] != std::string("-n350") 
        && argv[i + 1] != std::string("-n250") 
        && argv[i + 1] != std::string("-w500") 
        && argv[i + 1] != std::string("-w350") 
        && argv[i + 1] != std::string("-w250")) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-i500") 
              && argv[i + j] != std::string("-i350")
              && argv[i + j] != std::string("-n500") 
              && argv[i + j] != std::string("-n350") 
              && argv[i + j] != std::string("-n250")
              && argv[i + j] != std::string("-w500") 
              && argv[i + j] != std::string("-w350") 
              && argv[i + j] != std::string("-w250")) {
            pair250_inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_pair250FILES);
        pair250_inputfiles_set = true;
      } else {
        std::cerr << "You didn't give an argument for the 250GeV pair inputfile(s)!" << std::endl;
      }
    }
    if (argv[i] == std::string("-w500")){ 
      if(argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-w500") 
          && argv[i + 1] != std::string("-w350") 
          && argv[i + 1] != std::string("-w250") 
          && argv[i + 1] != std::string("-n500") 
          && argv[i + 1] != std::string("-n350") 
          && argv[i + 1] != std::string("-n250") 
          && argv[i + 1] != std::string("-i500") 
          && argv[i + 1] != std::string("-i350") 
          && argv[i + 1] != std::string("-i250")) {
        std::string str = argv[i + 1];
        weight_bunches_pair500 = std::atof(str.c_str());
        std::cout << "weight_bunches_pair500 = " << weight_bunches_pair500 << std::endl; 
        weight_bunches_pair500_set = true;
      } else {
        std::cerr << "You didn't give an argument for the 500GeV pairs weight!" << std::endl;
        exit(1);
      }
    }
    if (argv[i] == std::string("-w350")){ 
      if(argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-w500") 
          && argv[i + 1] != std::string("-w350") 
          && argv[i + 1] != std::string("-w250") 
          && argv[i + 1] != std::string("-n500") 
          && argv[i + 1] != std::string("-n350") 
          && argv[i + 1] != std::string("-n250") 
          && argv[i + 1] != std::string("-i500") 
          && argv[i + 1] != std::string("-i350") 
          && argv[i + 1] != std::string("-i250")) {
        std::string str = argv[i + 1];
        weight_bunches_pair350 = std::atof(str.c_str());
        std::cout << "weight_bunches_pair350 = " << weight_bunches_pair350 << std::endl; 
        weight_bunches_pair350_set = true;

      } else {
        std::cerr << "You didn't give an argument for the 350GeV pairs weight!" << std::endl;
        exit(1);
      }
    }
    if (argv[i] == std::string("-w250")){ 
      if(argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-w500") 
          && argv[i + 1] != std::string("-w350") 
          && argv[i + 1] != std::string("-n500") 
          && argv[i + 1] != std::string("-n350") 
          && argv[i + 1] != std::string("-n250") 
          && argv[i + 1] != std::string("-i500") 
          && argv[i + 1] != std::string("-i350") 
          && argv[i + 1] != std::string("-i250")) {
        std::string str = argv[i + 1];
        weight_bunches_pair250 = std::atof(str.c_str());
        std::cout << "weight_bunches_pair250 = " << weight_bunches_pair250 << std::endl; 
        weight_bunches_pair250_set = true;
      } else {
        std::cerr << "You didn't give an argument for the 250GeV pairs weight!" << std::endl;
        exit(1);
      }
    }
  }
  if (!pair500_inputfiles_set || !NUMBER_OF_pair500FILES_set || !weight_bunches_pair500_set ||
      !pair350_inputfiles_set || !NUMBER_OF_pair350FILES_set || !weight_bunches_pair350_set ||
      !pair250_inputfiles_set || !NUMBER_OF_pair250FILES_set || !weight_bunches_pair250_set) {
    std::cerr
      << "You didn't give the name for the inputfiles, the amount of pair inputfiles, or the bunch weight. Please try again!"
      << std::endl;
    return(false);
  }
  else return(true);
}
