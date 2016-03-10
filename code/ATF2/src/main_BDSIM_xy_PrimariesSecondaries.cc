#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TROOT.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "UsefulFunctions.h"
#include "Style.h"

using namespace std;
void Draw_Print_Beam(TH1* histo, TCanvas* canvas, std::string outputname, std::string PDFfilename);

int main(int const argc, char const * const * const argv) {
  UsePhDStyle();
  //The input is a TTree ROOT file(s)
  //The output is .pdf and .C files

  std::vector < std::string > *inputfilenames = new std::vector<std::string>();
  std::string sample;
  std::string outputfilename = "xyBeamDistribution_";

  int NUMBER_OF_FILES = 0;
  bool NUMBER_OF_FILES_set = false;
  bool inputfile_set = false;
  bool sample_set = false;

  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-n")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-s")
          && argv[i + 1] != std::string("-o")
          && argv[i + 1] != std::string("-i")) {
        NUMBER_OF_FILES = std::stoi(argv[i + 1]);
        std::cout << "Number of input files = " << NUMBER_OF_FILES
          << std::endl;
        NUMBER_OF_FILES_set = true;
      } else {
        std::cerr
          << "You didn't give an argument for the number of files!"
          << std::endl;
      }
    }
  }
  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-i") && argv[i + 1] != std::string("-n")
        && argv[i + 1] != std::string("-s")) {
      if (argv[i + 1] != NULL) {
        int j = 1;
        do {
          if (argv[i + j] != std::string("-n")
              && argv[i + j] != std::string("-s")
              && argv[i + j] != std::string("-o")) {
            inputfilenames->push_back(argv[i + j]);
            ++j;
          } else {
            break;
          }
        } while (j <= NUMBER_OF_FILES);
        inputfile_set = true;
      } else {
        std::cerr << "You didn't give an argument for the inputfile(s)!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-s")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-o")
          && argv[i + 1] != std::string("-n")
          && argv[i + 1] != std::string("-i")) {
        sample = argv[i + 1];
        sample_set = true;
      } else {
        std::cerr << "You didn't give an argument for the sample!"
          << std::endl;
      }
    }
    if (argv[i] == std::string("-o")) {
      if (argv[i + 1] != NULL 
          && argv[i + 1] != std::string("-s")
          && argv[i + 1] != std::string("-n")
          && argv[i + 1] != std::string("-i")) {
        outputfilename = argv[i + 1];
      } else {
        std::cerr << "You didn't give an argument for the sample!"
          << std::endl;
      }
    }

  }
  if (!inputfile_set || !sample_set || !NUMBER_OF_FILES_set) {
    std::cerr
      << "You didn't give the name for the sample, the inputfiles or the amount of files. Please try again!"
      << std::endl;
    exit(1);
  }

  //Plot the histogram and save it
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  //Make histogram for storing the information
  std::string const title1 = "xy distribution of primaries;x [m];y [m]";
  std::string const xtitle1 = "x distribution of primaries;x [m];Count";
  std::string const ytitle1 = "y distribution of primaries;y [m];Count";
  std::string const title2 = "xy distribution of secondaries;x [m];y [m]";
  std::string const xtitle2 = "x distribution of secondaries;x [m];Count";
  std::string const ytitle2 = "y distribution of secondaries;y [m];Count";
  //TH2D *histo_xy_Primaries = new TH2D("xy_Primaries", title1.c_str(), 1000, -0.5, 0.5, 50, -0.05, 0.05);
  //TH1D *histo_x_Primaries = new TH1D("x_Primaries", xtitle1.c_str(), 1000, -0.5, 0.5);
  //TH1D *histo_y_Primaries = new TH1D("y_Primaries", ytitle1.c_str(), 50, -0.05, 0.05);
  //TH2D *histo_xy_Secondaries = new TH2D("xy_Primaries", title2.c_str(), 1000, -0.5, 0.5, 50, -0.05, 0.05);
  //TH1D *histo_x_Secondaries = new TH1D("x_Secondaries", xtitle2.c_str(), 1000, -0.5, 0.5);
  //TH1D *histo_y_Secondaries = new TH1D("y_Secondaries", ytitle2.c_str(), 50, -0.05, 0.05);
  TH2D *histo_xy_Primaries = new TH2D("xy_Primaries", title1.c_str(), 100, -0.0005, 0.0005, 50, -0.000005, 0.000005);
  TH1D *histo_x_Primaries = new TH1D("x_Primaries", xtitle1.c_str(), 100, -0.0005, 0.0005);
  TH1D *histo_y_Primaries = new TH1D("y_Primaries", ytitle1.c_str(), 50, -0.000005, 0.000005);
  TH2D *histo_xy_Secondaries = new TH2D("xy_Primaries", title2.c_str(), 100, -0.0005, 0.0005, 50, -0.000005, 0.000005);
  TH1D *histo_x_Secondaries = new TH1D("x_Secondaries", xtitle2.c_str(), 100, -0.0005, 0.0005);
  TH1D *histo_y_Secondaries = new TH1D("y_Secondaries", ytitle2.c_str(), 50, -0.000005, 0.000005);

  std::string Sampler_Tree_name = "Sampler_" + sample;
  std::cout << Sampler_Tree_name << std::endl;

  for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES;
      ++file_iterator) {
    TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
    TTree *tree = nullptr;
    file->GetObject(Sampler_Tree_name.c_str(), tree);

    //Set the branches
    tree->SetBranchStatus("*", 0);
    //tree->SetBranchStatus("partID", 1);
    tree->SetBranchStatus("parentID", 1);
    tree->SetBranchStatus("x", 1);
    tree->SetBranchStatus("y", 1);
    //tree->SetBranchStatus("z", 1);
    //tree->SetBranchStatus("x_prod", 1);
    //tree->SetBranchStatus("y_prod", 1);
    //tree->SetBranchStatus("z_prod", 1);

    //int partID(0);
    int parentID(0);
    float x(0.0), y(0.0);
    //float z(0.0);
    //float x_prod(0.0), y_prod(0.0), z_prod(0.0);

    //tree->SetBranchAddress("partID", &partID);
    tree->SetBranchAddress("parentID", &parentID);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    //tree->SetBranchAddress("z", &z);
    //tree->SetBranchAddress("x_prod", &x_prod);
    //tree->SetBranchAddress("y_prod", &y_prod);
    //tree->SetBranchAddress("z_prod", &z_prod);

    //Now we loop through the tree
    long long int const entries = tree->GetEntries();
    for (long long int i = 0; i < entries; ++i) {
      tree->GetEntry(i);
      if (parentID == 0) {
        histo_xy_Primaries->Fill(x,y);
        histo_x_Primaries->Fill(x);
        histo_y_Primaries->Fill(y);
      } else {
        histo_xy_Secondaries->Fill(x,y);
        histo_x_Secondaries->Fill(x);
        histo_y_Secondaries->Fill(y);
      }
    } 
    file->Close();
  }

  TCanvas* PDF_Canvas = new TCanvas(); 
  std::string PDF_filename = "output/" + outputfilename + "_all_xyDistributions.pdf";
  std::string PDF_title1 = PDF_filename + "[";
  std::string PDF_title2 = PDF_filename + "]";
  PDF_Canvas->Print(PDF_title1.c_str());

  Draw_Print_Beam(histo_xy_Primaries, canvas, outputfilename + "_xy_Primaries", PDF_filename);
  Draw_Print_Beam(histo_x_Primaries, canvas, outputfilename + "_x_Primaries", PDF_filename);
  Draw_Print_Beam(histo_y_Primaries, canvas, outputfilename + "_y_Primaries", PDF_filename);
  Draw_Print_Beam(histo_xy_Secondaries, canvas, outputfilename + "_xy_Secondaries", PDF_filename);
  Draw_Print_Beam(histo_x_Secondaries, canvas, outputfilename + "_x_Secondaries", PDF_filename);
  Draw_Print_Beam(histo_y_Secondaries, canvas, outputfilename + "_y_Secondaries", PDF_filename);

  PDF_Canvas->Print(PDF_title2.c_str());
  delete PDF_Canvas;
  return 0;
}

void Draw_Print_Beam(TH1* histo, TCanvas* canvas, std::string outputname, std::string PDFfilename){
  canvas->cd();
  histo->Draw("colz");
  canvas->Update();
  TPaveStats *st =
    (TPaveStats*) histo->GetListOfFunctions()->FindObject(
        "stats");
  st->SetX1NDC(0.65); //new x start position
  st->SetX2NDC(0.8); //new x end position
  st->SetY1NDC(0.7); //new x start position
  st->SetY2NDC(0.9); //new x end position
 
  std::string outputname_pdf = "output/" + outputname + ".pdf";
  std::string outputname_cxx = "output/" + outputname + ".cxx";
  canvas->Print(outputname_pdf.c_str());
  canvas->Print(outputname_cxx.c_str());
  canvas->Print(PDFfilename.c_str());
}
