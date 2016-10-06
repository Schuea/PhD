#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TLine.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include "ConfigReaderAnalysis.h"
#include "UsefulFunctions.h"
#include "CellHits_class.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "Helix_class.h"

#include "Style.h"

using namespace std;

void Print_Origin_histo(TCanvas* canvas, TH2D* const histo, std::string const set_time, std::string const subdetectornames);

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();

				std::vector<std::string> *inputfilenames = new std::vector<std::string>();

				int NUMBER_OF_FILES = 0;
				bool NUMBER_OF_FILES_set = false;
				bool inputfile_set = false;

				for (int i = 1; i < argc; i++) {
								if (argv[i] == std::string("-n")) {
												if (argv[i + 1] != NULL 
																				&& argv[i + 1] != std::string("-s") 
																				&& argv[i + 1] != std::string("-n") 
																				&& argv[i + 1] != std::string("-t") 
																				&& argv[i + 1] != std::string("-i")) {
																NUMBER_OF_FILES = std::stoi(argv[i + 1]);
																std::cout << "Number of input files = " << NUMBER_OF_FILES << std::endl;
																NUMBER_OF_FILES_set = true;
												} else {
																std::cerr << "You didn't give an argument for the number of files!" << std::endl;
												}
								}
				}
				for (int i = 1; i < argc; i++) {
								if (argv[i] == std::string("-i")){
												if (argv[i + 1] != NULL) {
																int j = 1;
																do {
																				if (argv[i + j] != std::string("-n") 
																												&& argv[i + j] != std::string("-i") 
																												&& argv[i + j] != std::string("-n") 
																												&& argv[i + j] != std::string("-s") 
																												&& argv[i + j] != std::string("-t")) {
																								inputfilenames->push_back(argv[i + j]);
																								++j;
																				} else {
																								break;
																				}
																} while (j <= NUMBER_OF_FILES);
																inputfile_set = true;
												} else {
																std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
												}
								}
				}
				if (!inputfile_set || !NUMBER_OF_FILES_set) {
								std::cerr
												<< "You didn't give the name for the inputfiles or the amount of files. Please try again!"
												<< std::endl;
								exit(1);
				}

				//Make histogram for storing the information
				float zmin = 0.0;
				float zmax = 300.0;
				int zbin = 10000;

				std::string const title_x = "Pairs spiraling in the magnetic field;z [mm];x [mm];# of particles per (0.85mm x 0.25mm)";
				std::string const title_y = "Pairs spiraling in the magnetic field;z [mm];y [mm];# of particles per (0.85mm x 0.25mm)";
				TH2D* histo_x = new TH2D("Helix_tracks_xz", title_x.c_str(), zbin,zmin,zmax, 100,-27.5,27.5);
				TH2D* histo_y = new TH2D("Helix_tracks_yz", title_y.c_str(), zbin,zmin,zmax, 100,-27.5,27.5);

				float BField = 5.0;
				Helix helix(BField);

				for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
								TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
								TTree *tree = Get_TTree(file, "MCP");

								//Set the branches
								double vertex_x = 0.0;
								double vertex_y = 0.0;
								double vertex_z = 0.0;
								double momentum_x = 0.0;
								double momentum_y = 0.0;
								double momentum_z = 0.0;
								float charge = -99.0;
								bool CreatedInSimulation_Status = false;

								tree->SetBranchStatus("*", 0);
								tree->SetBranchStatus("Vertexx", 1);
								tree->SetBranchStatus("Vertexy", 1);
								tree->SetBranchStatus("Vertexz", 1);
								tree->SetBranchAddress("Vertexx", &vertex_x);
								tree->SetBranchAddress("Vertexy", &vertex_y);
								tree->SetBranchAddress("Vertexz", &vertex_z);
								tree->SetBranchStatus("Momentumx", 1);
								tree->SetBranchStatus("Momentumy", 1);
								tree->SetBranchStatus("Momentumz", 1);
								tree->SetBranchAddress("Momentumx", &momentum_x);
								tree->SetBranchAddress("Momentumy", &momentum_y);
								tree->SetBranchAddress("Momentumz", &momentum_z);
								tree->SetBranchStatus("Charge", 1);
								tree->SetBranchAddress("Charge", &charge);
								tree->SetBranchStatus("CreatedInSimulation_Status", 1);
								tree->SetBranchAddress("CreatedInSimulation_Status", &CreatedInSimulation_Status);

								std::vector< double > vertex;
								std::vector< double > momentum;
								double z = zmin;

								long long int const entries = tree->GetEntries();
								for (long long int i = 0; i < entries; ++i) {
												tree->GetEntry(i);
												if (CreatedInSimulation_Status == 1) continue;
												vertex = { vertex_x, vertex_y, vertex_z };
												if (abs(vertex.at(0)) > 0.1 || abs(vertex.at(1)) > 0.1 ) continue;
												momentum = { momentum_x, momentum_y, momentum_z };
												//vertex = { 0,0,0 };
												//momentum = { 0.1,0.2,3 };
												for (int step = 1; step <= zbin; ++step){
																double new_x = helix.Get_position(momentum, charge, vertex, z).at(0)*1000.0; // to convert from m to mm
																double new_y = helix.Get_position(momentum, charge, vertex, z).at(1)*1000.0; // to convert from m to mm
																//double new_x = helix.Get_position(momentum, 1, vertex, z).at(0)*1000.0;
																//double new_y = helix.Get_position(momentum, 1, vertex, z).at(1)*1000.0;
																//std::cout << "z = " << z << std::endl;
																//std::cout << "new_x = " << new_x << std::endl;
																//std::cout << "new_y = " << new_y << std::endl;
																histo_x->Fill(z, new_x);
																histo_y->Fill(z, new_y);
																z = step*(zmax-zmin)/zbin;
												}
								}
								file->Close();
				}
				TLine *line = new TLine(0,12,62.5,12);
				TLine *nline = new TLine(0,-12,62.5,-12);
				TLine *line2 = new TLine(62.5,12,200,20);
				TLine *nline2 = new TLine(62.5,-12,200,-20);
				TLine *line3 = new TLine(200,20,300,27.5);
				TLine *nline3 = new TLine(200,-20,300,-27.5);
				line->SetLineColor(2);
				nline->SetLineColor(2);
				line2->SetLineColor(2);
				nline2->SetLineColor(2);
				line3->SetLineColor(2);
				nline3->SetLineColor(2);
				
				gStyle->SetOptStat(0);
				//gStyle->SetOptStat(111111);

				//Plot the histogram and save it

				TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
				canvas->cd();

				canvas->SetLogz();
				histo_x->Draw("colz");
				//canvas->Update();
				//TPaveStats *st = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
				//st->SetX1NDC(0.6); //new x start position
				//st->SetX2NDC(0.85); //new x end position
				//st->SetY1NDC(0.6); //new x start position
				//st->SetY2NDC(0.9); //new x end position

				line->Draw();
				nline->Draw();
				line2->Draw();
				nline2->Draw();
				line3->Draw();
				nline3->Draw();

				std::string histoname_x(histo_x->GetName());

				canvas->Print(("output/"+histoname_x+".pdf").c_str());
				canvas->Print(("output/"+histoname_x+".cxx").c_str());

				canvas->SetLogz();
				histo_y->Draw("colz");
				line->Draw();
				nline->Draw();
				line2->Draw();
				nline2->Draw();
				line3->Draw();
				nline3->Draw();

				std::string histoname_y(histo_y->GetName());

				canvas->Print(("output/"+histoname_y+".pdf").c_str());
				canvas->Print(("output/"+histoname_y+".cxx").c_str());
				
				return 0;
}
