#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include "UsefulFunctions.h"
#include "CellHits_class.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "Time_class.h"

#include "Style.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();

				//The input is a TTree ROOT file(s)

				std::vector<std::string> *inputfilenames = new std::vector<std::string>();
				std::vector<std::string> argument_subdetectors;
				std::string set_time;

				int NUMBER_OF_FILES = 0;
				bool NUMBER_OF_FILES_set = false;
				bool inputfile_set = false;
				bool subdetector_set = false;
				bool time_set = false;

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
								if (argv[i] == std::string("-s")) {
												if (argv[i + 1] != NULL){ 
																int j = 1;
																while (argv[i + j] != NULL 
																								&& argv[i + j] != std::string("-t") 
																								&& argv[i + j] != std::string("-i") 
																								&& argv[i + j] != std::string("-n") 
																								&& argv[i + j] != std::string("-s")){
																				argument_subdetectors.push_back(argv[i + j]);
																				++j;
																}
																subdetector_set = true;
												} else {
																std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
												}
								}
								if (argv[i] == std::string("-t")) {
												if (argv[i + 1] != NULL 
																				&& argv[i + 1] != std::string("-t") 
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-n") 
																				&& argv[i + 1] != std::string("-s")){
																set_time = argv[i + 1];

																time_set = true;
												} else {
																std::cerr << "You didn't give an argument for the time!" << std::endl;
												}
								}
				}
				if (!inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set || !time_set) {
								std::cerr
												<< "You didn't give the name for the subdector, the time, the inputfiles or the amount of files. Please try again!"
												<< std::endl;
								exit(1);
				}

				std::vector<Subdetector*> * SubDetectors = new std::vector<Subdetector*>();
				std::string* subdetector_name = new std::string("");
				SetupSubDetectorsVector(SubDetectors, subdetector_name, argument_subdetectors);
				std::string subdetectornames = (*subdetector_name);

				std::stringstream subdetector_names;

				TCanvas* canvas = new TCanvas();
				std::string histo_title;
				histo_title = "P_T vs hit time of primaries hitting the " + subdetectornames + ";Hit time [ns];P_T [GeV]";
				TH2D* PT_Hittime = new TH2D("PT_HitTime",histo_title.c_str(),60,0,20,60,0,2);

				int total_no_primaries = 0;
				int late_coming_primaries = 0;
				float late_primaries_starttime = 5.0;//ns

				for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {
								subdetector_names << SubDetectors->at(subdetector_iterator)->GetName();

								for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
												TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
												TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

												//Set the branches
												float actualtime = 0.0;
												float creationtime = 0.0;
												double vertex_x = 0.0;
												double vertex_y = 0.0;
												double vertex_z = 0.0;
												float momentum_x = 0.0;
												float momentum_y = 0.0;
												bool createdInSimulation_Status = false;
												
												tree->SetBranchStatus("*", 0);

												if (tree->GetName() == std::string("Tree_EcalBarrel") 
																				|| tree->GetName() == std::string("Tree_EcalEndcap")
																				|| tree->GetName() == std::string("Tree_HcalBarrel") 
																				|| tree->GetName() == std::string("Tree_HcalEndcap")
																				|| tree->GetName() == std::string("Tree_MuonBarrel") 
																				|| tree->GetName() == std::string("Tree_MuonEndcap")
																				|| tree->GetName() == std::string("Tree_BeamCal") 
																				|| tree->GetName() == std::string("Tree_LumiCal")) {
																tree->SetBranchStatus("HitContrTime", 1);
																tree->SetBranchStatus("HitMotherVertex_x", 1);
																tree->SetBranchStatus("HitMotherVertex_y", 1);
																tree->SetBranchStatus("HitMotherVertex_z", 1);
																tree->SetBranchAddress("HitContrTime", &actualtime);
																tree->SetBranchAddress("HitMotherVertex_x", &vertex_x);
																tree->SetBranchAddress("HitMotherVertex_y", &vertex_y);
																tree->SetBranchAddress("HitMotherVertex_z", &vertex_z);
												}
												else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
																				|| tree->GetName() == std::string("Tree_SiVertexEndcap")
																				|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
																				|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
																				|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
																tree->SetBranchStatus("HitTime", 1);
																tree->SetBranchStatus("HitParticleCreationTime", 1);
																tree->SetBranchStatus("HitParticleVertex_x", 1);
																tree->SetBranchStatus("HitParticleVertex_y", 1);
																tree->SetBranchStatus("HitParticleVertex_z", 1);
																tree->SetBranchStatus("HitMomentum_x", 1);
																tree->SetBranchStatus("HitMomentum_y", 1);
																tree->SetBranchStatus("HitParticle_CreatedInSimulation_Status", 1);
																tree->SetBranchAddress("HitTime", &actualtime);
																tree->SetBranchAddress("HitParticleCreationTime", &creationtime);
																tree->SetBranchAddress("HitParticleVertex_x", &vertex_x);
																tree->SetBranchAddress("HitParticleVertex_y", &vertex_y);
																tree->SetBranchAddress("HitParticleVertex_z", &vertex_z);
																tree->SetBranchAddress("HitMomentum_x", &momentum_x);
																tree->SetBranchAddress("HitMomentum_y", &momentum_y);
																tree->SetBranchAddress("HitParticle_CreatedInSimulation_Status", &createdInSimulation_Status);
												} else {
																std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
																std::terminate();
												}

												float* time = nullptr;
												if (set_time == "creation") time = &creationtime;
												else if (set_time == "hit") time = &actualtime;
												else{
																std::cerr << "The time variable wasn't set properly! You have the choice between 'creation' and 'hit' time." << std::endl;
																exit(1);
												}

												std::array<double, 3> vertex;

												//Now we loop through the tree

												long long int const entries = tree->GetEntries();
												for (long long int i = 0; i < entries; ++i) {
																tree->GetEntry(i);
																vertex = { vertex_x, vertex_y, vertex_z };

																//Only interested in primary particles from IP:
																if ( (vertex[0]<-5 || vertex[0]>5) && 
																								(vertex[1]<-5 || vertex[1]>5) && 
																								(vertex[2]<-5 || vertex[2]>5) &&
																								createdInSimulation_Status == 1){
																				continue;
																}

																PT_Hittime->Fill(*time,std::sqrt(momentum_x*momentum_x+momentum_y*momentum_y));
																
																//Count how many primaries there are in total
																total_no_primaries++;	
																//Count how many primaries arrive late at the subdetector
																if (*time > late_primaries_starttime){
																				late_coming_primaries++;
																}

												}
												file->Close();
								}
				}
				std::cout << "Total number of primaries hitting the subdetector " << subdetector_names.str() << ": " << total_no_primaries << std::endl;
				std::cout << "Number of primaries hitting the subdetector " << subdetector_names.str() << " later than " << late_primaries_starttime << "ns: " << late_coming_primaries << std::endl;

				canvas->SetLogz();
				PT_Hittime->Draw("colz");

				canvas->Update();
				TPaveStats *st = (TPaveStats*)PT_Hittime->GetListOfFunctions()->FindObject("stats");
				st->SetX1NDC(0.6); //new x start position
				st->SetX2NDC(0.85); //new x end position
				st->SetY1NDC(0.6); //new x start position
				st->SetY2NDC(0.9); //new x end position

				std::stringstream output;
				output << "output/PT_hittime_primaries_" << subdetector_names.str();
				canvas->Print((output.str() + ".pdf").c_str());
				canvas->Print((output.str() + ".cxx").c_str());

				return 0;
}

