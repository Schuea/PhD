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
#include "CellHits_class_new.h"
#include "GeneralFunctions_SiDBkgSim_new.h"
#include "Time_class.h"

#include "Style.h"

using namespace std;

void Print_Origin_histo(TCanvas* canvas, TH2D* const histo, std::string const set_time, std::string const subdetectornames);

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();

				//ConfigReaderAnalysis config(argv[1]);
				//config.setUp();
				//cout << config.getConfigName() << endl;
				//cout << config.getDoEventLooper() << endl;
				//cout << config.getMaxEvents() << endl;
				//exit(1);
				//Three occupancy plots
				//Two that are 1D histograms, y-axis is average occupancy and x-axis is radius/phi
				//The difficult plot is buffer depth plot: y-axis is probability of a specific cell occupancy occuring and x-axis is occupancy

				//The input is a TTree ROOT file(s)
				//The output is .pdf and .C files

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

        double weight = 1.0;
        //weight = (double)(1.0/double(NUMBER_OF_FILES));
				std::string subdetectornames = argument_subdetectors.at(0);
				//std::vector<float> range_array;

				//range_array = SubDetectors->at(0)->GetROOTHisto_binning3D(); //SOMETHING HARD CODED!!
				float rmax = 225; 
				//float rmax = 1.5*sqrt(pow(range_array[5], 2) + pow(range_array[8], 2));//Make the plot a big bigger than the data
				float rmin = 0.;
				float rrange = rmax - rmin;
				float zmax = 3500;
				//float zmax = range_array[2];
				float zmin = -zmax;
				float zrange = rmax - zmin;
				std::array<float, 6> axis_vector = { float(zrange / 20.0), zmin, zmax, float(rrange / 7.0), rmin, rmax};

				//Make histogram for storing the information
				std::vector< TH2D* > histo_vector;
				std::string const title1 = "Time < 10ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo1 = new TH2D("0ns-10ns", title1.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo1);
				std::string const title2 = "10ns <= Time < 20ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo2 = new TH2D("10ns-20ns", title2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2);
				std::string const title2_2 = "10ns <= Time < 11ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo2_2 = new TH2D("10ns-11ns", title2_2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2_2);
				std::string const title2_3 = "11ns <= Time < 13ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo2_3 = new TH2D("11ns-13ns", title2_3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2_3);
				std::string const title2_4 = "13ns <= Time < 20ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo2_4 = new TH2D("13ns-20ns", title2_4.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2_4);
				std::string const title3 = "20ns <= Time < 30ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo3 = new TH2D("20ns-30ns", title3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo3);
				std::string const title3_2 = "20ns <= Time < 23ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo3_2 = new TH2D("20ns-23ns", title3_2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo3_2);
				std::string const title3_3 = "23ns <= Time < 30ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo3_3 = new TH2D("23ns-30ns", title3_3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo3_3);
				std::string const title4 = "30ns <= Time < 50ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo4 = new TH2D("30ns-50ns", title4.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo4);
				std::string const title5 = "50ns <= Time < 1000ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];Number of origins";
				TH2D* histo5 = new TH2D("50ns-1000ns", title5.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo5);

				for(int i = 0; i < histo_vector.size(); ++i){
								histo_vector.at(i)->GetZaxis()->SetRangeUser(1,2*std::pow(10,6));
				}

				std::stringstream subdetector_names;

				for (size_t subdetector_iterator = 0; subdetector_iterator < argument_subdetectors.size(); ++subdetector_iterator) {

                Subdetector det( argument_subdetectors.at(subdetector_iterator) );
				        subdetector_names<< det.getName();
                std::cout << det.getName() << std::endl;
								for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
												TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
												TTree *tree = Get_TTree(file, det.getName());

												//Set the branches
												float actualtime = 0.0;
												float creationtime = 0.0;
												double vertex_x = 0.0;
												double vertex_y = 0.0;
												double vertex_z = 0.0;
												float momentum_x = 0.0;
												float momentum_y = 0.0;
												float momentum_z = 0.0;

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
																tree->SetBranchStatus("HitMomentum_z", 1);
																tree->SetBranchAddress("HitTime", &actualtime);
																tree->SetBranchAddress("HitParticleCreationTime", &creationtime);
																tree->SetBranchAddress("HitParticleVertex_x", &vertex_x);
																tree->SetBranchAddress("HitParticleVertex_y", &vertex_y);
																tree->SetBranchAddress("HitParticleVertex_z", &vertex_z);
																tree->SetBranchAddress("HitMomentum_x", &momentum_x);
																tree->SetBranchAddress("HitMomentum_y", &momentum_y);
																tree->SetBranchAddress("HitMomentum_z", &momentum_z);
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
												std::pair<int, int> Number_train_bunch = Set_train_bunch_number(file_iterator);

												//Now we loop through the tree
												//Combine the two Cell ID's into a single new Cell ID
												//See how often the new Cell ID occurs in total, this is the occupancy

												long long int const entries = tree->GetEntries();
												for (long long int i = 0; i < entries; ++i) {
																tree->GetEntry(i);
																vertex = { vertex_x, vertex_y, vertex_z };
																if (*time < 10.0) histo1->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																else if (*time >= 10.0 && *time < 20.0){
																				histo2->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																				if (*time >= 10.0 && *time < 11.0) histo2_2->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																				if (*time >= 11.0 && *time < 13.0) histo2_3->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																				if (*time >= 13.0 && *time < 20.0) histo2_4->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																}
																else if (*time >= 20.0 && *time < 30.0){
																				histo3->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																				if (*time >= 20.0 && *time < 23.0) histo3_2->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																				if (*time >= 23.0 && *time < 30.0) histo3_3->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																}
																else if (*time >= 30.0 && *time < 50.0)   histo4->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
																else if (*time >= 50.0 && *time < 1000.0) histo5->Fill(vertex[2], sqrt(pow(vertex[0], 2.0) + pow(vertex[1], 2.0)),weight);
												}
												file->Close();
								}
				}
				gStyle->SetOptStat(111111);

				//Plot the histogram and save it

				TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);

				Print_Origin_histo(canvas, histo1, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo2, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo2_2, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo2_3, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo2_4, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo3, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo3_2, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo3_3, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo4, set_time, subdetector_names.str());
				Print_Origin_histo(canvas, histo5, set_time, subdetector_names.str());


				return 0;
}

void Print_Origin_histo(TCanvas* canvas, TH2D* const histo, std::string const set_time, std::string const subdetectornames){
				canvas->cd();
				canvas->SetLogz();
				histo->Draw("colz");
				canvas->Update();
				TPaveStats *st = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
				st->SetX1NDC(0.6); //new x start position
				st->SetX2NDC(0.85); //new x end position
				st->SetY1NDC(0.6); //new x start position
				st->SetY2NDC(0.9); //new x end position

				std::string histoname(histo->GetName());

				canvas->Print(("output/hitmaps_particleorigins_"+set_time+"time_"+histoname+"_"+subdetectornames+"_oldLStar.pdf").c_str());
				canvas->Print(("output/hitmaps_particleorigins_"+set_time+"time_"+histoname+"_"+subdetectornames+"_oldLStar.cxx").c_str());
				//canvas->Print(("output/hitmaps_particleorigins_"+set_time+"time_"+histoname+"_"+subdetectornames+"_newLStar.pdf").c_str());
				//canvas->Print(("output/hitmaps_particleorigins_"+set_time+"time_"+histoname+"_"+subdetectornames+"_newLStar.cxx").c_str());
}
