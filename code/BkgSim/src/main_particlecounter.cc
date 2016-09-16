#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveStats.h"


#include <bitset>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

#include "ConfigReaderAnalysis.h"
#include "UsefulFunctions.h"
#include "CellHits_class.h"
#include "GeneralFunctions_SiDBkgSim.h"

#include "Style.h"

using namespace std;

void Print_Origin_histo(TCanvas* canvas, TH2D* const histo, std::string const set_time, std::string const subdetectornames);

int num_hits = 2;

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();
				
				std::vector<std::string> *inputfilenames = new std::vector<std::string>();
				std::vector<std::string> argument_subdetectors;

				int NUMBER_OF_FILES = 0;
				bool NUMBER_OF_FILES_set = false;
				bool inputfile_set = false;
				bool subdetector_set = false;

				for (int i = 1; i < argc; i++) {
								if (argv[i] == std::string("-n")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-s") && argv[i + 1] != std::string("-i")) {
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
																} while (j <= NUMBER_OF_FILES);
																inputfile_set = true;
												} else {
																std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
												}
								}
								if (argv[i] == std::string("-s")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-i")) {
																if (argv[i + 1] != NULL) {
																				int j = 1;
																				do {
																								argument_subdetectors.push_back(argv[i + j]);
																								++j;
																				} while (argv[i + j] != NULL && argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s"));
																				subdetector_set = true;
																} else {
																				std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
																}
												}
								}
				}
				if (!inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set) {
								std::cerr
												<< "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
												<< std::endl;
								exit(1);
				}

				std::vector<Subdetector*> * SubDetectors = new std::vector<Subdetector*>();
				std::string* subdetector_name = new std::string("");
				SetupSubDetectorsVector(SubDetectors, subdetector_name, argument_subdetectors);
				std::string subdetectornames = (*subdetector_name);

				std::vector<CellHits*> AllCounter;

				bool IsCalo = false;
				bool IsEndcap = false;

				for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {
								CellHits * Counter = new CellHits(SubDetectors->at(subdetector_iterator));
								AllCounter.push_back(Counter);

								for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
												TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
												TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

												int HitCellID0(0), HitCellID1(0);
												int LCIO_ID(0);
												float HitCalPosition_x(0.0), HitCalPosition_y(0.0), HitCalPosition_z(0.0);
												double HitTrackerPosition_x(0.0), HitTrackerPosition_y(0.0), HitTrackerPosition_z(0.0);
												float HitMomentum_x(0.0), HitMomentum_y(0.0), HitMomentum_z(0.0);
												double HitMotherMomentum_x(0.0), HitMotherMomentum_y(0.0), HitMotherMomentum_z(0.0);
												double HitMotherVertex_x(0.0), HitMotherVertex_y(0.0), HitMotherVertex_z(0.0);
												double HitParticleVertex_x(0.0), HitParticleVertex_y(0.0), HitParticleVertex_z(0.0);
												float Time(0.0);

												std::vector< double* > HitPosition;
												std::vector< double* > HitParticleVertex;
												std::vector< double* > HitMomentum;

												//Set the branches
												tree->SetBranchStatus("*", 0);
												tree->SetBranchStatus("HitPosition_x", 1);
												tree->SetBranchStatus("HitPosition_y", 1);
												tree->SetBranchStatus("HitPosition_z", 1);

												if (tree->GetName() == std::string("Tree_EcalBarrel") 
																				|| tree->GetName() == std::string("Tree_EcalEndcap")
																				|| tree->GetName() == std::string("Tree_HcalBarrel") 
																				|| tree->GetName() == std::string("Tree_HcalEndcap")
																				|| tree->GetName() == std::string("Tree_MuonBarrel") 
																				|| tree->GetName() == std::string("Tree_MuonEndcap")
																				|| tree->GetName() == std::string("Tree_BeamCal") 
																				|| tree->GetName() == std::string("Tree_LumiCal")) {

																IsCalo = true;
																std::string treename(tree->GetName());
																if (treename.find( "Endcap" ) != std::string::npos) {
																				IsEndcap = true;
																}

																tree->SetBranchStatus("HitCellID0", 1);
																tree->SetBranchAddress("HitCellID0", &HitCellID0);
																tree->SetBranchStatus("HitCellID1", 1);
																tree->SetBranchAddress("HitCellID1", &HitCellID1);
																tree->SetBranchStatus("HitMotherLCIO_id", 1);
																tree->SetBranchAddress("HitMotherLCIO_id", &LCIO_ID);
																tree->SetBranchAddress("HitPosition_x", &HitCalPosition_x);
																tree->SetBranchAddress("HitPosition_y", &HitCalPosition_y);
																tree->SetBranchAddress("HitPosition_z", &HitCalPosition_z);
																tree->SetBranchStatus("HitMotherVertex_x", 1);
																tree->SetBranchStatus("HitMotherVertex_y", 1);
																tree->SetBranchStatus("HitMotherVertex_z", 1);
																tree->SetBranchAddress("HitMotherVertex_x", &HitMotherVertex_x);
																tree->SetBranchAddress("HitMotherVertex_y", &HitMotherVertex_y);
																tree->SetBranchAddress("HitMotherVertex_z", &HitMotherVertex_z);
																tree->SetBranchStatus("HitMotherMomentum_x", 1);
																tree->SetBranchStatus("HitMotherMomentum_y", 1);
																tree->SetBranchStatus("HitMotherMomentum_z", 1);
																tree->SetBranchAddress("HitMotherMomentum_x", &HitMotherMomentum_x);
																tree->SetBranchAddress("HitMotherMomentum_y", &HitMotherMomentum_y);
																tree->SetBranchAddress("HitMotherMomentum_z", &HitMotherMomentum_z);
																tree->SetBranchStatus("HitContrTime", 1);
																tree->SetBranchAddress("HitContrTime", &Time);

																HitPosition = {(double*)&HitCalPosition_x,(double*)&HitCalPosition_y,(double*)&HitCalPosition_z};
																HitParticleVertex = {(double*)&HitMotherVertex_x,(double*)&HitMotherVertex_y,(double*)&HitMotherVertex_z};
																HitMomentum = {(double*)&HitMotherMomentum_x,(double*)&HitMotherMomentum_y,(double*)&HitMotherMomentum_z};
												}
												else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
																				|| tree->GetName() == std::string("Tree_SiVertexEndcap")
																				|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
																				|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
																				|| tree->GetName() == std::string("Tree_SiTrackerForward")) {

																std::string treename = tree->GetName();
																if (treename.find( "Endcap" ) != std::string::npos) {
																				IsEndcap = true;
																}

																tree->SetBranchStatus("HitCellID", 1);
																tree->SetBranchAddress("HitCellID", &HitCellID0);
																//tree->SetBranchStatus("HitCellID0", 1);
																//tree->SetBranchAddress("HitCellID0", &HitCellID0);
																tree->SetBranchStatus("HitParticleLCIO_ID", 1);
																tree->SetBranchAddress("HitParticleLCIO_ID", &LCIO_ID);
																tree->SetBranchAddress("HitPosition_x", &HitTrackerPosition_x);
																tree->SetBranchAddress("HitPosition_y", &HitTrackerPosition_y);
																tree->SetBranchAddress("HitPosition_z", &HitTrackerPosition_z);
																tree->SetBranchStatus("HitParticleVertex_x", 1);
																tree->SetBranchStatus("HitParticleVertex_y", 1);
																tree->SetBranchStatus("HitParticleVertex_z", 1);
																tree->SetBranchAddress("HitParticleVertex_x", &HitParticleVertex_x);
																tree->SetBranchAddress("HitParticleVertex_y", &HitParticleVertex_y);
																tree->SetBranchAddress("HitParticleVertex_z", &HitParticleVertex_z);
																tree->SetBranchStatus("HitMomentum_x", 1);
																tree->SetBranchStatus("HitMomentum_y", 1);
																tree->SetBranchStatus("HitMomentum_z", 1);
																tree->SetBranchAddress("HitMomentum_x", &HitMomentum_x);
																tree->SetBranchAddress("HitMomentum_y", &HitMomentum_y);
																tree->SetBranchAddress("HitMomentum_z", &HitMomentum_z);
																tree->SetBranchStatus("HitTime", 1);
																tree->SetBranchAddress("HitTime", &Time);

																HitPosition = {(double*)&HitTrackerPosition_x,(double*)&HitTrackerPosition_y,(double*)&HitTrackerPosition_z};
																HitParticleVertex = {(double*)&HitParticleVertex_x,(double*)&HitParticleVertex_y,(double*)&HitParticleVertex_z};
																HitMomentum = {(double*)&HitMomentum_x,(double*)&HitMomentum_y,(double*)&HitMomentum_z};
												}

												//Now we loop through the tree
												//Combine the two Cell ID's into a single new Cell ID

												long long int const entries = tree->GetEntries();
												for (long long int i = 0; i < entries; ++i) {
																tree->GetEntry(i);
																//if (IsEndcap && HitPosition.at(2) < 0) continue;

																//Make a combined cell ID
																long long int combined_cell_id = 0;
																if (IsCalo) {
																				combined_cell_id = (long long) HitCellID1 << 32 | HitCellID0;
																}
																else combined_cell_id = (long long) HitCellID0;

																//Use the CellHits class for storing the hit cells and their particlecounts
																Counter->Check_ParticleID(combined_cell_id, LCIO_ID, Time, *HitParticleVertex.at(0), *HitParticleVertex.at(1), *HitParticleVertex.at(2), *HitPosition.at(0), *HitPosition.at(1), *HitPosition.at(2), std::sqrt((*HitMomentum.at(0))* (*HitMomentum.at(0))+(*HitMomentum.at(1))* (*HitMomentum.at(1))), (*HitMomentum.at(2)));
												}
												file->Close();
								}
				}
				
				std::vector<float> range_array;
				range_array = SubDetectors->at(0)->GetROOTHisto_binning3D(); //SOMETHING HARD CODED!!
				float rmax = 1.5*sqrt(pow(range_array[5], 2) + pow(range_array[8], 2));//Make the plot a big bigger than the data
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
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo1 = new TH2D("histo1", title1.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo1);
				std::string const title2 = "10ns <= Time < 20ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo2 = new TH2D("histo2", title2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2);
				std::string const title2_2 = "10ns <= Time < 11ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo2_2 = new TH2D("histo2_2", title2_2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2_2);
				std::string const title2_3 = "11ns <= Time < 13ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo2_3 = new TH2D("histo2_3", title2_3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2_3);
				std::string const title2_4 = "13ns <= Time < 20ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo2_4 = new TH2D("histo2_4", title2_4.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo2_4);
				std::string const title3 = "20ns <= Time < 30ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo3 = new TH2D("histo3", title3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo3);
				std::string const title3_2 = "20ns <= Time < 23ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo3_2 = new TH2D("histo3_2", title3_2.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo3_2);
				std::string const title3_3 = "23ns <= Time < 30ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo3_3 = new TH2D("histo3_3", title3_3.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo3_3);
				std::string const title4 = "30ns <= Time < 50ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo4 = new TH2D("histo4", title4.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo4);
				std::string const title5 = "50ns <= Time < 1000ns, Particle origins for those hitting "
								+ subdetectornames + ";z [mm];r [mm];# of origins";
				TH2D* histo5 = new TH2D("histo5", title5.c_str(), axis_vector[0], axis_vector[1], axis_vector[2], axis_vector[3],
												axis_vector[4], axis_vector[5]);
				histo_vector.push_back(histo5);

				for(int i = 0; i < histo_vector.size(); ++i){
								histo_vector.at(i)->GetZaxis()->SetRangeUser(1,2*std::pow(10,6));
				}

				//Make histogram for storing the information
				std::string const title = "Hits per particle, for subdetector " + subdetectornames +";Number of hits;Number of particles";
				//std::string const title = "Normalized buffer depth for subdetector " + subdetectornames;
				std::vector< TH1D* > histos;
				TH1D* All_Layers_histo = new TH1D("All layers", title.c_str(), 50, 0, 100);
				std::vector< TPaveStats* > stats;

				int tot_no_hits = 0;

				for (size_t subdetector_it = 0; subdetector_it < SubDetectors->size(); ++subdetector_it) {
								int max_num_layers = 0;
								//	if(max_num_layers < SubDetectors->at(subdetector_it)->GetNumberOfLayers()){
								//		max_num_layers = SubDetectors->at(subdetector_it)->GetNumberOfLayers();
								//	}
								max_num_layers = SubDetectors->at(subdetector_it)->GetNumberOfLayers();
								for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
												std::stringstream layername;
												layername << SubDetectors->at(subdetector_it)->GetName() << " Layer " << number_layer;
												histos.emplace_back(new TH1D(layername.str().c_str(), title.c_str(), 50, 0, 100));
								}
				}
				for (size_t allcounters = 0; allcounters < AllCounter.size(); ++allcounters) {
								for (size_t particlecounts = 0; particlecounts < AllCounter.at(allcounters)->Get_ParticleCount().size(); ++particlecounts) {
												if(AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts) > 0){
																if (IsCalo){
																				histos.at(AllCounter.at(allcounters)->Get_Layer().at(particlecounts) + allcounters*SubDetectors->at(allcounters)->GetNumberOfLayers())->Fill(AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts));
																}
																else {//Fill the particle counts into the appropriate histogram: for each subdetector and their layers there is one histogram
																				// -1 in layer number because the vertex and tracker layers start from number 1
																				if(allcounters == 0) histos.at(AllCounter.at(allcounters)->Get_Layer().at(particlecounts) - 1)->Fill(AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts));
																				else histos.at(AllCounter.at(allcounters)->Get_Layer().at(particlecounts) - 1 + SubDetectors->at(allcounters-1)->GetNumberOfLayers())->Fill(AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts));
																}
																All_Layers_histo->Fill(AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts));
																tot_no_hits += AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts);
																//---Only for producing vertex maps again for particles hitting < num_hits times
																/*
																if(AllCounter.at(allcounters)->Get_ParticleCount().at(particlecounts)<num_hits){
																				float time =  AllCounter.at(allcounters)->Get_HitTime().at(particlecounts);
																				std::vector< double > vertex = {
																								AllCounter.at(allcounters)->Get_HitParticleVertex('x').at(particlecounts),
																								AllCounter.at(allcounters)->Get_HitParticleVertex('y').at(particlecounts),
																								AllCounter.at(allcounters)->Get_HitParticleVertex('z').at(particlecounts)
																				};
																				if (time < 10.0) histo1->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																				else if (time >= 10.0 && time < 20.0){
																								histo2->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																								if (time >= 10.0 && time < 11.0) histo2_2->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																								if (time >= 11.0 && time < 13.0) histo2_3->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																								if (time >= 13.0 && time < 20.0) histo2_4->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																				}
																				else if (time >= 20.0 && time < 30.0){
																								histo3->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																								if (time >= 20.0 && time < 23.0) histo3_2->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																								if (time >= 23.0 && time < 30.0) histo3_3->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																				}
																				else if (time >= 30.0 && time < 50.0) histo4->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																				else if (time >= 50.0 && time < 1000.0) histo5->Fill(vertex[2], sqrt(pow(vertex[0], 2) + pow(vertex[1], 2)));
																}
																*/			
												}
								}
				}

				std::cout<< "---------------" <<std::endl;
				std::cout<< "Total number of hits counted for this histogram: " << tot_no_hits <<std::endl;
				std::cout<< "---------------" <<std::endl;

				double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
				for(size_t iterator = 0; iterator < histos.size(); ++iterator){
								histos.at(iterator)->SetMinimum(0.5);
								histos.at(iterator)->SetMaximum(max);
								histos.at(iterator)->Sumw2();
				}
				//NormalizeHistogram(histo, 1.0);
				//Plot the histogram and save it
				TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
				canvas->SetLogy(1);
				float boxsize = 0.0;
				int color = 2; // Very first histogram will be drawn with the color 2, then counted up
				int marker = 20; // Very first histogram will be drawn with the marker style 20, then counted up
				for (size_t number_histo = 0; number_histo< histos.size(); ++number_histo) {
								if(number_histo == 0){
												histos.at(number_histo)->SetLineColor(color);
												histos.at(number_histo)->SetMarkerColor(color);
												histos.at(number_histo)->SetMarkerStyle(marker);
												histos.at(number_histo)->Draw();
												canvas->Update();
												TPaveStats* st =  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats");
												st->SetTextColor(color);
												st->SetX1NDC(0.6); //new x start position
												st->SetX2NDC(0.75); //new x end position
												st->SetY1NDC(0.8); //new y start position
												st->SetY2NDC(0.9); //new y end position
												stats.push_back(st);
												boxsize = stats.at(number_histo)->GetY2NDC() - stats.at(number_histo)->GetY1NDC();
								}
								else if(number_histo > 0){
												color++;
												marker++;
												if(color == 5 || color == 10) color += 1; // 5 would be yellow, 10 would be very light gray 
												histos.at(number_histo)->SetLineColor(color);
												histos.at(number_histo)->SetMarkerColor(color);
												histos.at(number_histo)->SetMarkerStyle(marker);
												histos.at(number_histo)->Draw("SAMES");
												canvas->Update();
												TPaveStats* st =  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats");
												stats.push_back(st);
												stats.at(number_histo)->SetTextColor(color);
												if(number_histo >= 5){
																stats.at(number_histo)->SetX1NDC(0.75); //new x start position
																stats.at(number_histo)->SetX2NDC(0.9); //new x end position
												}
												else {
																stats.at(number_histo)->SetX1NDC(0.6); //new x start position
																stats.at(number_histo)->SetX2NDC(0.75); //new x end position
												}
												if(number_histo == 5) {
																stats.at(number_histo)->SetY1NDC(0.8); //new y end position
																stats.at(number_histo)->SetY2NDC(0.9); //new y end position
												}
												else {
																stats.at(number_histo)->SetY2NDC(stats.at(number_histo-1)->GetY1NDC()); //new y end position
																stats.at(number_histo)->SetY1NDC(stats.at(number_histo)->GetY2NDC()-boxsize); //new y start position
												}
								}
				}
				std::stringstream output;
				output << "output/particlecounts_" << subdetectornames;
				canvas->Print((output.str() + ".pdf").c_str());
				canvas->Print((output.str() + ".cxx").c_str());

				All_Layers_histo->SetMinimum(0.1);
				All_Layers_histo->SetLineColor(2);
				All_Layers_histo->Draw();
				canvas->Update();
				TPaveStats* st =  (TPaveStats*)All_Layers_histo->GetListOfFunctions()->FindObject("stats");
				st->SetX1NDC(0.65); //new x start position
				st->SetX2NDC(0.85); //new x end position
				st->SetY1NDC(0.8); //new y start position
				st->SetY2NDC(0.9); //new y end position

				std::stringstream output2;
				output2 << "output/particlecounts_all_layers_" << subdetectornames;
				canvas->Print((output2.str() + ".pdf").c_str());
				canvas->Print((output2.str() + ".cxx").c_str());
/*
				std::string set_time = "hit";
				Print_Origin_histo(canvas, histo1, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo2, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo2_2, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo2_3, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo2_4, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo3, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo3_2, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo3_3, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo4, set_time, subdetectornames);
				Print_Origin_histo(canvas, histo5, set_time, subdetectornames);
*/
				return 0;
}
void Print_Origin_histo(TCanvas* canvas, TH2D* const histo, std::string const set_time, std::string const subdetectornames){
				canvas->cd();
				canvas->SetLogy(0);
				canvas->SetLogz();
				histo->Draw("colz");
				canvas->Update();
				TPaveStats *st = (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
				st->SetX1NDC(0.6); //new x start position
				st->SetX2NDC(0.85); //new x end position
				st->SetY1NDC(0.6); //new x start position
				st->SetY2NDC(0.9); //new x end position

				std::string histoname(histo->GetName());
				std::stringstream numhits;
				numhits << num_hits;

				canvas->Print(("output/hitmaps_particleorigins_less"+numhits.str()+"hits_"+set_time+"time_"+histoname+"_"+subdetectornames+".pdf").c_str());
				canvas->Print(("output/hitmaps_particleorigins_less"+numhits.str()+"hits_"+set_time+"time_"+histoname+"_"+subdetectornames+".cxx").c_str());
}
