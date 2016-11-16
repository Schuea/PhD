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

#include "ConfigReaderAnalysis.h"
#include "UsefulFunctions.h"
#include "CellHits_class.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "Time_class.h"

#include "Style.h"

using namespace std;

int main(int const argc, char const * const * const argv) {
				UsePhDStyle();
				TH1::SetDefaultSumw2();

				std::vector<std::string> *pair_inputfilenames = new std::vector<std::string>();
				std::string ewpw_inputfilename;
				std::string ewpb_inputfilename;
				std::string ebpw_inputfilename;
				std::string ebpb_inputfilename;
				std::vector<std::string> argument_subdetectors;

				int NUMBER_OF_pairFILES = 0;
				bool NUMBER_OF_pairFILES_set = false;
				float weight_bunches = 0;
				bool weight_set = false;
				bool pair_inputfile_set = false;
				bool ewpw_inputfile_set = false;
				bool ewpb_inputfile_set = false;
				bool ebpw_inputfile_set = false;
				bool ebpb_inputfile_set = false;
				bool subdetector_set = false;

				for (int i = 1; i < argc; i++) {
								if (argv[i] == std::string("-n")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-s") 
																				&& argv[i + 1] != std::string("-w")
																				&& argv[i + 1] != std::string("-i")
																				&& argv[i + 1] != std::string("-ewpw") 
																				&& argv[i + 1] != std::string("-ewpb") 
																				&& argv[i + 1] != std::string("-ebpw") 
																				&& argv[i + 1] != std::string("-ebpb")) {
																NUMBER_OF_pairFILES = std::stoi(argv[i + 1]);
																std::cout << "Number of pair input files = " << NUMBER_OF_pairFILES << std::endl;
																NUMBER_OF_pairFILES_set = true;
												} else {
																std::cerr << "You didn't give an argument for the number of pair files!" << std::endl;
												}
								}
				}
				for (int i = 1; i < argc; i++) {
								if (argv[i] == std::string("-i") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-w") 
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL) {
																int j = 1;
																do {
																				if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
																												&& argv[i + j] != std::string("-w") 
																												&& argv[i + j] != std::string("-ewpw") 
																												&& argv[i + j] != std::string("-ewpb") 
																												&& argv[i + j] != std::string("-ebpw") 
																												&& argv[i + j] != std::string("-ebpb")) {
																								pair_inputfilenames->push_back(argv[i + j]);
																								++j;
																				} else {
																								break;
																				}
																} while (j <= NUMBER_OF_pairFILES);
																pair_inputfile_set = true;
												} else {
																std::cerr << "You didn't give an argument for the pair inputfile(s)!" << std::endl;
												}
								}
								if (argv[i] == std::string("-w") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-ewpw") 
																				&& argv[i + 1] != std::string("-ewpb") 
																				&& argv[i + 1] != std::string("-ebpw") 
																				&& argv[i + 1] != std::string("-ebpb")) {
																std::string str = argv[i + 1];
																weight_bunches = std::atof(str.c_str());
																std::cout << "weight_bunches = " << weight_bunches << std::endl; 
												} else {
																std::cerr << "You didn't give an argument for the weight!" << std::endl;
																exit(1);
												}
												weight_set = true;
								}

								if (argv[i] == std::string("-ewpw") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-w") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-w") 
																				&& argv[i + 1] != std::string("-ewpb") 
																				&& argv[i + 1] != std::string("-ebpw") 
																				&& argv[i + 1] != std::string("-ebpb")) {
																ewpw_inputfilename = argv[i + 1];
												} else {
																std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
																exit(1);
												}
												ewpw_inputfile_set = true;
								}
								if (argv[i] == std::string("-ewpb") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-w") 
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-w") 
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-ewpw") 
																				&& argv[i + 1] != std::string("-ebpw") 
																				&& argv[i + 1] != std::string("-ebpb")) {
																ewpb_inputfilename = argv[i + 1];
												} else {
																std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
																exit(1);
												}
												ewpb_inputfile_set = true;
								}	
								if (argv[i] == std::string("-ebpw") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-w") 
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-w") 
																				&& argv[i + 1] != std::string("-ewpb") 
																				&& argv[i + 1] != std::string("-ewpw") 
																				&& argv[i + 1] != std::string("-ebpb")) {
																ebpw_inputfilename = argv[i + 1];
												} else {
																std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
																exit(1);
												}
												ebpw_inputfile_set = true;
								}
								if (argv[i] == std::string("-ebpb") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-w") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ewpw")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-w") 
																				&& argv[i + 1] != std::string("-ewpb") 
																				&& argv[i + 1] != std::string("-ebpw") 
																				&& argv[i + 1] != std::string("-ewpw")) {
																ebpb_inputfilename = argv[i + 1];
												} else {
																std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
																exit(1);
												}
												ebpb_inputfile_set = true;
								}
								if (argv[i] == std::string("-s")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") 
																				&& argv[i + 1] != std::string("-i") 
																				&& argv[i + 1] != std::string("-w") 
																				&& argv[i + 1] != std::string("-ewpw") 
																				&& argv[i + 1] != std::string("-ewpb") 
																				&& argv[i + 1] != std::string("-ebpw") 
																				&& argv[i + 1] != std::string("-ebpb")) {
																if (argv[i + 1] != NULL) {
																				int j = 1;
																				do {
																								argument_subdetectors.push_back(argv[i + j]);
																								++j;
																				} while (argv[i + j] != NULL && argv[i + j] != std::string("-n") 
																												&& argv[i + j] != std::string("-i")
																												&& argv[i + j] != std::string("-w") 
																												&& argv[i + j] != std::string("-ewpw") 
																												&& argv[i + j] != std::string("-ewpb") 
																												&& argv[i + j] != std::string("-ebpw") 
																												&& argv[i + j] != std::string("-ebpb"));
																				subdetector_set = true;
																} else {
																				std::cerr << "You didn't give an argument for the subdetector!" << std::endl;
																}
												}
								}
				}
				if (!pair_inputfile_set || !ewpw_inputfile_set || !ewpb_inputfile_set || !ebpw_inputfile_set || !ebpb_inputfile_set || !subdetector_set || !NUMBER_OF_pairFILES_set || !weight_set) {
								std::cerr
												<< "You didn't give the name for the subdector, the inputfiles, the amount of pair inputfiles, or the bunch weight. Please try again!"
												<< std::endl;
								exit(1);
				}

	std::vector<std::string> *inputfiles = new std::vector<std::string>();
	inputfiles->push_back(ewpw_inputfilename);
	inputfiles->push_back(ewpb_inputfilename);
	inputfiles->push_back(ebpw_inputfilename);
	inputfiles->push_back(ebpb_inputfilename);
	inputfiles->insert(inputfiles->end(), pair_inputfilenames->begin(), pair_inputfilenames->end());
	

	std::vector<Subdetector*> * SubDetectors = new std::vector<Subdetector*>();
	std::string* subdetector_name = new std::string("");
	SetupSubDetectorsVector(SubDetectors, subdetector_name, argument_subdetectors);
	std::string subdetectornames = (*subdetector_name);

	//Make histogram for storing the information
	std::string const histo_title = "P_T of particles from #gamma#gamma #rightarrow hadrons events and pair background, " + subdetectornames + ";P_T [GeV]; #hits";
	TH1D* Hits_P_T_gg = new TH1D("P_T_gammagamma_hadrons",histo_title.c_str(),30,0,5);
	TH1D* Hits_P_T_pairs = new TH1D("P_T_pairs",histo_title.c_str(),30,0,5);
	
	//std::string const histo_title_x = "P_x of particles from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";P_x [GeV]; #events";
	//std::string const histo_title_y = "P_y of particles from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";P_y [GeV]; #events";
	//std::string const histo_title_z = "P_z of particles from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";P_z [GeV]; #events";
	std::string const histo_title_x = "P_x of all except e+ e- (leaving detector) from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";P_x [GeV]; #events";
	std::string const histo_title_y = "P_y of all except e+ e- (leaving detector) from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";P_y [GeV]; #events";
	std::string const histo_title_z = "P_z of all except e+ e- (leaving detector) from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";P_z [GeV]; #events";
	std::string const histo_title_E = "E_sum of e+ e- (leaving detector) from #gamma#gamma #rightarrow hadrons events , " + subdetectornames + ";E_sum [GeV]; #events";
	TH1D* Hits_P_x_gg = new TH1D("P_x_gammagamma_hadrons",histo_title_x.c_str(),40,-3,3);
	TH1D* Hits_P_y_gg = new TH1D("P_y_gammagamma_hadrons",histo_title_y.c_str(),40,-3,3);
	TH1D* Hits_P_z_gg = new TH1D("P_z_gammagamma_hadrons",histo_title_z.c_str(),180,-270,270);
	TH1D* Hits_Esum_gg = new TH1D("E_sum_gammagamma_hadrons",histo_title_E.c_str(),60,0,40);

	//std::map< int, std::pair< std::vector<int>, std::vector<double> > > energy_map; // int: event_id, int: particle_id doubles particle energy 
	std::vector< int > eventid_vec;
	std::vector< int > particleid_vec;
	std::vector< double > energy_vec;

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {

		for (size_t file_iterator = 0; file_iterator < inputfiles->size(); ++file_iterator) {
		
			TFile *file = TFile::Open(inputfiles->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			double mom_x = 0.0;
			double mom_y = 0.0;
			double mom_z = 0.0;
			float charge = -1.0;
			double energy = -999.0;
			bool CreatedInSimulation_Status = 0;
			bool DecayedInTracker_Status = 0;
			bool hasLeftDetector_Status = 0;
			int parents = -1;
			int particle_pdg = -999;
			int particle_id = -999;
			int event_id = -999;
			tree->SetBranchStatus("*", 0);
			
			if (tree->GetName() == std::string("Tree_MCP")){
							tree->SetBranchStatus("NumberOfParents", kTRUE);
							tree->SetBranchAddress("NumberOfParents", &parents);
							tree->SetBranchStatus("CreatedInSimulation_Status", kTRUE);
							tree->SetBranchAddress("CreatedInSimulation_Status", &CreatedInSimulation_Status);
							tree->SetBranchStatus("DecayedInTracker_Status", kTRUE);
							tree->SetBranchAddress("DecayedInTracker_Status", &DecayedInTracker_Status);
							tree->SetBranchStatus("hasLeftDetector_Status", kTRUE);
							tree->SetBranchAddress("hasLeftDetector_Status", &hasLeftDetector_Status);
							if(file_iterator<4){
											tree->SetBranchStatus("Particle_PDG", kTRUE);
											tree->SetBranchAddress("Particle_PDG", &particle_pdg);
											tree->SetBranchStatus("Particle_ID", kTRUE);
											tree->SetBranchAddress("Particle_ID", &particle_id);
							}
							else{
											tree->SetBranchStatus("Particle_PDG", kTRUE);
											tree->SetBranchAddress("Particle_PDG", &particle_pdg);
											tree->SetBranchStatus("Particle_ID", kTRUE);
											tree->SetBranchAddress("Particle_ID", &particle_id);
							}
							tree->SetBranchStatus("Event_ID", kTRUE);
							tree->SetBranchAddress("Event_ID", &event_id);
							tree->SetBranchStatus("Momentumx", 1);
							tree->SetBranchAddress("Momentumx", &mom_x);
							tree->SetBranchStatus("Momentumy", 1);
							tree->SetBranchAddress("Momentumy", &mom_y);
							tree->SetBranchStatus("Momentumz", 1);
							tree->SetBranchAddress("Momentumz", &mom_z);
							tree->SetBranchStatus("Charge", 1);
							tree->SetBranchAddress("Charge", &charge);
							tree->SetBranchStatus("Energy", 1);
							tree->SetBranchAddress("Energy", &energy);
			}
			else if (tree->GetName() == std::string("Tree_EcalBarrel") 
											|| tree->GetName() == std::string("Tree_EcalEndcap")
											|| tree->GetName() == std::string("Tree_HcalBarrel") 
											|| tree->GetName() == std::string("Tree_HcalEndcap")
											|| tree->GetName() == std::string("Tree_MuonBarrel") 
											|| tree->GetName() == std::string("Tree_MuonEndcap")
											|| tree->GetName() == std::string("Tree_BeamCal") 
											|| tree->GetName() == std::string("Tree_LumiCal")) {
							tree->SetBranchStatus("HitMother_CreatedInSimulation_Status", kTRUE);
							tree->SetBranchAddress("HitMother_CreatedInSimulation_Status", &CreatedInSimulation_Status);
							tree->SetBranchStatus("HitMother_DecayedInTracker_Status", kTRUE);
							tree->SetBranchAddress("HitMother_DecayedInTracker_Status", &DecayedInTracker_Status);
							tree->SetBranchStatus("HitMother_hasLeftDetector_Status", kTRUE);
							tree->SetBranchAddress("HitMother_hasLeftDetector_Status", &hasLeftDetector_Status);
							tree->SetBranchStatus("event_id", kTRUE);
							tree->SetBranchAddress("event_id", &event_id);
							tree->SetBranchStatus("HitMotherMomentum_x", 1);
							tree->SetBranchAddress("HitMotherMomentum_x", &mom_x);
							tree->SetBranchStatus("HitMotherMomentum_y", 1);
							tree->SetBranchAddress("HitMotherMomentum_y", &mom_y);
							tree->SetBranchStatus("HitMotherMomentum_z", 1);
							tree->SetBranchAddress("HitMotherMomentum_z", &mom_z);
							tree->SetBranchStatus("HitMotherCharge", 1);
							tree->SetBranchAddress("HitMotherCharge", &charge);
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
							tree->SetBranchStatus("HitParticle_CreatedInSimulation_Status", kTRUE);
							tree->SetBranchAddress("HitParticle_CreatedInSimulation_Status", &CreatedInSimulation_Status);
							tree->SetBranchStatus("HitParticle_DecayedInTracker_Status", kTRUE);
							tree->SetBranchAddress("HitParticle_DecayedInTracker_Status", &DecayedInTracker_Status);
							tree->SetBranchStatus("HitParticle_hasLeftDetector_Status", kTRUE);
							tree->SetBranchAddress("HitParticle_hasLeftDetector_Status", &hasLeftDetector_Status);
							tree->SetBranchStatus("event_id", kTRUE);
							tree->SetBranchAddress("event_id", &event_id);
							tree->SetBranchStatus("HitParticleMomentum_x", 1);
							tree->SetBranchAddress("HitParticleMomentum_x", &mom_x);
							tree->SetBranchStatus("HitParticleMomentum_y", 1);
							tree->SetBranchAddress("HitParticleMomentum_y", &mom_y);
							tree->SetBranchStatus("HitParticleMomentum_z", 1);
							tree->SetBranchAddress("HitParticleMomentum_z", &mom_z);
							tree->SetBranchStatus("HitParticleCharge", 1);
							tree->SetBranchAddress("HitParticleCharge", &charge);
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}

			//Now we loop through the tree
			long long int const entries = tree->GetEntries();
			double weight(1);
			if( file_iterator == 0) weight = (291.0/10000.0)*(weight_bunches/1312.0);
			if( file_iterator == 1) weight = (501.0/10000.0)*(weight_bunches/1312.0);
			if( file_iterator == 2) weight = (502.0/10000.0)*(weight_bunches/1312.0);
			if( file_iterator == 3) weight = (829.0/10000.0)*(weight_bunches/1312.0);
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				//if(sqrt(mom_x*mom_x+mom_y*mom_y)<1.6 || sqrt(mom_x*mom_x+mom_y*mom_y)>1.8) continue;
				if ( particle_pdg != 11 && particle_pdg != -11 && hasLeftDetector_Status != 1) continue;
				//if ( (particle_pdg == 11 || particle_pdg == -11) && hasLeftDetector_Status == 1) continue;
				//if(charge == 0) continue;
				//if(DecayedInTracker_Status == 1 && hasLeftDetector_Status == 0 && CreatedInSimulation_Status == 1) continue;
				if( file_iterator < 4){
								Hits_P_T_gg->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight);
								Hits_P_x_gg->Fill(mom_x);
								Hits_P_y_gg->Fill(mom_y);
								Hits_P_z_gg->Fill(mom_z);
								if (tree->GetName() == std::string("Tree_MCP")){
									if ( (particle_pdg == 11 || particle_pdg==-11) && hasLeftDetector_Status==1){
										//energy_map[event_id].first.push_back(particle_id);		  
										//energy_map[event_id].second.push_back(energy);		  
										//std::cout << "event_id = " << event_id << std::endl;
										//std::cout << "energy = " << energy << std::endl;
										//std::cout << "particle_id = " << particle_id << std::endl;
										if( eventid_vec.size()==0 ){
										//std::cout << "event_id = " << event_id << std::endl;
										//std::cout << "energy = " << energy << std::endl;
										//std::cout << "particle_id = " << particle_id << std::endl;
														eventid_vec.push_back(event_id);
														energy_vec.push_back(energy);
														particleid_vec.push_back(particle_id);
										}
										else{
														if(event_id == eventid_vec.back() && (abs(particle_id - particleid_vec.back())==1) ){
										//std::cout << "event_id = " << event_id << std::endl;
										//std::cout << "energy = " << energy << std::endl;
										//std::cout << "particle_id = " << particle_id << std::endl;
																		energy_vec.push_back(energy);
																		Hits_Esum_gg->Fill( abs(energy_vec.at(0) - energy_vec.at(1) ));
														}
														else{
																		eventid_vec.clear();
																		energy_vec.clear();
																		particleid_vec.clear();
														}
										}
									}
								}
				}
				else if( file_iterator >= 4){
					Hits_P_T_pairs->Fill(sqrt(mom_x*mom_x+mom_y*mom_y),weight_bunches/float(NUMBER_OF_pairFILES));
				}
			}
						file->Close();
		}
	}
	//gStyle->SetOptStat(11);

	//Plot the histogram and save it
	std::vector< TH1D* > histos;
	histos.push_back(Hits_P_T_pairs);
	histos.push_back(Hits_P_T_gg);
	double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
	for(size_t iterator = 0; iterator < histos.size(); ++iterator){
					histos.at(iterator)->SetMinimum(0.5);
					histos.at(iterator)->SetMaximum(max);
	}

	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);

	canvas1->cd();
	canvas1->SetLogy(1);

	Hits_P_T_pairs->SetMarkerColor(0+1);
	Hits_P_T_pairs->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec;
	st_vec.push_back(new TPaveStats());
	st_vec.at(0) = (TPaveStats*)Hits_P_T_pairs->GetListOfFunctions()->FindObject("stats");
	st_vec.at(0)->SetLineColor(0+1);
	st_vec.at(0)->SetX1NDC(0.65); //new x start position
	st_vec.at(0)->SetX2NDC(0.85); //new x end position
	st_vec.at(0)->SetY1NDC(0.83); //new x start position
	st_vec.at(0)->SetY2NDC(0.9); //new x end position
	float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

	Hits_P_T_gg->SetMarkerColor(2);
	Hits_P_T_gg->SetLineColor(2);
	Hits_P_T_gg->Draw("SAMES");
	canvas1->Update();
	st_vec.push_back(new TPaveStats());
	st_vec.at(1)= (TPaveStats*)Hits_P_T_gg->GetListOfFunctions()->FindObject("stats");
	st_vec.at(1)->SetLineColor(2);
	st_vec.at(1)->SetX1NDC(0.65); //new x start position
	st_vec.at(1)->SetX2NDC(0.85); //new x end position
	st_vec.at(1)->SetY2NDC(st_vec.at(0)->GetY1NDC()); //new x end position
	st_vec.at(1)->SetY1NDC(st_vec.at(1)->GetY2NDC()-boxsize); //new x start position

	//TLegend *leg = new TLegend(0.65,0.76,0.85,0.9);
	//leg->AddEntry((TObject*)0,"P_T_pairs","C");
	//stringstream entriestimepairs;
	//entriestimepairs << "Entries = " << Hits_P_T_pairs->GetEntries();
	//leg->AddEntry((TObject*)0,entriestimepairs.str().c_str(),"C");
	//leg->AddEntry((TObject*)0,"P_T_#gamma#gamma#rightarrow hadrons","C");
	//stringstream entriestimegg;
	//entriestimegg << "Entries = " << Hits_P_T_gg->Integral();
	//leg->AddEntry((TObject*)0,entriestimegg.str().c_str(),"C");
	//
	//TList *list = (TList*)leg->GetListOfPrimitives();
	//TIter next(list);
	//TObject* object = 0;
	//while ((object = next())){
	//	TLegendEntry * entry = (TLegendEntry*)object;
	//	entry->SetTextAlign(22);
	//}
	//leg->Draw();

	TPaveText *text1 = new TPaveText(0.65,0.83,0.85,0.9,"brNDC");
	text1->SetTextFont(62);
	text1->SetTextColor(1);
	text1->SetFillColor(0);
	text1->AddText("P_T_pairs");
	text1->AddLine(0,0.5,1,0.5);
	std::stringstream entriespairs;
	entriespairs << "Entries = " << Hits_P_T_pairs->Integral();
	text1->AddText(entriespairs.str().c_str());
	TPaveText *text2 = new TPaveText(0.65,0.76,0.85,0.83,"brNDC");
	text2->SetTextFont(62);
	text2->SetTextColor(2);
	text2->SetFillColor(0);
	//text->AddLine(0,0.5,1,0.5);
	text2->AddText("P_T_#gamma#gamma#rightarrow hadrons");
	text2->AddLine(0,0.5,1,0.5);
	std::stringstream entriesgg;
	entriesgg << "Entries = " << Hits_P_T_gg->Integral();
	text2->AddText(entriesgg.str().c_str());
	text1->Draw();
	text2->Draw();

	canvas1->Print(("output/gg-hadrons_pairs_comparison_elecposi_leavindetector_PT_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_pairs_comparison_elecposi_leavindetector_PT_"+subdetectornames+".cxx").c_str());

	Hits_P_x_gg->Draw();
	canvas1->Print(("output/gg-hadrons_elecposi_leavindetector_Px_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_elecposi_leavindetector_Px_"+subdetectornames+".cxx").c_str());
	Hits_P_y_gg->Draw();
	canvas1->Print(("output/gg-hadrons_elecposi_leavindetector_Py_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_elecposi_leavindetector_Py_"+subdetectornames+".cxx").c_str());
	Hits_P_z_gg->Draw();
	canvas1->Print(("output/gg-hadrons_elecposi_leavindetector_Pz_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_elecposi_leavindetector_Pz_"+subdetectornames+".cxx").c_str());
	if( Hits_Esum_gg->GetEntries() != 0 ){
					Hits_Esum_gg->Draw();
					canvas1->Print(("output/gg-hadrons_elecposi_leavingdetector_Esum_"+subdetectornames+".pdf").c_str());
					canvas1->Print(("output/gg-hadrons_elecposi_leavingdetector_Esum_"+subdetectornames+".cxx").c_str());
	}

	return 0;
}

