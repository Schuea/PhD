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

	std::vector<std::string> *ewpw_inputfilenames = new std::vector<std::string>();
	std::vector<std::string> *ewpb_inputfilenames = new std::vector<std::string>();
	std::vector<std::string> *ebpw_inputfilenames = new std::vector<std::string>();
	std::vector<std::string> *ebpb_inputfilenames = new std::vector<std::string>();
	std::vector<std::string> argument_subdetectors;

	int NUMBER_OF_FILES = 0;
	bool NUMBER_OF_FILES_set = false;
	bool ewpw_inputfile_set = false;
	bool ewpb_inputfile_set = false;
	bool ebpw_inputfile_set = false;
	bool ebpb_inputfile_set = false;
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
		if (argv[i] == std::string("-ewpw") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
											&& argv[i + 1] != std::string("-ewpb") 
											&& argv[i + 1] != std::string("-ebpw") 
											&& argv[i + 1] != std::string("-ebpb")) {
			if (argv[i + 1] != NULL) {
				int j = 1;
				do {
					if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
											&& argv[i + j] != std::string("-ewpb") 
											&& argv[i + j] != std::string("-ebpw") 
											&& argv[i + j] != std::string("-ebpb")) {
						ewpw_inputfilenames->push_back(argv[i + j]);
						++j;
					} else {
						break;
					}
				} while (j <= NUMBER_OF_FILES);
				ewpw_inputfile_set = true;
			} else {
				std::cerr << "You didn't give an argument for the ewpw inputfile(s)!" << std::endl;
			}
		}
		if (argv[i] == std::string("-ewpb") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
											&& argv[i + 1] != std::string("-ewpw") 
											&& argv[i + 1] != std::string("-ebpw") 
											&& argv[i + 1] != std::string("-ebpb")) {
			if (argv[i + 1] != NULL) {
				int j = 1;
				do {
					if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
											&& argv[i + j] != std::string("-ewpw") 
											&& argv[i + j] != std::string("-ebpw") 
											&& argv[i + j] != std::string("-ebpb")) {
						ewpb_inputfilenames->push_back(argv[i + j]);
						++j;
					} else {
						break;
					}
				} while (j <= NUMBER_OF_FILES);
				ewpb_inputfile_set = true;
			} else {
				std::cerr << "You didn't give an argument for the ewpb inputfile(s)!" << std::endl;
			}
		}
		if (argv[i] == std::string("-ebpw") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
											&& argv[i + 1] != std::string("-ewpb") 
											&& argv[i + 1] != std::string("-ewpw") 
											&& argv[i + 1] != std::string("-ebpb")) {
			if (argv[i + 1] != NULL) {
				int j = 1;
				do {
					if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
											&& argv[i + j] != std::string("-ewpb") 
											&& argv[i + j] != std::string("-ewpw") 
											&& argv[i + j] != std::string("-ebpb")) {
						ebpw_inputfilenames->push_back(argv[i + j]);
						++j;
					} else {
						break;
					}
				} while (j <= NUMBER_OF_FILES);
				ebpw_inputfile_set = true;
			} else {
				std::cerr << "You didn't give an argument for the ebpw inputfile(s)!" << std::endl;
			}
		}
		if (argv[i] == std::string("-ebpb") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
											&& argv[i + 1] != std::string("-ewpb") 
											&& argv[i + 1] != std::string("-ebpw") 
											&& argv[i + 1] != std::string("-ewpw")) {
			if (argv[i + 1] != NULL) {
				int j = 1;
				do {
					if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
											&& argv[i + j] != std::string("-ewpb") 
											&& argv[i + j] != std::string("-ebpw") 
											&& argv[i + j] != std::string("-ewpw")) {
						ebpb_inputfilenames->push_back(argv[i + j]);
						++j;
					} else {
						break;
					}
				} while (j <= NUMBER_OF_FILES);
				ebpb_inputfile_set = true;
			} else {
				std::cerr << "You didn't give an argument for the ebpb inputfile(s)!" << std::endl;
			}
		}
		if (argv[i] == std::string("-s")) {
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") 
											&& argv[i + 1] != std::string("-ewpw") 
											&& argv[i + 1] != std::string("-ewpb") 
											&& argv[i + 1] != std::string("-ebpw") 
											&& argv[i + 1] != std::string("-ebpb")) {
				if (argv[i + 1] != NULL) {
					int j = 1;
					do {
						argument_subdetectors.push_back(argv[i + j]);
						++j;
					} while (argv[i + j] != NULL && argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
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
	if (!ewpw_inputfile_set || !ewpb_inputfile_set || !ebpw_inputfile_set || !ebpb_inputfile_set || !subdetector_set || !NUMBER_OF_FILES_set) {
		std::cerr
				<< "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
				<< std::endl;
		exit(1);
	}

	std::vector<std::string> *inputfiles = new std::vector<std::string>();
	inputfiles->insert(inputfiles->begin(),ewpw_inputfilenames->begin(),ewpw_inputfilenames->end());
	inputfiles->insert(inputfiles->end(),  ewpb_inputfilenames->begin(),ewpb_inputfilenames->end());
	inputfiles->insert(inputfiles->end(),  ebpw_inputfilenames->begin(),ebpw_inputfilenames->end());
	inputfiles->insert(inputfiles->end(),  ebpb_inputfilenames->begin(),ebpb_inputfilenames->end());
	

	std::vector<Subdetector*> * SubDetectors = new std::vector<Subdetector*>();
	std::string* subdetector_name = new std::string("");
	SetupSubDetectorsVector(SubDetectors, subdetector_name, argument_subdetectors);
	std::string subdetectornames = (*subdetector_name);
	std::vector<float> range_array;

	//Make histogram for storing the information
	//std::vector< TH2D* > Hits_Time_rtime_2D_;
	range_array = SubDetectors->at(0)->GetROOTHisto_binning3D(); //SOMETHING HARD CODED!!
	float rmax = 0.7*sqrt(pow(range_array[5], 2) + pow(range_array[8], 2));//Make the plot a big smaller than the actual measurement of the subdetector
	float rmin = 0.;
	float rrange = rmax - rmin;
	float timemax = 0.0; //ns (one bunch spacing is 554 ns)
	if(NUMBER_OF_FILES>1) timemax = NUMBER_OF_FILES * 600.0; //ns (one bunch spacing is 554 ns)
	else timemax = 70.0; //ns (one bunch spacing is 554 ns)
	float timemin = 0.;
	float timerange = timemax - timemin;
	std::array<float, 6> axis_vector = { int(timerange), timemin, timemax, int(rrange), rmin, rmax };

	std::string const histo_title1D = "Hit time for particles from #gamma#gamma #rightarrow hadrons events hitting the " + subdetectornames + ";Hit time [ns];Hits";
	std::string const histo_title2D = "Radial position of hits over hit time for particles from #gamma#gamma #rightarrow hadrons events hitting the " + subdetectornames + ";Hit time [ns];r [mm]";
	std::vector<TH2D*> Hits_Time_rtime_;
	std::vector<TH1D*> Hits_Time_;

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {
		Time PassedTime;

		for (size_t file_iterator = 0; file_iterator < inputfiles->size(); ++file_iterator) {
			std::cout << inputfiles->at(file_iterator) << std::endl;
		
			if (NUMBER_OF_FILES > 1 && file_iterator % NUMBER_OF_FILES == 0){
							std::stringstream histoname1D, histoname2D;
							histoname1D << "HitTime_" << subdetectornames << "_train#" << file_iterator+1;
							histoname2D << "HitTime_radialPos_" << subdetectornames << "_train#" << file_iterator+1;
							std::string histo_name1D = histoname1D.str();
							std::string histo_name2D = histoname2D.str();
							Hits_Time_.emplace_back(new TH1D(histo_name1D.c_str(),  histo_title1D.c_str(),20,0,100));
							Hits_Time_rtime_.emplace_back(new TH2D(histo_name2D.c_str(), histo_title2D.c_str(), 
																			axis_vector.at(0), axis_vector.at(1), axis_vector.at(2),
																			axis_vector.at(3), axis_vector.at(4), axis_vector.at(5)));
			}
			else if (file_iterator == 0){//For one bunch only, emplace only one histogram back
							std::stringstream histoname1D, histoname2D;
							histoname1D << "HitTime_" << subdetectornames << "_train#" << 1;
							histoname2D << "HitTime_radialPos_" << subdetectornames << "_train#" << 1;
							std::string histo_name1D = histoname1D.str();
							std::string histo_name2D = histoname2D.str();
							Hits_Time_.emplace_back(new TH1D(histo_name1D.c_str(),  histo_title1D.c_str(),20,0,100));
							Hits_Time_rtime_.emplace_back(new TH2D(histo_name2D.c_str(), histo_title2D.c_str(), 
																			axis_vector.at(0), axis_vector.at(1), axis_vector.at(2),
																			axis_vector.at(3), axis_vector.at(4), axis_vector.at(5)));
			}
			TFile *file = TFile::Open(inputfiles->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			float actualtime = 0.0;
			float creationtime = 0.0;
			float x_cal = 0.0;
			float y_cal = 0.0;
			double x_tracker = 0.0;
			double y_tracker = 0.0;
			int event_id = -999;
			tree->SetBranchStatus("*", 0);
			tree->SetBranchStatus("HitPosition_x", kTRUE); 
			tree->SetBranchStatus("HitPosition_y", kTRUE);
			tree->SetBranchStatus("event_id", kTRUE);
			tree->SetBranchAddress("event_id", &event_id);
			
			if (tree->GetName() == std::string("Tree_EcalBarrel") 
											|| tree->GetName() == std::string("Tree_EcalEndcap")
											|| tree->GetName() == std::string("Tree_HcalBarrel") 
											|| tree->GetName() == std::string("Tree_HcalEndcap")
											|| tree->GetName() == std::string("Tree_MuonBarrel") 
											|| tree->GetName() == std::string("Tree_MuonEndcap")
											|| tree->GetName() == std::string("Tree_BeamCal") 
											|| tree->GetName() == std::string("Tree_LumiCal")) {
							tree->SetBranchStatus("HitMotherCreationTime", 1);
							tree->SetBranchAddress("HitMotherCreationTime", &creationtime);
							tree->SetBranchStatus("HitContrTime", 1);
							tree->SetBranchAddress("HitContrTime", &actualtime);
							tree->SetBranchAddress("HitPosition_x", &x_cal);
							tree->SetBranchAddress("HitPosition_y", &y_cal);
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
							tree->SetBranchStatus("HitParticleCreationTime", 1);
							tree->SetBranchAddress("HitParticleCreationTime", &creationtime);
							tree->SetBranchStatus("HitTime", 1);
							tree->SetBranchAddress("HitTime", &actualtime);
							tree->SetBranchAddress("HitPosition_x", &x_tracker);
							tree->SetBranchAddress("HitPosition_y", &y_tracker);
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}
			float x = 0.0;
			float y = 0.0;

			//Figure out which bunch in which train we are, in order to calculate the time that's already passed:
			std::pair<int, int> Number_train_bunch(0,0);
			if(NUMBER_OF_FILES == 1) Number_train_bunch = Set_train_bunch_number(0);
			else{
							if((int)file_iterator%2 == 0) Number_train_bunch = Set_train_bunch_number(0);
							else Number_train_bunch = Set_train_bunch_number(1);
			}
			PassedTime.Calculate_passedbytime(Number_train_bunch.first, Number_train_bunch.second);//first: number of train, second: number of bunch
  		float absolutetime = 0.0;

			//Have the markers for different bunches coloured differently
			Hits_Time_rtime_.at(file_iterator % NUMBER_OF_FILES)->SetMarkerColor(Number_train_bunch.second);

			//Now we loop through the tree

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				if (tree->GetName() == std::string("Tree_EcalBarrel") 
											|| tree->GetName() == std::string("Tree_EcalEndcap")
											|| tree->GetName() == std::string("Tree_HcalBarrel") 
											|| tree->GetName() == std::string("Tree_HcalEndcap")
											|| tree->GetName() == std::string("Tree_MuonBarrel") 
											|| tree->GetName() == std::string("Tree_MuonEndcap")
											|| tree->GetName() == std::string("Tree_BeamCal") 
											|| tree->GetName() == std::string("Tree_LumiCal")) {
								x = x_cal;
								y = y_cal;
			}
			else if (tree->GetName() == std::string("Tree_SiVertexBarrel")
											|| tree->GetName() == std::string("Tree_SiVertexEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerBarrel")
											|| tree->GetName() == std::string("Tree_SiTrackerEndcap")
											|| tree->GetName() == std::string("Tree_SiTrackerForward")) {
								x = x_tracker;
								y = y_tracker;
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}
			//absolute time = time in respect to the current bunch interaction + time passed by since first bunch interaction
			  absolutetime = actualtime + PassedTime.Get_passedbytime();
				//double w = 1.0;
				//if((int)file_iterator < 1*NUMBER_OF_FILES) w = 0.029158296; //ewpw
				//else if(1*NUMBER_OF_FILES <= (int)file_iterator && (int)file_iterator < 2*NUMBER_OF_FILES) w = 0.033760669; //ewpb
				//else if(2*NUMBER_OF_FILES <= (int)file_iterator && (int)file_iterator < 3*NUMBER_OF_FILES) w = 0.033777845; //ebpw
				//else if(3*NUMBER_OF_FILES <= (int)file_iterator && (int)file_iterator < 4*NUMBER_OF_FILES) w = 0.0483693499; //ebpb
				if((int)file_iterator < 1*NUMBER_OF_FILES  
												&& event_id <= 291){ //ewpw
					Hits_Time_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime); 
					Hits_Time_rtime_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime, sqrt(pow(x, 2) + pow(y, 2)));
				}
				else if(1*NUMBER_OF_FILES <= (int)file_iterator && (int)file_iterator < 2*NUMBER_OF_FILES
												&& event_id <= 501){ //ewpb
					Hits_Time_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime); 
					Hits_Time_rtime_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime, sqrt(pow(x, 2) + pow(y, 2)));
				}
				else if(2*NUMBER_OF_FILES <= (int)file_iterator && (int)file_iterator < 3*NUMBER_OF_FILES
												&& event_id <= 502){ //ebpw
					Hits_Time_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime); 
					Hits_Time_rtime_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime, sqrt(pow(x, 2) + pow(y, 2)));
				}
				else if(3*NUMBER_OF_FILES <= (int)file_iterator && (int)file_iterator < 4*NUMBER_OF_FILES
												&& event_id <= 829){ //ebpb
					Hits_Time_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime); 
					Hits_Time_rtime_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime, sqrt(pow(x, 2) + pow(y, 2)));
				}
				//The first file of every process is from bunch 1 and needs to be filled into histo 0, the second file of every process need to be filled into histo 1 and so on...: 
				//Hits_Time_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime,w); 
				//Hits_Time_rtime_.at(file_iterator % NUMBER_OF_FILES)->Fill(absolutetime, sqrt(pow(x, 2) + pow(y, 2)),w);
			}
			file->Close();
		}
	}
	gStyle->SetOptStat(11);

	//Plot the histogram and save it
	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);

	canvas1->cd();
	canvas1->SetLogy(0);
	//histo1->Draw("colz");
	
	if(NUMBER_OF_FILES==1){
					canvas1->SetLogz(1);
					Hits_Time_rtime_.at(0)->Draw("colz");
	}
	else Hits_Time_rtime_.at(0)->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec;
	st_vec.push_back(new TPaveStats());
	st_vec.at(0) = (TPaveStats*)Hits_Time_rtime_.at(0)->GetListOfFunctions()->FindObject("stats");
	st_vec.at(0)->SetLineColor(0+1);
	st_vec.at(0)->SetX1NDC(0.65); //new x start position
	st_vec.at(0)->SetX2NDC(0.85); //new x end position
	st_vec.at(0)->SetY1NDC(0.83); //new x start position
	st_vec.at(0)->SetY2NDC(0.9); //new x end position
	float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

	for (size_t iterator = 1; iterator < Hits_Time_rtime_.size(); ++iterator) {
		if(NUMBER_OF_FILES==1) Hits_Time_rtime_.at(iterator)->Draw("colz,SAMES");
		else Hits_Time_rtime_.at(iterator)->Draw("SAMES");
		canvas1->Update();
		st_vec.push_back(new TPaveStats());
		st_vec.at(iterator)= (TPaveStats*)Hits_Time_rtime_.at(iterator)->GetListOfFunctions()->FindObject("stats");
		st_vec.at(iterator)->SetLineColor(iterator+1);
		st_vec.at(iterator)->SetX1NDC(0.65); //new x start position
		st_vec.at(iterator)->SetX2NDC(0.85); //new x end position
		st_vec.at(iterator)->SetY2NDC(st_vec.at(iterator-1)->GetY1NDC()); //new x end position
		st_vec.at(iterator)->SetY1NDC(st_vec.at(iterator)->GetY2NDC()-boxsize); //new x start position
	}


	canvas1->Print(("output/gg-hadrons_hittime_radialpos_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_hittime_radialpos_"+subdetectornames+".cxx").c_str());

	canvas1->SetLogy(1);
	Hits_Time_.at(0)->SetMarkerColor(0+1);
	Hits_Time_.at(0)->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec2;
	st_vec2.push_back(new TPaveStats());
	st_vec2.at(0) = (TPaveStats*)Hits_Time_.at(0)->GetListOfFunctions()->FindObject("stats");
	st_vec2.at(0)->SetLineColor(0+1);
	st_vec2.at(0)->SetX1NDC(0.65); //new x start position
	st_vec2.at(0)->SetX2NDC(0.85); //new x end position
	st_vec2.at(0)->SetY1NDC(0.83); //new x start position
	st_vec2.at(0)->SetY2NDC(0.9); //new x end position
	float boxsize2 = st_vec2.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

	for (size_t iterator = 1; iterator < Hits_Time_.size(); ++iterator) {
		Hits_Time_.at(iterator)->SetMarkerColor(iterator+1);
		Hits_Time_.at(iterator)->SetLineColor(iterator+1);
		Hits_Time_.at(iterator)->Draw("SAMES");
		canvas1->Update();
		st_vec2.push_back(new TPaveStats());
		st_vec2.at(iterator)= (TPaveStats*)Hits_Time_.at(iterator)->GetListOfFunctions()->FindObject("stats");
		st_vec2.at(iterator)->SetLineColor(iterator+1);
		st_vec2.at(iterator)->SetX1NDC(0.65); //new x start position
		st_vec2.at(iterator)->SetX2NDC(0.85); //new x end position
		st_vec2.at(iterator)->SetY2NDC(st_vec2.at(iterator-1)->GetY1NDC()); //new x end position
		st_vec2.at(iterator)->SetY1NDC(st_vec2.at(iterator)->GetY2NDC()-boxsize2); //new x start position
	}

	canvas1->Print(("output/gg-hadrons_hittime_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_hittime_"+subdetectornames+".cxx").c_str());
	
	return 0;
}

