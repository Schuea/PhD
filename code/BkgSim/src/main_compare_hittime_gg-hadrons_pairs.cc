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
				bool pair_inputfile_set = false;
				bool ewpw_inputfile_set = false;
				bool ewpb_inputfile_set = false;
				bool ebpw_inputfile_set = false;
				bool ebpb_inputfile_set = false;
				bool subdetector_set = false;

				for (int i = 1; i < argc; i++) {
								if (argv[i] == std::string("-n")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-s") 
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
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL) {
																int j = 1;
																do {
																				if (argv[i + j] != std::string("-n") && argv[i + j] != std::string("-s")
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
								if (argv[i] == std::string("-ewpw") && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
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
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
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
																&& argv[i + 1] != std::string("-i") 
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ewpw") 
																&& argv[i + 1] != std::string("-ebpb")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
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
																&& argv[i + 1] != std::string("-ewpb") 
																&& argv[i + 1] != std::string("-ebpw") 
																&& argv[i + 1] != std::string("-ewpw")) {
												if (argv[i + 1] != NULL && argv[i + 1] != std::string("-n") && argv[i + 1] != std::string("-s")
																				&& argv[i + 1] != std::string("-i") 
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
				if (!pair_inputfile_set || !ewpw_inputfile_set || !ewpb_inputfile_set || !ebpw_inputfile_set || !ebpb_inputfile_set || !subdetector_set || !NUMBER_OF_pairFILES_set) {
								std::cerr
												<< "You didn't give the name for the subdector, the inputfiles or the amount of files. Please try again!"
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
	std::vector<float> range_array;

	//Make histogram for storing the information
	std::string const histo_title_creationtime = "Creation time of particles from #gamma#gamma #rightarrow hadrons events and pair background hitting the " + subdetectornames + ";Hit time [ns]; #hits";
	TH1D* Hits_Time_gg = new TH1D("Time_gammagamma_hadrons",histo_title_creationtime.c_str(),20,0,100);
	TH1D* Hits_Time_pairs = new TH1D("Time_pairs",histo_title_creationtime.c_str(),20,0,100);

	for (size_t subdetector_iterator = 0; subdetector_iterator < SubDetectors->size(); ++subdetector_iterator) {

		for (size_t file_iterator = 0; file_iterator < inputfiles->size(); ++file_iterator) {
		
			TFile *file = TFile::Open(inputfiles->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, SubDetectors->at(subdetector_iterator)->GetName());

			//Set the branches
			float actualtime = 0.0;
			float creationtime = 0.0;
			int event_id = -999;
			tree->SetBranchStatus("*", 0);
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
			} else {
							std::cerr << "The given TTree name does not match any TTree in the inputfile!" << std::endl;
							std::terminate();
			}

			//Now we loop through the tree

			long long int const entries = tree->GetEntries();
			double weight(1);
			if( file_iterator == 0) weight = 291.0/10000.0;
			if( file_iterator == 1) weight = 501.0/10000.0;
			if( file_iterator == 2) weight = 502.0/10000.0;
			if( file_iterator == 3) weight = 829.0/10000.0;
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				if( file_iterator < 4){
								Hits_Time_gg->Fill(actualtime,weight);
				}
				//if( file_iterator == 0	&& event_id <= 291){ //ewpw
				//	Hits_Time_gg->Fill(actualtime); 
				//}
				//else if( file_iterator == 1 && event_id <= 501){ //ewpb
				//	Hits_Time_gg->Fill(actualtime); 
				//}
				//else if( file_iterator == 2 && event_id <= 502){ //ebpw
				//	Hits_Time_gg->Fill(actualtime); 
				//}
				//else if( file_iterator == 3 && event_id <= 829){ //ebpb
				//	Hits_Time_gg->Fill(actualtime); 
				//}
				else if( file_iterator >= 4){
					Hits_Time_pairs->Fill(actualtime);
				}
			}
			file->Close();
		}
	}
	gStyle->SetOptStat(11);

	//Plot the histogram and save it
	Hits_Time_pairs->SetMinimum(0.5);
	Hits_Time_gg->SetMinimum(0.5);

	
	TCanvas *canvas1 = new TCanvas("canvas1", "canvas", 800, 600);

	canvas1->cd();
	canvas1->SetLogy(1);
	Hits_Time_pairs->SetMarkerColor(0+1);
	Hits_Time_pairs->Draw();
	canvas1->Update();
	std::vector<TPaveStats*> st_vec;
	st_vec.push_back(new TPaveStats());
	st_vec.at(0) = (TPaveStats*)Hits_Time_pairs->GetListOfFunctions()->FindObject("stats");
	st_vec.at(0)->SetLineColor(0+1);
	st_vec.at(0)->SetX1NDC(0.65); //new x start position
	st_vec.at(0)->SetX2NDC(0.85); //new x end position
	st_vec.at(0)->SetY1NDC(0.83); //new x start position
	st_vec.at(0)->SetY2NDC(0.9); //new x end position
	float boxsize = st_vec.at(0)->GetY2NDC()-st_vec.at(0)->GetY1NDC();

	Hits_Time_gg->SetMarkerColor(2);
	Hits_Time_gg->SetLineColor(2);
	Hits_Time_gg->Draw("SAMES");
	canvas1->Update();
	st_vec.push_back(new TPaveStats());
	st_vec.at(1)= (TPaveStats*)Hits_Time_gg->GetListOfFunctions()->FindObject("stats");
	st_vec.at(1)->SetLineColor(2);
	st_vec.at(1)->SetX1NDC(0.65); //new x start position
	st_vec.at(1)->SetX2NDC(0.85); //new x end position
	st_vec.at(1)->SetY2NDC(st_vec.at(0)->GetY1NDC()); //new x end position
	st_vec.at(1)->SetY1NDC(st_vec.at(1)->GetY2NDC()-boxsize); //new x start position

	//TLegend *leg = new TLegend(0.65,0.76,0.85,0.9);
	//leg->AddEntry((TObject*)0,"Time_pairs","C");
	//stringstream entriestimepairs;
	//entriestimepairs << "Entries = " << Hits_Time_pairs->GetEntries();
	//leg->AddEntry((TObject*)0,entriestimepairs.str().c_str(),"C");
	//leg->AddEntry((TObject*)0,"Time_#gamma#gamma#rightarrow hadrons","C");
	//stringstream entriestimegg;
	//entriestimegg << "Entries = " << Hits_Time_gg->Integral();
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
	text1->AddText("Time_pairs");
	text1->AddLine(0,0.5,1,0.5);
	std::stringstream entriestimepairs;
	entriestimepairs << "Entries = " << Hits_Time_pairs->GetEntries();
	text1->AddText(entriestimepairs.str().c_str());
	TPaveText *text2 = new TPaveText(0.65,0.76,0.85,0.83,"brNDC");
	text2->SetTextFont(62);
	text2->SetTextColor(2);
	text2->SetFillColor(0);
	//text->AddLine(0,0.5,1,0.5);
	text2->AddText("Time_#gamma#gamma#rightarrow hadrons");
	text2->AddLine(0,0.5,1,0.5);
	std::stringstream entriestimegg;
	entriestimegg << "Entries = " << Hits_Time_gg->Integral();
	text2->AddText(entriestimegg.str().c_str());
	text1->Draw();
	text2->Draw();

	canvas1->Print(("output/gg-hadrons_pairs_comparison_time_"+subdetectornames+".pdf").c_str());
	canvas1->Print(("output/gg-hadrons_pairs_comparison_time_"+subdetectornames+".cxx").c_str());
	
	return 0;
}

