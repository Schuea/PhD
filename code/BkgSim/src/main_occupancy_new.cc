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

long long int MakeNewCellID(double const x, double const y, Subdetector const & component);
void Draw_single_plots ( TH1D* histo, TCanvas* canvas, bool normalize);
void Draw_multiple_plots ( std::vector< TH1D* > histos, TCanvas* canvas, bool normalize);

float weight = 1;

int main(int const argc, char const * const * const argv) {
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
        } while (j <= NUMBER_OF_FILES);
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

	CellHits * HitCount = new CellHits( &det );
		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
			TFile *file = TFile::Open(inputfilenames->at(file_iterator).c_str());
			TTree *tree = Get_TTree(file, det.getName());

			//Set the branches
			tree->SetBranchStatus("*", 0);
			//tree->SetBranchStatus("HitCellID", 1);
			tree->SetBranchStatus("HitCellID0", 1);
			//tree->SetBranchStatus("HitCellID1", 1);
			tree->SetBranchStatus("HitPosition_x", 1);
			tree->SetBranchStatus("HitPosition_y", 1);
			tree->SetBranchStatus("HitPosition_z", 1);

			int HitCellID0(0);
      //int HitCellID1(0);
			//float HitPosition_x(0.0), HitPosition_y(0.0), HitPosition_z(0.0);
			double HitPosition_x(0.0), HitPosition_y(0.0), HitPosition_z(0.0);

			//tree->SetBranchAddress("HitCellID", &HitCellID0);
			tree->SetBranchAddress("HitCellID0", &HitCellID0);
			//tree->SetBranchAddress("HitCellID1", &HitCellID1);
			tree->SetBranchAddress("HitPosition_x", &HitPosition_x);
			tree->SetBranchAddress("HitPosition_y", &HitPosition_y);
			tree->SetBranchAddress("HitPosition_z", &HitPosition_z);

			//Now we loop through the tree
			//Combine the two Cell ID's into a single new Cell ID
			//See how often the new Cell ID occurs in total, this is the occupancy

			long long int const entries = tree->GetEntries();
			for (long long int i = 0; i < entries; ++i) {
				tree->GetEntry(i);
				//if (HitPosition_z < 0) continue;
				//Make a combined cell ID
        long long int HitCellID1 = MakeNewCellID(HitPosition_x,HitPosition_y,det);
				long long int const combined_cell_id = (long long) HitCellID1 << 32 | HitCellID0;
				//Use the CellHits class for storing the hit cells and their hitcounts
				HitCount->Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y, HitPosition_z);
			}
			file->Close();
		}

	//Make histogram for storing the information
  std::string subdetectorname = det.getName();
	std::string const title = "Occupancy for subdetector " + subdetectorname + ";Number of hits per cell;Number of cells";
	std::string const title2 = "Occupancy for subdetector " + subdetectorname + " wrt to total number of cells;Number of hits per cell ;Number of cells";
	std::string const title3 = "Number of hits lost for a given buffer depth for subdetector " + subdetectorname + ";Assumend buffer depth;Number of hits lost";
	std::string const title4 = "Number of dead cells for a given buffer depth for subdetector " + subdetectorname + ";Assumend buffer depth;Number of dead cells";
	std::vector< TH1D* > histos;
	std::vector< TH1D* > histos_numcells;
	std::vector< TH1D* > histos_bufferdepth;
	std::vector< TH1D* > histos_deadcells;

	int tot_no_hits = 0;
	int max_no_hits = 0;
  for (size_t vecpos = 0; vecpos < HitCount->Get_HitCount().size(); ++vecpos) {
			if (HitCount->Get_HitCount().at(vecpos) > max_no_hits) max_no_hits = HitCount->Get_HitCount().at(vecpos);
      tot_no_hits += HitCount->Get_HitCount().at(vecpos);
	}

	int xrange = max_no_hits + max_no_hits/10;
	TH1D* All_Layers_histo = new TH1D("All_layers", title.c_str(), xrange, 0, xrange);
	TH1D* All_Layers_histo_numcells = new TH1D("All_layers_wrt_#cells", title2.c_str(), xrange, 0, xrange);
	TH1D* All_Layers_histo_bufferdepth = new TH1D("All_layers_bufferdepth", title3.c_str(), xrange, 0, xrange);
	TH1D* All_Layers_histo_deadcells = new TH1D("All_layers_deadcells", title4.c_str(), xrange, 0, xrange);
	int max_num_layers = det.getNumberOfLayers();
	for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
		std::stringstream layername, layername2, layername3, layername4;
		layername << "Layer_" << number_layer;
		layername2 << "Layer_" << number_layer << "_numcells";
		layername3 << "Layer_" << number_layer << "_bufferdepth";
		layername4 << "Layer_" << number_layer << "_deadcells";
		histos.emplace_back(new TH1D(layername.str().c_str(), title.c_str(), xrange, 0, xrange));
		histos_numcells.emplace_back(new TH1D(layername2.str().c_str(), title2.c_str(), xrange, 0, xrange));
		histos_bufferdepth.emplace_back(new TH1D(layername3.str().c_str(), title3.c_str(), xrange, 0, xrange));
		histos_deadcells.emplace_back(new TH1D(layername4.str().c_str(), title4.c_str(), xrange, 0, xrange));
	}
  for (size_t vecpos = 0; vecpos < HitCount->Get_HitCount().size(); ++vecpos) {
    if(HitCount->Get_HitCount().at(vecpos) > 0){
			//std::cout << "Layer: " << HitCount->Get_Layer().at(vecpos) << std::endl;
      histos.at(HitCount->Get_Layer().at(vecpos) -1 )->Fill(HitCount->Get_HitCount().at(vecpos),weight);//-1 for SiTrackerBarrel, because layer count starts from 1
      All_Layers_histo->Fill(HitCount->Get_HitCount().at(vecpos),weight);
    }
  }
	//Filling numcells plots:
	long long int tot_num_cells = 0;
	long long int tot_num_hitcells = 0;
	for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
					tot_num_cells += det.getNumberOfCells().at(number_layer);
					for (int bin = 2; bin < histos.at(number_layer)->GetNbinsX(); ++bin) {
									tot_num_hitcells += histos.at(number_layer)->GetBinContent(bin);
									histos_numcells.at(number_layer)->SetBinContent(bin, histos.at(number_layer)->GetBinContent(bin) );
					}
					histos.at( number_layer )->SetBinContent(1, det.getNumberOfCells().at(number_layer) - tot_num_hitcells );
					histos_numcells.at( number_layer )->SetBinContent(1, det.getNumberOfCells().at(number_layer));
					tot_num_hitcells = 0;
	}
	for (int bin = 2; bin < All_Layers_histo->GetNbinsX(); ++bin) {
					tot_num_hitcells += All_Layers_histo->GetBinContent(bin);
					All_Layers_histo_numcells->SetBinContent(bin, All_Layers_histo->GetBinContent(bin));
	}
	All_Layers_histo->SetBinContent(1, tot_num_cells - tot_num_hitcells);
	All_Layers_histo_numcells->SetBinContent(1, tot_num_cells);

	//Filling bufferdepth plots:
	for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
					for (int i = 0; i <= max_no_hits; ++i){//For each bufferdepth
									int tot = 0;
									int deadcells = 0;
									for (int bin = i+1; bin < histos.at(number_layer)->GetNbinsX(); ++bin) {//go through the histo from bufferdepth value onwards
													tot += histos.at(number_layer)->GetBinContent(bin) * (histos.at(number_layer)->GetBinLowEdge(bin) - i);//Sum the total number of hits in each of these bins
													deadcells += histos.at(number_layer)->GetBinContent(bin);
									}
									histos_bufferdepth.at(number_layer)->SetBinContent(i+1, tot);
									histos_deadcells.at(number_layer)->SetBinContent(i+1, deadcells);
					}
	}
	for (int i = 0; i <= max_no_hits; ++i){//For each bufferdepth
					int tot = 0;
					int deadcells = 0;
					for (int bin = i+1; bin < All_Layers_histo->GetNbinsX(); ++bin) {//go through the histo from bufferdepth value onwards
									tot += All_Layers_histo->GetBinContent(bin) * (All_Layers_histo->GetBinLowEdge(bin) - i);//Sum the total number of hits in each of these bins
									deadcells += All_Layers_histo->GetBinContent(bin);
					}
					All_Layers_histo_bufferdepth->SetBinContent(i+1, tot);
					All_Layers_histo_deadcells->SetBinContent(i+1, deadcells);
	}

	std::cout<< "---------------" <<std::endl;
	std::cout<< "Total number of hits counted for this histogram: " << tot_no_hits <<std::endl;
	std::cout<< "---------------" <<std::endl;

	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	canvas->SetLogy(1);

	Draw_multiple_plots(histos, canvas, false);
  std::stringstream output;
  output << "output/muon_occupancy_" << subdetectorname;
	canvas->Print((output.str() + ".pdf").c_str());
	canvas->Print((output.str() + ".cxx").c_str());

	Draw_multiple_plots(histos_numcells, canvas, true);
  std::stringstream output2;
  output2 << "output/muon_occupancy_numcells_" << subdetectorname;
	canvas->Print((output2.str() + ".pdf").c_str());
	canvas->Print((output2.str() + ".cxx").c_str());

	Draw_multiple_plots(histos_bufferdepth, canvas, false);
  std::stringstream output3;
  output3 << "output/muon_occupancy_bufferdepth_" << subdetectorname;
	canvas->Print((output3.str() + ".pdf").c_str());
	canvas->Print((output3.str() + ".cxx").c_str());

	Draw_multiple_plots(histos_deadcells, canvas, false);
  std::stringstream output4;
  output4 << "output/muon_occupancy_deadcells_" << subdetectorname;
	canvas->Print((output4.str() + ".pdf").c_str());
	canvas->Print((output4.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_histo,canvas, false);
  std::stringstream All_output;
  All_output << "output/muon_occupancy_all_layers_" << subdetectorname;
	canvas->Print((All_output.str() + ".pdf").c_str());
	canvas->Print((All_output.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_histo_numcells,canvas, true); 
	std::stringstream All_output2;
  All_output2 << "output/muon_occupancy_numcells_all_layers_" << subdetectorname;
	canvas->Print((All_output2.str() + ".pdf").c_str());
	canvas->Print((All_output2.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_histo_bufferdepth,canvas, false);
  std::stringstream All_output3;
  All_output3 << "output/muon_occupancy_bufferdepth_all_layers_" << subdetectorname;
	canvas->Print((All_output3.str() + ".pdf").c_str());
	canvas->Print((All_output3.str() + ".cxx").c_str());

  Draw_single_plots ( All_Layers_histo_deadcells,canvas, false);
  std::stringstream All_output4;
  All_output4 << "output/muon_occupancy_deadcells_all_layers_" << subdetectorname;
	canvas->Print((All_output4.str() + ".pdf").c_str());
	canvas->Print((All_output4.str() + ".cxx").c_str());

	return 0;
}

long long int MakeNewCellID(double const x, double const y, Subdetector const & component){
  int newX = static_cast<int>(x/component.getCellSizeArea()); //Check if Cell Size Area is the same as Cell Dimension
  int newY = static_cast<int>(y/component.getCellSizeArea()); //Check if Cell Size Area is the same as Cell Dimension
  if(x >= 0) ++newX;
  if(y >= 0) ++newY;
  bitset<32> bitY = newY;
  newY = 0;
  for(int i = 0; i < 31; ++i){
    newY += bitY[i]*pow(2,i);
  }
  return newX << 32 | newY;
}
void Draw_single_plots ( TH1D* histo, TCanvas* canvas, bool normalize){
	histo->SetMinimum(0.1);
	if(normalize == true){
					histo->Scale(1.0/histo->GetBinContent(1));
					histo->SetMinimum( pow(10,-12) );
	}
	histo->SetLineColor(2);
	histo->Draw();
	canvas->Update();
	TPaveStats* st =  (TPaveStats*)histo->GetListOfFunctions()->FindObject("stats");
	st->SetX1NDC(0.75); //new x start position
	st->SetX2NDC(0.9); //new x end position
	st->SetY1NDC(0.8); //new y start position
	st->SetY2NDC(0.9); //new y end position

}
void Draw_multiple_plots ( std::vector< TH1D* > histos, TCanvas* canvas, bool normalize){
	std::vector< TPaveStats* > stats;
	int color = 2; // Very first histogram will be drawn with the color 2, then counted up
	int marker = 20; // Very first histogram will be drawn with the marker style 20, then counted up
	float boxsize = 0.0;
	double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
	for (size_t number_histo = 0; number_histo< histos.size(); ++number_histo) {
					histos.at(number_histo)->SetMinimum(0.1);
					if(normalize == true){
									histos.at(number_histo)->Scale(1.0/histos.at(number_histo)->GetBinContent(1));
									histos.at(number_histo)->SetMinimum( pow(10,-12) );
					}
					histos.at(number_histo)->SetMinimum(0.1);
					histos.at(number_histo)->SetMaximum(max);
					histos.at(number_histo)->Sumw2();
					if(number_histo == 0){
									histos.at(number_histo)->SetLineColor(color);
									histos.at(number_histo)->SetMarkerColor(color);
									histos.at(number_histo)->SetMarkerStyle(marker);
									histos.at(number_histo)->Draw();
									canvas->Update();
									TPaveStats* st =  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats");
									if (histos.size() > 5){
													st->SetX1NDC(0.6); //new x start position
													st->SetX2NDC(0.75); //new x end position
									}
									else{
													st->SetX1NDC(0.75); //new x start position
													st->SetX2NDC(0.9); //new x end position
									}
									st->SetY1NDC(0.8); //new y start position
									st->SetY2NDC(0.9); //new y end position
									stats.push_back(st);
									boxsize = stats.at(number_histo)->GetY2NDC() - stats.at(number_histo)->GetY1NDC();
					}
					if(number_histo > 0){
									color++;
									marker++;
									if(color == 5 || color == 10) color += 1; // 5 would be yellow, 10 would be very light gray 
									histos.at(number_histo)->SetLineColor(color);
									histos.at(number_histo)->SetMarkerColor(color);
									histos.at(number_histo)->SetMarkerStyle(marker);
									histos.at(number_histo)->Draw("SAMES");
									canvas->Update();
									stats.push_back(  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats") );
									stats.at(number_histo)->SetTextColor(color);
									if (histos.size() > 5) {
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
									else{
													stats.at(number_histo)->SetX1NDC(0.75); //new x start position
													stats.at(number_histo)->SetX2NDC(0.9); //new x end position
													stats.at(number_histo)->SetY2NDC(stats.at(number_histo-1)->GetY1NDC()); //new y end position
													stats.at(number_histo)->SetY1NDC(stats.at(number_histo)->GetY2NDC()-boxsize); //new y start position

									}
					}
	}
}
