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
void Draw_All_Layer_plots_together ( std::vector< TH1D* > histo, TCanvas* canvas, bool normalize);
void Draw_single_plots ( TH1D* histo, TCanvas* canvas, bool normalize);
void Draw_multiple_plots (int num_layers, int start_layer, std::vector< TH1D* > histos, std::vector < TPaveStats* > &stats, TCanvas* canvas, bool normalize);
void Print_multiple_plots_from_same_vec (int num_layers, std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::string output);

float weight = 0.08846; //The weight is the same for both scenarios, for both, the electron and the positron line, e.g.: 5sp+wall, elec: (4321/10155 * 898.34)/4321

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
        } while (j <= NUMBER_OF_FILES*2);//2 files per scenario
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
  std::vector< CellHits* > CellHits_vec;

		for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES*2; ++file_iterator) {
      if (file_iterator == 0 || file_iterator == NUMBER_OF_FILES) CellHits_vec.emplace_back(new CellHits( &det ));
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
			float HitPosition_x(0.0), HitPosition_y(0.0), HitPosition_z(0.0);
			//double HitPosition_x(0.0), HitPosition_y(0.0), HitPosition_z(0.0);

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
				CellHits_vec.at( file_iterator/NUMBER_OF_FILES )->Check_CellID(combined_cell_id, HitPosition_x, HitPosition_y, HitPosition_z);
			}
			file->Close();
		}

	//Make histogram vectors for storing the histograms
  std::string subdetectorname = det.getName();
  std::vector < std::string > all_name_sp;
  std::vector < std::string > all_name_spwall ;
	all_name_sp.emplace_back( "All_layers_5sp" );
	all_name_sp.emplace_back( "All_layers_wrt_#cells_5sp" );
	all_name_sp.emplace_back( "All_layers_losthits_5sp" );
	all_name_sp.emplace_back( "All_layers_deadcells_5sp" );
	all_name_spwall.emplace_back( "All_layers_5sp_wall" );
	all_name_spwall.emplace_back( "All_layers_wrt_#cells_5sp_wall" );
	all_name_spwall.emplace_back( "All_layers_losthits_5sp_wall" );
	all_name_spwall.emplace_back( "All_layers_deadcells_5sp_wall" );
  std::vector < std::string > title_sp;
  std::vector < std::string > title_spwall ;
	title_sp.emplace_back( "Occupancy for " + subdetectorname + " - 5 spoilers;Number of hits per cell;Number of cells" );
	title_sp.emplace_back( "Occupancy for " + subdetectorname + " wrt to tot # cells - 5 spoilers;Number of hits per cell;Number of cells" );
	title_sp.emplace_back( "Number of hits lost for a given buffer depth for " + subdetectorname + " - 5 spoilers;Assumend buffer depth;Number of hits lost" );
	title_sp.emplace_back( "Number of dead cells for a given buffer depth for " + subdetectorname + " - 5 spoilers;Assumend buffer depth;Number of dead cells" );
	title_spwall.emplace_back( "Occupancy for " + subdetectorname + " - 5 spoilers+wall;Number of hits per cell;Number of cells" );
	title_spwall.emplace_back( "Occupancy for " + subdetectorname + " wrt to tot # cells - 5 spoilers+wall;Number of hits per cell;Number of cells" );
	title_spwall.emplace_back( "Number of hits lost for a given buffer depth for " + subdetectorname + " - 5 spoilers+wall;Assumend buffer depth;Number of hits lost" );
	title_spwall.emplace_back( "Number of dead cells for a given buffer depth for " + subdetectorname + " - 5 spoilers+wall;Assumend buffer depth;Number of dead cells" );
  std::vector< TH1D* > All_Layers_histo;
	std::vector< TH1D* > All_Layers_histo_numcells;
	std::vector< TH1D* > All_Layers_histo_bufferdepth;
	std::vector< TH1D* > All_Layers_histo_deadcells;
	std::vector< TH1D* > histos;
	std::vector< TH1D* > histos_numcells;
	std::vector< TH1D* > histos_bufferdepth;
	std::vector< TH1D* > histos_deadcells;

  //Find out the maximum number of hits per cell and the total number of hits overall
  std::vector< int > tot_no_hits;
	int max_no_hits = 0;
  for (size_t num_hitcount_classes = 0; num_hitcount_classes < CellHits_vec.size(); ++ num_hitcount_classes){
    tot_no_hits.push_back(0);
    for (size_t vecpos = 0; vecpos < CellHits_vec.at(num_hitcount_classes)->Get_HitCount().size(); ++vecpos) {
      if (CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos) > max_no_hits){
        max_no_hits = CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos);
      }
      tot_no_hits.at(num_hitcount_classes) += CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos);
    }
  }
	int xrange = max_no_hits + max_no_hits/10;
  int max_num_layers = det.getNumberOfLayers();

  //Loop through the CellHits_vec vector
  for (size_t num_hitcount_classes = 0; num_hitcount_classes < CellHits_vec.size(); ++ num_hitcount_classes){
    //Make the histograms
    std::vector< std::string > title;
    std::vector< std::string > name;
    if (num_hitcount_classes == 0){
      title = title_sp;
      name = all_name_sp;
    }
    else{
      title = title_spwall;
      name = all_name_spwall;
    }
    All_Layers_histo.emplace_back( new TH1D(name.at(0).c_str(), title.at(0).c_str(), xrange, 0, xrange) );
    All_Layers_histo_numcells.emplace_back( new TH1D(name.at(1).c_str(), title.at(1).c_str(), xrange, 0, xrange) );
    All_Layers_histo_bufferdepth.emplace_back( new TH1D(name.at(2).c_str(), title.at(2).c_str(), xrange, 0, xrange) );
    All_Layers_histo_deadcells.emplace_back( new TH1D(name.at(3).c_str(), title.at(3).c_str(), xrange, 0, xrange) );

    for (int number_layer = 0; number_layer < max_num_layers; ++number_layer) {
      std::stringstream layername, layername2, layername3, layername4;
      if (num_hitcount_classes == 0){
      layername << "Layer_" << number_layer << "_5sp";
      layername2 << "Layer_" << number_layer << "_numcells_5sp";
      layername3 << "Layer_" << number_layer << "_losthits_5sp";
      layername4 << "Layer_" << number_layer << "_deadcells_5sp";
      }
      else{
      layername << "Layer_" << number_layer << "_5sp_wall";
      layername2 << "Layer_" << number_layer << "_numcells_5sp_wall";
      layername3 << "Layer_" << number_layer << "_losthits_5sp_wall";
      layername4 << "Layer_" << number_layer << "_deadcells_5sp_wall";
      }
      histos.emplace_back(new TH1D(layername.str().c_str(), title.at(0).c_str(), xrange, 0, xrange));
      histos_numcells.emplace_back(new TH1D(layername2.str().c_str(), title.at(1).c_str(), xrange, 0, xrange));
      histos_bufferdepth.emplace_back(new TH1D(layername3.str().c_str(), title.at(2).c_str(), xrange, 0, xrange));
      histos_deadcells.emplace_back(new TH1D(layername4.str().c_str(), title.at(3).c_str(), xrange, 0, xrange));
    }
    //Filling the primary histograms with the entries from the CellHits_vec vectors
    for (size_t vecpos = 0; vecpos < CellHits_vec.at(num_hitcount_classes)->Get_HitCount().size(); ++vecpos) {
      if(CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos) > 0){
        //std::cout << "Layer: " << CellHits_vec.at(num_hitcount_classes)->Get_Layer().at(vecpos) << std::endl;
        int current_layer = CellHits_vec.at(num_hitcount_classes)->Get_Layer().at(vecpos);
        if (num_hitcount_classes == 1) current_layer += max_num_layers;//For the second loop of the num_hitcount_classes loop, the plots in the second half of the histogramms vector shall be filled
        //histos.at(current_layer -1 )->Fill(CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos),weight);//-1 for Silicon detectors only, because layer count starts from 1
        histos.at(current_layer)->Fill(CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos),weight);
        All_Layers_histo.at(num_hitcount_classes)->Fill( CellHits_vec.at(num_hitcount_classes)->Get_HitCount().at(vecpos),weight );
      }
    }
    //Filling numcells plots:
    long long int tot_num_cells = 0;
    long long int tot_num_hitcells = 0;
    for (int number_layer = num_hitcount_classes*max_num_layers; number_layer < (num_hitcount_classes+1)*max_num_layers; ++number_layer) {//For the second loop of the num_hitcount_classes loop, the plots in the second half of the histogramms vector shall be filled
      tot_num_cells += det.getNumberOfCells().at(number_layer - num_hitcount_classes*max_num_layers);//The layers in the subdetector class only go to max_num_layers
      for (int bin = 2; bin < histos.at(number_layer)->GetNbinsX(); ++bin) {
        tot_num_hitcells += histos.at(number_layer)->GetBinContent(bin);
        histos_numcells.at(number_layer)->SetBinContent(bin, histos.at(number_layer)->GetBinContent(bin) );
      }
      histos.at( number_layer )->SetBinContent(1, det.getNumberOfCells().at(number_layer - num_hitcount_classes*max_num_layers) - tot_num_hitcells );
      histos_numcells.at( number_layer )->SetBinContent(1, det.getNumberOfCells().at(number_layer - num_hitcount_classes*max_num_layers));
      tot_num_hitcells = 0;
    }
    for (int bin = 2; bin < All_Layers_histo.at(num_hitcount_classes)->GetNbinsX(); ++bin) {
      tot_num_hitcells += All_Layers_histo.at(num_hitcount_classes)->GetBinContent(bin);
      All_Layers_histo_numcells.at(num_hitcount_classes)->SetBinContent(bin, All_Layers_histo.at(num_hitcount_classes)->GetBinContent(bin));
    }
    All_Layers_histo.at(num_hitcount_classes)->SetBinContent(1, tot_num_cells - tot_num_hitcells);
    All_Layers_histo_numcells.at(num_hitcount_classes)->SetBinContent(1, tot_num_cells);

    //Filling bufferdepth plots:
    for (int number_layer = num_hitcount_classes*max_num_layers; number_layer < (num_hitcount_classes+1)*max_num_layers; ++number_layer) {//For the second loop of the num_hitcount_classes loop, the plots in the second half of the histogramms vector shall be filled
      for (int i = 0; i <= max_no_hits; ++i){//For each bufferdepth
        long long int tot = 0;
        long long int deadcells = 0;
        for (int bin = i+1; bin < histos.at(number_layer)->GetNbinsX(); ++bin) {//go through the histo from bufferdepth value onwards
          tot += histos.at(number_layer)->GetBinContent(bin) * (histos.at(number_layer)->GetBinLowEdge(bin) - i);//Sum the total number of hits in each of these bins
          deadcells += histos.at(number_layer)->GetBinContent(bin);
        }
        histos_bufferdepth.at(number_layer)->SetBinContent(i+1, tot);
        histos_deadcells.at(number_layer)->SetBinContent(i+1, deadcells);
      }
    }
    for (int i = 0; i <= max_no_hits; ++i){//For each bufferdepth
      long long int tot = 0;
      long long int deadcells = 0;
      for (int bin = i+1; bin < All_Layers_histo.at(num_hitcount_classes)->GetNbinsX(); ++bin) {//go through the histo from bufferdepth value onwards
        tot += All_Layers_histo.at(num_hitcount_classes)->GetBinContent(bin) * (All_Layers_histo.at(num_hitcount_classes)->GetBinLowEdge(bin) - i);//Sum the total number of hits in each of these bins
        deadcells += All_Layers_histo.at(num_hitcount_classes)->GetBinContent(bin);
      }
      All_Layers_histo_bufferdepth.at(num_hitcount_classes)->SetBinContent(i+1, tot);
      All_Layers_histo_deadcells.at(num_hitcount_classes)->SetBinContent(i+1, deadcells);
    }
  }//end for loop over CellHits_vec entries

  std::cout<< "---------------" <<std::endl;
  std::cout<< "Total number of hits counted for subdetector "<< subdetectorname <<std::endl;
  for (size_t num_hitcount_classes = 0; num_hitcount_classes < CellHits_vec.size(); ++ num_hitcount_classes){
    std::cout<< "for scenario "<< num_hitcount_classes <<": " << tot_no_hits.at(num_hitcount_classes) <<std::endl;
  }
  std::cout<< "---------------" <<std::endl;

	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	canvas->SetLogy(1);

  std::stringstream output;
  output << "output/muon_occupancy_" << subdetectorname;
  Print_multiple_plots_from_same_vec (max_num_layers, histos, canvas, false, output.str());

  std::stringstream output2;
  output2 << "output/muon_occupancy_numcells_" << subdetectorname;
  Print_multiple_plots_from_same_vec (max_num_layers, histos_numcells, canvas, true, output2.str());
  

  std::stringstream output3;
  output3 << "output/muon_occupancy_bufferdepth_" << subdetectorname;
  
  Print_multiple_plots_from_same_vec (max_num_layers, histos_bufferdepth, canvas, false, output3.str());

  std::stringstream output4;
  output4 << "output/muon_occupancy_deadcells_" << subdetectorname;
  Print_multiple_plots_from_same_vec (max_num_layers, histos_deadcells, canvas, true, output4.str());


  Draw_All_Layer_plots_together( All_Layers_histo,canvas, false); 
  std::stringstream All_output;
  All_output << "output/muon_occupancy_all_layers_" << subdetectorname;
  canvas->Print((All_output.str() + ".pdf").c_str());
  canvas->Print((All_output.str() + ".cxx").c_str());
  
  Draw_All_Layer_plots_together ( All_Layers_histo_numcells,canvas, true); 
  std::stringstream All_numcells_output1;
  All_numcells_output1 << "output/muon_occupancy_numcells_all_layers_" << subdetectorname;
  canvas->Print((All_numcells_output1.str() + ".pdf").c_str());
  canvas->Print((All_numcells_output1.str() + ".cxx").c_str());
  
  Draw_All_Layer_plots_together ( All_Layers_histo_bufferdepth,canvas, false);
  std::stringstream All_bufferdepth_output1;
  All_bufferdepth_output1 << "output/muon_occupancy_bufferdepth_all_layers_" << subdetectorname;
  canvas->Print((All_bufferdepth_output1.str() + ".pdf").c_str());
  canvas->Print((All_bufferdepth_output1.str() + ".cxx").c_str());
  
  Draw_All_Layer_plots_together ( All_Layers_histo_deadcells,canvas, true);
  std::stringstream All_deadcells_output1;
  All_deadcells_output1 << "output/muon_occupancy_deadcells_all_layers_" << subdetectorname;
  canvas->Print((All_deadcells_output1.str() + ".pdf").c_str());
  canvas->Print((All_deadcells_output1.str() + ".cxx").c_str());
//  Draw_single_plots ( All_Layers_histo.at(0),canvas, false);
//  std::stringstream All_output;
//  All_output << "output/muon_occupancy_all_layers_" << subdetectorname << "_5spoilers";
//	canvas->Print((All_output.str() + ".pdf").c_str());
//	canvas->Print((All_output.str() + ".cxx").c_str());
//  Draw_single_plots ( All_Layers_histo.at(1),canvas, false);
//  std::stringstream All_output2;
//  All_output2 << "output/muon_occupancy_all_layers_" << subdetectorname << "_5spoilers_wall";
//	canvas->Print((All_output2.str() + ".pdf").c_str());
//	canvas->Print((All_output2.str() + ".cxx").c_str());
//
//  Draw_single_plots ( All_Layers_histo_numcells.at(0),canvas, true); 
//	std::stringstream All_numcells_output1;
//  All_numcells_output1 << "output/muon_occupancy_numcells_all_layers_" << subdetectorname << "_5spoilers";
//	canvas->Print((All_numcells_output1.str() + ".pdf").c_str());
//	canvas->Print((All_numcells_output1.str() + ".cxx").c_str());
//  Draw_single_plots ( All_Layers_histo_numcells.at(1),canvas, true); 
//	std::stringstream All_numcells_output2;
//  All_numcells_output2 << "output/muon_occupancy_numcells_all_layers_" << subdetectorname << "_5spoilers_wall";
//	canvas->Print((All_numcells_output2.str() + ".pdf").c_str());
//	canvas->Print((All_numcells_output2.str() + ".cxx").c_str());
//
//  Draw_single_plots ( All_Layers_histo_bufferdepth.at(0),canvas, false);
//  std::stringstream All_bufferdepth_output1;
//  All_bufferdepth_output1 << "output/muon_occupancy_bufferdepth_all_layers_" << subdetectorname << "_5spoilers";
//	canvas->Print((All_bufferdepth_output1.str() + ".pdf").c_str());
//	canvas->Print((All_bufferdepth_output1.str() + ".cxx").c_str());
//  Draw_single_plots ( All_Layers_histo_bufferdepth.at(1),canvas, false);
//  std::stringstream All_bufferdepth_output2;
//  All_bufferdepth_output2 << "output/muon_occupancy_bufferdepth_all_layers_" << subdetectorname << "_5spoilers_wall";
//	canvas->Print((All_bufferdepth_output2.str() + ".pdf").c_str());
//	canvas->Print((All_bufferdepth_output2.str() + ".cxx").c_str());
//
//  Draw_single_plots ( All_Layers_histo_deadcells.at(0),canvas, true);
//  std::stringstream All_deadcells_output1;
//  All_deadcells_output1 << "output/muon_occupancy_deadcells_all_layers_" << subdetectorname << "_5spoilers";
//	canvas->Print((All_deadcells_output1.str() + ".pdf").c_str());
//	canvas->Print((All_deadcells_output1.str() + ".cxx").c_str());
//  Draw_single_plots ( All_Layers_histo_deadcells.at(1),canvas, true);
//  std::stringstream All_deadcells_output2;
//  All_deadcells_output2 << "output/muon_occupancy_deadcells_all_layers_" << subdetectorname << "_5spoilers_wall";
//	canvas->Print((All_deadcells_output2.str() + ".pdf").c_str());
//	canvas->Print((All_deadcells_output2.str() + ".cxx").c_str());

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
void Draw_All_Layer_plots_together ( std::vector< TH1D* > histo, TCanvas* canvas, bool normalize){
	std::vector< TPaveStats* > stats;
	float boxsize = 0.0;
	double max=GetMinMaxForMultipleOverlappingHistograms(histo,true).second;
  for (size_t vec_entry = 0; vec_entry < histo.size(); ++vec_entry){
    histo.at(vec_entry)->SetMinimum(max);
    histo.at(vec_entry)->SetMinimum(0.1);
    if(normalize == true){
      histo.at(vec_entry)->Scale(1.0/histo.at(vec_entry)->GetBinContent(1));
      histo.at(vec_entry)->SetMinimum( pow(10,-12) );
    }
		histo.at(vec_entry)->Sumw2(0);
    if(vec_entry == 0){
      histo.at(vec_entry)->SetLineColor(2);
      histo.at(vec_entry)->Draw();
      canvas->Update();
      stats.push_back( (TPaveStats*)histo.at(vec_entry)->GetListOfFunctions()->FindObject("stats") );
      stats.at(vec_entry)->SetTextColor(2);
      stats.at(vec_entry)->SetX1NDC(0.75); //new x start position
      stats.at(vec_entry)->SetX2NDC(0.9); //new x end position
      stats.at(vec_entry)->SetY1NDC(0.8); //new y start position
      stats.at(vec_entry)->SetY2NDC(0.9); //new y end position
			boxsize = stats.at(vec_entry)->GetY2NDC() - stats.at(vec_entry)->GetY1NDC();
    }
    else{
      histo.at(vec_entry)->SetLineColor(4);
      histo.at(vec_entry)->Draw("SAMES");
      canvas->Update();
      stats.push_back( (TPaveStats*)histo.at(vec_entry)->GetListOfFunctions()->FindObject("stats") );
      stats.at(vec_entry)->SetTextColor(4);
      stats.at(vec_entry)->SetX1NDC(0.75); //new x start position
      stats.at(vec_entry)->SetX2NDC(0.9); //new x end position
      stats.at(vec_entry)->SetY2NDC(stats.at(vec_entry-1)->GetY1NDC()); //new y end position
      stats.at(vec_entry)->SetY1NDC(stats.at(vec_entry)->GetY2NDC()-boxsize); //new y start position
    }
  }
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
void Print_multiple_plots_from_same_vec (int num_layers, std::vector< TH1D* > histos, TCanvas* canvas, bool normalize, std::string output){
	std::vector< TPaveStats* > stats;
  int start_layer = 0;
	Draw_multiple_plots(num_layers, start_layer, histos, stats, canvas, normalize);
  std::stringstream output1;
  output1 << output << "_5spoilers";
	canvas->Print((output1.str() + ".pdf").c_str());
	canvas->Print((output1.str() + ".cxx").c_str());
  

  start_layer = num_layers;
	Draw_multiple_plots(num_layers, start_layer, histos, stats, canvas, normalize);
  std::stringstream output2;
  output2 << output << "_5spoilers_wall";
	canvas->Print((output2.str() + ".pdf").c_str());
	canvas->Print((output2.str() + ".cxx").c_str());
}
void Draw_multiple_plots (int num_layers, int start_layer, std::vector< TH1D* > histos, std::vector < TPaveStats* > &stats, TCanvas* canvas, bool normalize){
	float boxsize = 0.0;
  //int stat_count = 0;
	int color = 2; // Very first histogram will be drawn with the color 2, then counted up
	int marker = 20; // Very first histogram will be drawn with the marker style 20, then counted up
	double max=GetMinMaxForMultipleOverlappingHistograms(histos,true).second;
	for (int number_histo = start_layer; number_histo< start_layer+num_layers; ++number_histo) {
					histos.at(number_histo)->Sumw2(1);
					if(normalize == true){
									histos.at(number_histo)->Scale(1.0/histos.at(number_histo)->GetBinContent(1));
									histos.at(number_histo)->SetMinimum( pow(10,-12) );
					}
					else{
									histos.at(number_histo)->SetMaximum(max);
									histos.at(number_histo)->SetMinimum(0.1);
					}
					if(number_histo == start_layer){
									histos.at(number_histo)->SetLineColor(color);
									histos.at(number_histo)->SetMarkerColor(color);
									histos.at(number_histo)->SetMarkerStyle(marker);
									histos.at(number_histo)->Draw("P");
									canvas->Update();
									TPaveStats* st =  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats");
									st->SetTextColor(color);
									if (num_layers-start_layer > 5){
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
                  //stat_count++;
									boxsize = stats.back()->GetY2NDC() - stats.back()->GetY1NDC();
					}
					if(number_histo > start_layer){
									color++;
									marker++;
									if(color == 5 || color == 10) color += 1; // 5 would be yellow, 10 would be very light gray 
									histos.at(number_histo)->SetLineColor(color);
									histos.at(number_histo)->SetMarkerColor(color);
									histos.at(number_histo)->SetMarkerStyle(marker);
									histos.at(number_histo)->Draw("P,SAMES");
									canvas->Update();
									stats.push_back(  (TPaveStats*)histos.at(number_histo)->GetListOfFunctions()->FindObject("stats") );
                  //stat_count++;
									stats.back()->SetTextColor(color);
									if (num_layers > 5) {
													if(number_histo >= 5+start_layer){
																	stats.back()->SetX1NDC(0.75); //new x start position
																	stats.back()->SetX2NDC(0.9); //new x end position
													}
													else {
																	stats.back()->SetX1NDC(0.6); //new x start position
																	stats.back()->SetX2NDC(0.75); //new x end position
													}
													if(number_histo == 5+start_layer) {
																	stats.back()->SetY1NDC(0.8); //new y end position
																	stats.back()->SetY2NDC(0.9); //new y end position
													}
													else {
																	stats.back()->SetY2NDC(stats.at(number_histo-1)->GetY1NDC()); //new y end position
																	stats.back()->SetY1NDC(stats.back()->GetY2NDC()-boxsize); //new y start position
													}
									} 
									else{
													stats.back()->SetX1NDC(0.75); //new x start position
													stats.back()->SetX2NDC(0.9); //new x end position
													stats.back()->SetY2NDC(stats.at(number_histo-1)->GetY1NDC()); //new y end position
													stats.back()->SetY1NDC(stats.back()->GetY2NDC()-boxsize); //new y start position
									}
					}
	}
}
