#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TLegend.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

#include "Helix_class.h"
#include "GeneralFunctions_SiDBkgSim.h"
#include "UsefulFunctions.h"

#include "Style.h"

using namespace std;
void Print_ProjectionY_plots(std::vector< TH1D* > ProjectionY, bool normalize);

int main(int const argc, char const * const * const argv) {
	UsePhDStyle();

	std::vector<std::string> *inputfilenames = new std::vector<std::string>();
	std::string specialname;
	float BField = 0.0;

	int NUMBER_OF_FILES = 0;
	int NUMBER_OF_BUNCHES = 0;
	bool NUMBER_OF_FILES_set = false;
	bool NUMBER_OF_BUNCHES_set = false;
	bool inputfile_set = false;
	bool outputfilename_set = false;
	bool bfield_set = false;

	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-nb")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-n") 
					&& argv[i + 1] != std::string("-b") 
					&& argv[i + 1] != std::string("-o") 
					&& argv[i + 1] != std::string("-i")) {
				NUMBER_OF_BUNCHES = std::stoi(argv[i + 1]);
				std::cout << "Number of bunches = " << NUMBER_OF_BUNCHES << std::endl;
				NUMBER_OF_BUNCHES_set = true;
			} else {
				std::cerr << "You didn't give an argument for the number of bunches!" << std::endl;
			}
		}
    if (argv[i] == std::string("-n")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-nb") 
					&& argv[i + 1] != std::string("-b") 
					&& argv[i + 1] != std::string("-o") 
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
							&& argv[i + j] != std::string("-b")
							&& argv[i + j] != std::string("-nb")
							&& argv[i + j] != std::string("-o")) {
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
		if (argv[i] == std::string("-o")){
			if (argv[i + 1] != NULL 
				&& argv[i + 1] != std::string("-b") 
				&& argv[i + 1] != std::string("-n") 
				&& argv[i + 1] != std::string("-nb") 
				&& argv[i + 1] != std::string("-i")) {
				specialname = argv[i + 1];
				outputfilename_set = true;
			} else {
				std::cerr << "You didn't give an argument for the output filename!" << std::endl;
			}
		}
		if (argv[i] == std::string("-b")){
			if (argv[i + 1] != NULL 
				&& argv[i + 1] != std::string("-n") 
				&& argv[i + 1] != std::string("-nb") 
				&& argv[i + 1] != std::string("-o") 
				&& argv[i + 1] != std::string("-i")) {
				BField = std::stof(argv[i + 1]);
				bfield_set = true;
			} else {
				std::cerr << "You didn't give an argument for the B-field!" << std::endl;
			}
		}
	}
	if (!inputfile_set || !outputfilename_set || !NUMBER_OF_FILES_set || !NUMBER_OF_BUNCHES_set || !bfield_set) {
		std::cerr
			<< "You didn't give the name for the inputfiles, the amount of inputfiles, the name for the outputfile, the number of bunches, or the B-field. Please try again!"
			<< std::endl;
		exit(1);
	}

	TFile* Outputfile = new TFile(("output/Helix_in_beampipe_"+specialname+".root").c_str(),"RECREATE");
	TTree *outputtree = new TTree("Helix_Tracks","Helix_Tracks");

	//Make histogram for storing the information
	float const zmin = -300.0;
	//float const zmin = 0.0;
	float const zmax = 300.0;
	int const zbin = 800;
	//int const zbin = 10000;
	float const ymin = -29;
	float const ymax = 29;
	int const ybin = 200;
	float const xmin = -29;
	float const xmax = 29;
	int const xbin = 200;

	//std::string const title_x = "Pairs spiraling in the magnetic field;z [mm];x [mm];# of particles per (0.375mm x 0.29mm)";
	std::string const title_x = "Pairs spiraling in the magnetic field;z [mm];x [mm];# of particles per (0.75mm x 0.58mm)";
	//std::string const title_y = "Pairs spiraling in the magnetic field;z [mm];y [mm];# of particles per (0.375mm x 0.29mm)";
	std::string const title_y = "Pairs spiraling in the magnetic field;z [mm];y [mm];# of particles per (0.75mm x 0.58mm)";
	TH2D* histo_x = new TH2D("Helix_tracks_xz", title_x.c_str(), zbin,zmin,zmax, xbin, xmin, xmax);
	TH2D* histo_y = new TH2D("Helix_tracks_yz", title_y.c_str(), zbin,zmin,zmax, ybin, ymin, ymax);
	//TTree for outputfile -> store new x, y and z positions of the helixes in there 
	double tree_x(0),tree_y(0),tree_z(0);
	outputtree->Branch("x",&tree_x);
	outputtree->Branch("y",&tree_y);
	outputtree->Branch("z",&tree_z);

	//For counting all particles drawn, and the ones leaving the beam pipe:
	long long int particles_outside_beampipe = 0;
	long long int particles_inside_beampipe = 0;

	//Initializing the Helix class:
	Helix helix(BField);
	//Looping through the root file(s):
	for (int file_iterator = 0; file_iterator < NUMBER_OF_FILES; ++file_iterator) {
    std::cout << inputfilenames->at(file_iterator) << std::endl;
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
		int particleID = 0;
		bool CreatedInSimulation_Status = false;
		//float creationtime = 0.0;

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
		//tree->SetBranchStatus("Particle_ID", 1);
		//tree->SetBranchAddress("Particle_ID", &particleID);
		//tree->SetBranchStatus("Particle_PDG", 1);
		//tree->SetBranchAddress("Particle_PDG", &particleID);
		tree->SetBranchStatus("CreatedInSimulation_Status", 1);
		tree->SetBranchAddress("CreatedInSimulation_Status", &CreatedInSimulation_Status);
		//tree->SetBranchStatus("CreationTime", 1);
		//tree->SetBranchAddress("CreationTime", &creationtime);

		std::vector< double > vertex;
		std::vector< double > momentum;
		double z = zmin;
		double* helix_positions = nullptr;

		float const beampipe_angle1 = 3.266*M_PI/180;
		float const beampipe_angle2 = 5.329*M_PI/180;
		double const tan_beampipe_angle1( tan(beampipe_angle1) );
		double const tan_beampipe_angle2( tan(beampipe_angle2) );

		bool particle_went_outside_beampipe = false;

		//For trying out FillN:
		double* x_array = new double[zbin];
		double* y_array = new double[zbin];
		double* z_array = new double[zbin];
		double* w_array = new double[zbin];

		long long int const entries = tree->GetEntries();
		for (long long int i = 0; i < entries; ++i) {
			tree->GetEntry(i);
			//if (particleID != -11) continue;
			if (momentum_z < 0 /*|| momentum_z>0.1*/) continue;
			if (CreatedInSimulation_Status == 1) continue;
			vertex = { vertex_x, vertex_y, vertex_z };
			momentum = { momentum_x, momentum_y, momentum_z };

			helix.Set_particlevalues(momentum, charge, vertex); // setting the constant values for the current particle in the helix class

			particle_went_outside_beampipe = false; //assume first that every particle will stay inside the beampipe

			//Looping through the histogramm bins in z:
			for (int step = 1; step <= zbin; ++step){
				helix_positions = helix.Get_position(z);
				x_array[step-1]=helix_positions[0]*1000.0;
				y_array[step-1]=helix_positions[1]*1000.0;
				z_array[step-1]=helix_positions[2]*1000.0;
        w_array[step-1]=1.0/(double)NUMBER_OF_BUNCHES;
        if(i == 1 && step==1) std::cout << w_array[step-1] << std::endl;
				//Fill the output TTree:
				tree_x = x_array[step-1];
				tree_y = y_array[step-1];
				tree_z = z_array[step-1];
				outputtree->Fill();
				//Check if particle leaves the beam pipe, and if yes set boolian to true -> after that the if loop should not be accessed again
				if ( particle_went_outside_beampipe == false &&
						((std::abs(tree_x) > 12 && z <= 62.5) || //beam pipe inside vertex barrel: cylinder with 12mm radius, 32.5mm long
						 (std::abs(tree_y) > 12 && z <= 62.5) ||
						 (z > 62.5 && z <= 205 && std::abs(tree_x) > tan_beampipe_angle1*(z-62.5)+12) ||  //beam pipe outside vertex barrel: cone with half angle of beampipe_angle1, length 205-62.5mm, starting at z=62.5mm
						 (z > 62.5 && z <= 205 && std::abs(tree_y) > tan_beampipe_angle1*(z-62.5)+12) ||
						 (z > 205 && z <= zmax && std::abs(tree_x) > tan_beampipe_angle2*(z-205)+20.13)|| //beam pipe outside vertex barrel: cone with half angle of beampipe_angle2, length 205-62.5mm, starting at z=205mm
						 (z > 205 && z <= zmax && std::abs(tree_y) > tan_beampipe_angle2*(z-205)+20.13)
						) ) {
					particle_went_outside_beampipe = true;
				}
				z = step*(zmax-zmin)/zbin +zmin;
				helix_positions = nullptr;
			}
			histo_x->FillN(zbin, z_array, x_array, w_array,1);//number of entries in arrays, array for x, array for y, array for weights (if NULL then weight=1),step size through arrays
			histo_y->FillN(zbin, z_array, y_array, w_array,1);

			if (particle_went_outside_beampipe == true){
				particles_outside_beampipe++;
			}
			else{
				particles_inside_beampipe++;
			}
		}
		file->Close();
	}

  Outputfile->cd();
	std::vector< TH1D* > ProjectionY_xz;
  TH1D* temp1 = (TH1D*)histo_x->ProjectionY("Projection_xz_1Kink",62.5*zbin/300.,62.5*zbin/300.+1);
  TH1D* temp2 = (TH1D*)histo_x->ProjectionY("Projection_xz_2Kink",205.0*zbin/300.,205.0*zbin/300.+1);
  TH1D* temp3 = (TH1D*)histo_x->ProjectionY("Projection_xz_3Kink",295.0*zbin/300.,295.0*zbin/300.+1);
	ProjectionY_xz.push_back( temp1 );
	ProjectionY_xz.push_back( temp2 );
	ProjectionY_xz.push_back( temp3 );
	std::vector< TH1D* > ProjectionY_yz;
  TH1D* temp4 = (TH1D*)histo_y->ProjectionY("Projection_yz_1Kink",62.5*zbin/300.,62.5*zbin/300.+1);
  TH1D* temp5 = (TH1D*)histo_y->ProjectionY("Projection_yz_2Kink",205.0*zbin/300.,205.0*zbin/300.+1);
  TH1D* temp6 = (TH1D*)histo_y->ProjectionY("Projection_yz_3Kink",295.0*zbin/300.,295.0*zbin/300.+1);
	ProjectionY_yz.push_back( temp4 );
	ProjectionY_yz.push_back( temp5 );
	ProjectionY_yz.push_back( temp6 );

	double original_binsize = (double)(xmax - xmin)/(double)xbin;
	std::vector< double > variable_bins_vec;
	variable_bins_vec.push_back(xmin); 
	double temp = xmin;
	while(temp < xmax){
		if(temp < -12.0 || temp > 12.0){ 
			variable_bins_vec.push_back( temp + 4.0*original_binsize ); 
			temp += 4.0*original_binsize;
		}
		else if( (temp>-12.0 && temp < -2.0) || (temp > 2.0 && temp < 12.0) ){ 
			variable_bins_vec.push_back( temp + 2.0*original_binsize ); 
			temp += 2.0*original_binsize;
		}
		else{
			variable_bins_vec.push_back( temp + 1.0*original_binsize ); 
			temp += 1.0*original_binsize;
		}
	}
	double* variable_bins = &variable_bins_vec[0];
	std::vector< TH1D* > Rebinned_ProjectionY_xz;
	std::vector< TH1D* > Rebinned_ProjectionY_yz;
  TH1D* rebinned_temp1 = dynamic_cast<TH1D*>( ProjectionY_xz.at(0)->Rebin(variable_bins_vec.size()-1,"Rebinned_Projection_xz_1Kink",variable_bins) );
  TH1D* rebinned_temp2 = dynamic_cast<TH1D*>( ProjectionY_xz.at(1)->Rebin(variable_bins_vec.size()-1,"Rebinned_Projection_xz_2Kink",variable_bins) );
  TH1D* rebinned_temp3 = dynamic_cast<TH1D*>( ProjectionY_xz.at(2)->Rebin(variable_bins_vec.size()-1,"Rebinned_Projection_xz_3Kink",variable_bins) );
  TH1D* rebinned_temp4 = dynamic_cast<TH1D*>( ProjectionY_yz.at(0)->Rebin(variable_bins_vec.size()-1,"Rebinned_Projection_yz_1Kink",variable_bins) );
  TH1D* rebinned_temp5 = dynamic_cast<TH1D*>( ProjectionY_yz.at(1)->Rebin(variable_bins_vec.size()-1,"Rebinned_Projection_yz_2Kink",variable_bins) );
  TH1D* rebinned_temp6 = dynamic_cast<TH1D*>( ProjectionY_yz.at(2)->Rebin(variable_bins_vec.size()-1,"Rebinned_Projection_yz_3Kink",variable_bins) );
	Rebinned_ProjectionY_xz.push_back( rebinned_temp1 );
	Rebinned_ProjectionY_xz.push_back( rebinned_temp2 );
	Rebinned_ProjectionY_xz.push_back( rebinned_temp3 );
	Rebinned_ProjectionY_yz.push_back( rebinned_temp4 );
	Rebinned_ProjectionY_yz.push_back( rebinned_temp5 );
	Rebinned_ProjectionY_yz.push_back( rebinned_temp6 );
	for(size_t histo_iterator = 0; histo_iterator < Rebinned_ProjectionY_xz.size(); ++histo_iterator){
		Rebinned_ProjectionY_xz.at(histo_iterator)->SetLineColor(((int)histo_iterator+1));//*2);
		Rebinned_ProjectionY_yz.at(histo_iterator)->SetLineColor(((int)histo_iterator+1));//*2);
	}

	std::cout << "-----------------" << std::endl;
	std::cout << "All particles drawn: " << particles_outside_beampipe + particles_inside_beampipe <<  std::endl;
	std::cout << "All particles outside the beam pipe: " << particles_outside_beampipe <<  std::endl;
	std::cout << "Ratio: "  << std::fixed << std::setprecision(3) << ((float)particles_outside_beampipe)/((float)(particles_outside_beampipe + particles_inside_beampipe)) <<  std::endl;
	std::cout << "-----------------" << std::endl;


	//TLine *line = new TLine(0,12,62.5,12);
	//TLine *nline = new TLine(0,-12,62.5,-12);
	//TLine *line2 = new TLine(62.5,12,205,20.13);
	//TLine *nline2 = new TLine(62.5,-12,205,-20.13);
	//TLine *line3 = new TLine(205,20.13,300,28.99);
	//TLine *nline3 = new TLine(205,-20.13,300,-28.99);
	TLine *line = new TLine(-62.5,12,62.5,12);
	TLine *nline = new TLine(-62.5,-12,62.5,-12);
	TLine *line2 = new TLine(62.5,12,205,20.13);
	TLine *nline2 = new TLine(62.5,-12,205,-20.13);
	TLine *line3 = new TLine(-62.5,12,-205,20.13);
	TLine *nline3 = new TLine(-62.5,-12,-205,-20.13);
	TLine *line4 = new TLine(205,20.13,300,28.99);
	TLine *nline4 = new TLine(205,-20.13,300,-28.99);
	TLine *line5 = new TLine(-205,20.13,-300,28.99);
	TLine *nline5 = new TLine(-205,-20.13,-300,-28.99);
	line->SetLineColor(2);
	nline->SetLineColor(2);
	line2->SetLineColor(2);
	nline2->SetLineColor(2);
	line3->SetLineColor(2);
	nline3->SetLineColor(2);
	line4->SetLineColor(2);
	nline4->SetLineColor(2);
	line5->SetLineColor(2);
	nline5->SetLineColor(2);
	
  line->SetLineWidth(6);
	nline->SetLineWidth(6);
	line2->SetLineWidth(6);
	nline2->SetLineWidth(6);
	line3->SetLineWidth(6);
	nline3->SetLineWidth(6);
	line4->SetLineWidth(6);
	nline4->SetLineWidth(6);
	line5->SetLineWidth(6);
	nline5->SetLineWidth(6);


	gStyle->SetOptStat(0);
	//gStyle->SetOptStat(111111);

	//Plot the histogram and save it
	TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
	canvas->cd();

	canvas->SetLogx(0);
	canvas->SetLogy(0);
	canvas->SetLogz();
	//histo_x->SetMinimum(1e-8);
	//histo_y->SetMinimum(1e-8);

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
	line4->Draw();
	nline4->Draw();
	line5->Draw();
	nline5->Draw();

	std::string histoname_x(histo_x->GetName());

	canvas->Print(("output/"+histoname_x+"_"+specialname+".pdf").c_str());
	canvas->Print(("output/"+histoname_x+"_"+specialname+".cxx").c_str());

	histo_y->Draw("colz");
	line->Draw();
	nline->Draw();
	line2->Draw();
	nline2->Draw();
	line3->Draw();
	nline3->Draw();
	line4->Draw();
	nline4->Draw();
	line5->Draw();
	nline5->Draw();

	std::string histoname_y(histo_y->GetName());

	canvas->Print(("output/"+histoname_y+"_"+specialname+".pdf").c_str());
	canvas->Print(("output/"+histoname_y+"_"+specialname+".cxx").c_str());

	canvas->SetLogy(1);
	Print_ProjectionY_plots( Rebinned_ProjectionY_xz, false );
	canvas->Print(("output/ProjectionY_xz_"+specialname+".pdf").c_str());
	Print_ProjectionY_plots( Rebinned_ProjectionY_yz, false );
	canvas->Print(("output/ProjectionY_yz_"+specialname+".pdf").c_str());

	Print_ProjectionY_plots( Rebinned_ProjectionY_xz, true );
	canvas->Print(("output/ProjectionY_xz_normalized_"+specialname+".pdf").c_str());
	Print_ProjectionY_plots( Rebinned_ProjectionY_yz, true );
	canvas->Print(("output/ProjectionY_yz_normalized_"+specialname+".pdf").c_str());

	Outputfile->Write();
	return 0;
}

void Print_ProjectionY_plots(std::vector< TH1D* > ProjectionY, bool normalize){
	for(size_t histo_iterator = 0; histo_iterator < ProjectionY.size(); ++histo_iterator){
		ProjectionY.at(histo_iterator)->GetXaxis()->SetTitleOffset(1);
		if (normalize == false) ProjectionY.at(histo_iterator)->GetYaxis()->SetTitle("# of particles per (0.03mm x 0.58mm)");
		if (normalize == true) NormalizeHistogram( ProjectionY.at(histo_iterator),1.0 );
		if (histo_iterator==0) ProjectionY.at(histo_iterator)->Draw();
		ProjectionY.at(histo_iterator)->Draw("same");
	}

	TLegend* leg = new TLegend(0.6,0.7,0.95,0.9);
	leg->SetTextSize(0.023);
	leg->SetHeader("Projection at different z-positions");
	leg->AddEntry(ProjectionY.at(0),"Beam pipe kink at z=62.5mm","l");
	leg->AddEntry(ProjectionY.at(1),"Beam pipe kink at z=205 mm","l");
	leg->AddEntry(ProjectionY.at(2),"At z=295 mm","l");
	leg->Draw();
}
