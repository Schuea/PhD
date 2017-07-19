#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <fstream>
#include <sstream>

int main(int const argc, char const * const * const argv){

	std::vector< std::string > inputfilenames;
	std::string outputfilename;

	if (argc < 5) {
		std::cerr << "Please provide the names of the input and output files: " << std::endl;
		std::cerr << "./LumiPlotter -i file1.txt file2.txt -o output" << std::endl;
		exit(2);
	}
	for (int i = 1; i < argc; i++) {
		if (argv[i] == std::string("-i")){
			if (argv[i + 1] != NULL && argv[i + 1] != std::string("-o")) {
				int j = 1;
				do {
					inputfilenames.push_back(argv[i + j]);
					++j;
				} while (argv[i + j] != NULL && argv[i + j] != std::string("-o") );
			} else {
				std::cerr << "You didn't give an argument for the inputfile(s)!" << std::endl;
				exit(2);
			}
		}
		if (argv[i] == std::string("-o")) {
			if (argv[i + 1] != NULL 
					&& argv[i + 1] != std::string("-i")) {
				outputfilename = argv[i + 1];
			} else {
				std::cerr << "You didn't give an argument for the outputfilename!"
					<< std::endl;
				exit(2);
			}
		}
	}

	std::string Buzz_word_START = "$HISTOGRAM::lumi_ee";
	std::string Buzz_word_END = "$HISTOGRAM::lumi_eg";

  std::stringstream new_outputfilename;
  new_outputfilename << "output/" << outputfilename << ".root";
	TFile* ROOTFile = new TFile(new_outputfilename.str().c_str(),"RECREATE","Luminosity histogram");
  std::vector< TH1F* > LumiHistos;

	for(size_t no_files = 0; no_files < inputfilenames.size(); ++ no_files){
	  
    std::vector< float > lumi;
    std::stringstream name;
    name << "Lumi_" << no_files+1;
    TH1F* temp = new TH1F(name.str().c_str(),"Luminosity vs. E_CM",200,0,250.025);
    temp->SetBit(TH1::kIsAverage);
		LumiHistos.push_back( temp );

		std::ifstream inputfile( inputfilenames.at(no_files) );
		std::string line;
	  bool READ = false;
	  bool STOP_READ = false;

		if(inputfile.is_open()){
			std::cout << "Opened the text file "<< inputfilenames.at(no_files) << std::endl;
			while (!inputfile.eof()){
				std::getline(inputfile, line);
				if ( Buzz_word_START.compare(line)==0){
					READ=true;
					std::getline(inputfile, line);
					std::getline(inputfile, line);
				}
				else if ( Buzz_word_END.compare(line)==0){
					READ=false;
					STOP_READ=true;
				}
				if ( READ == true ){	
					std::istringstream in(line);
					std::string col1, col2, col3, col4;
					in >> col1 >> col2 >> col3 >> col4;
          //std::cout << col1 << " - "  << col2 << " - " << col3 << " - " << col4 << std::endl;

					lumi.push_back(std::stof(col1));
					//std::cout << lumi.back() << std::endl;
					lumi.push_back(std::stof(col2));
					//std::cout << lumi.back() << std::endl;
					lumi.push_back(std::stof(col3));
					//std::cout << lumi.back() << std::endl;
					lumi.push_back(std::stof(col4));
					//std::cout << lumi.back() << std::endl;
				}
				if ( STOP_READ == true ) break;
			}
			inputfile.close();
		}
		else{
			std::cout<<"Error! File "<< inputfilenames.at(no_files) << " not found!";
			exit(1);
		}


		for(int bin = 1; bin <= 200; ++bin){
			LumiHistos.at(no_files)->SetBinContent( bin, lumi.at(bin-1) );
		}
	}
	TH1F* combined_LumiHistos = (TH1F*) LumiHistos.at(0)->Clone();
  combined_LumiHistos->SetName("Lumi");
  combined_LumiHistos->SetBit(TH1::kIsAverage);
	for(size_t no_files = 1; no_files < inputfilenames.size(); ++ no_files){
    combined_LumiHistos->Add(LumiHistos.at(no_files));
  }
  combined_LumiHistos->SetMinimum( std::pow(10.0,21.0) );

	std::cout << "Output will be created: " << outputfilename << std::endl;
	ROOTFile->Write();
	ROOTFile->Close();
	return 0;
}
