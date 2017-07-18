#include "TCanvas.h"
#include "TH1F.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "Style.h"

float Calculate_mean(std::vector<float> vec, std::vector<float> Evec, bool part);
int Calculate_sign(float num);

int main(int const argc, char const * const * const argv){
				UsePhDStyle();

				if(argc > 2){
								std::cerr << "Only give the pairs.dat file!" << std::endl;
								exit(1);
				}
				std::string inputfilename = argv[1]; 

				std::ifstream inputfile;
				inputfile.open(inputfilename);
				if ( ! inputfile ) {
								std::cerr << "Error: Can't open the file.\n";
								exit(1);
				}

				std::string line;
				std::vector< float > Energy;
				std::vector< float > xVertex;
				std::vector< float > yVertex;
				std::vector< float > zVertex;
				std::vector< float > New_zVertex;

				while ( std::getline(inputfile, line) ){
								std::istringstream input(line);
								float E(0.0), vx(0.0), vy(0.0), vz(0.0), x(0.0), y(0.0), z(0.0);

								input >> E >> vx >> vy >> vz >> x >> y >> z;
								Energy.push_back( E );
								xVertex.push_back( x*std::pow(10,-6) );
								yVertex.push_back( y*std::pow(10,-6) );
								zVertex.push_back( z*std::pow(10,-6) );
				}
				float mean_z_ele(0.0), mean_z_posi(0.0), mean_z(0.0);
				bool ele = true;
				bool not_ele = false;
				mean_z_ele = Calculate_mean(zVertex, Energy, ele);
				mean_z_posi = Calculate_mean(zVertex, Energy, not_ele);
				std::cout << mean_z_ele << std::endl;
				std::cout << mean_z_posi << std::endl;

				TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
				TH1F* histo_x = new TH1F("xVertex","Vertex x of the pairs;x [mm];# pairs",50,-0.4,0.4);
				TH1F* histo_y = new TH1F("yVertex","Vertex y of the pairs;y [mm];# pairs",50,-0.4,0.4);
				TH1F* histo_z = new TH1F("zVertex","Vertex z of the pairs;z [mm];# pairs",50,-2.8,2.8);
				for (size_t i = 0; i < zVertex.size(); ++i){
								int sign  = Calculate_sign(Energy.at(i));
								if (sign > 0) mean_z = mean_z_ele;
								else mean_z = mean_z_posi;
								New_zVertex.push_back( zVertex.at(i) + sign*mean_z );
								histo_x->Fill( xVertex.at(i) );
								histo_y->Fill( yVertex.at(i) );
								histo_z->Fill( zVertex.at(i) );
								//histo->Fill( New_zVertex.back() );
				}

				histo_x->Draw();
				canvas->Print("Pairs_xVertex.pdf");
				canvas->Print("Pairs_xVertex.cxx");
				histo_y->Draw();
				canvas->Print("Pairs_yVertex.pdf");
				canvas->Print("Pairs_yVertex.cxx");
				histo_z->Draw();
				canvas->Print("Pairs_zVertex.pdf");
				canvas->Print("Pairs_zVertex.cxx");

				return 0;
}

float Calculate_mean(std::vector<float> vec, std::vector<float> Evec, bool part){
				float mean(0.0);
				int count_part(0);
				for (size_t i = 0; i < vec.size(); ++i){
								if (part == false && Evec.at(i)<0){
												mean += vec.at(i);
												count_part++;
								}
								else if (part == true && Evec.at(i)>0){
												mean += vec.at(i);
												count_part++;
								}
								else continue;
				}
				mean = mean/count_part;
				return mean;
}
int Calculate_sign(float num){
				int sign = num/std::abs( num );
				return sign;
}
