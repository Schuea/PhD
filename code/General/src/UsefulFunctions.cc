#include "UsefulFunctions.h"

#include "TH1.h"
#include "TCanvas.h"

#include <limits>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

float FindMax(float const value, float max){
  if (value > max)  max = value;
  return max;
}
float FindMax(int const value, float max){
  return FindMax(float(value),max);
}
float FindMax(float const value, int max){
  return FindMax(value,float(max));
}
int FindMax(int const value, int max){
  return int(FindMax(float(value),float(max)));
}
float FindMin(float const value, float min){
  if (value < min)  min = value;
  return min;
}
float FindMin(float const value, int min){
  return FindMin(value,float(min));
}
float FindMin(int const value, float min){
  return FindMin(float(value),min);
}
int FindMin(int const value, int min){
  return int(FindMin(float(value), float(min)));
}

std::string Convert_FloatToString (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

void NormalizeHistogram(TH1 * const histo, float const area) {
	if (histo == NULL) {
		std::cerr << "Trying to normalize the histogram " << histo->GetName() << ", but the histo doesn't exist!"
				<< std::endl;
		throw std::exception();
	} else {
		if (histo->Integral() == 0) {
			std::cerr << "Histogram not filled in the x-axis range you specified" << std::endl;
			std::cerr << "Underflow = " << histo->GetBinContent(0) << ", Overflow = "
					<< histo->GetBinContent(histo->GetNbinsX() + 1) << std::endl;
			throw std::exception();
		} else {
			histo->Scale(area / histo->Integral());
		}
	}
}
void PRINT(TCanvas* canvas, std::string title, std::string PDFtitle){
  std::string outputname1 = "output/" + title + ".pdf";
  std::string outputname2 = "output/" + title + ".cxx";
  canvas->Print(outputname1.c_str());
  canvas->Print(outputname2.c_str());
  canvas->Print(("output/" + PDFtitle).c_str());
}

std::pair< double, double> GetMinMaxForMultipleOverlappingHistograms(std::vector< TH1D* > histos, bool const logscale){
  double min(std::numeric_limits< double >::max()),max(std::numeric_limits< double >::min());
  for(int i = 0; i < histos.size(); ++i){
    if(histos.at(i)->GetMinimum() < min) min = histos.at(i)->GetMinimum();
    if(histos.at(i)->GetMaximum() > max) max = histos.at(i)->GetMaximum();
  }
  if(logscale) max *= 10;
  else max *= 1.1;
  if(logscale) min *= 0.1;
  else if(min < 0) min *= 1.1;
  std::pair< double, double > result(min,max);
  return result;
}
std::pair< double, double> GetMinMaxForMultipleOverlappingHistograms(std::vector< TH1F* > histos, bool const logscale){
  double min(std::numeric_limits< double >::max()),max(std::numeric_limits< double >::min());
  for(int i = 0; i < histos.size(); ++i){
    if(histos.at(i)->GetMinimum() < min) min = histos.at(i)->GetMinimum();
    if(histos.at(i)->GetMaximum() > max) max = histos.at(i)->GetMaximum();
  }
  if(logscale) max *= 10;
  else max *= 1.1;
  if(logscale) min *= 0.1;
  else if(min < 0) min *= 1.1;
  std::pair< double, double > result(min,max);
  return result;
}


TTree* Get_TTree(TFile* inputfile, std::string subdetector_name) {
	std::stringstream temp;
	temp << "Tree_" << subdetector_name;
  if(!inputfile)
  {
    throw std::exception();
  }
	TTree* Tree = nullptr;
	inputfile->GetObject(temp.str().c_str(), Tree);
	if (!Tree) {
		throw std::exception();
	}
	return Tree;
}
