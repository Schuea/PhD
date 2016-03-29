#include "UsefulFunctions.h"

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
  canvas->Print(PDFtitle.c_str());
}
