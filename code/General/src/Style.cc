#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"

#include "Style.h"

void UsePhDStyle ()
{
  static TStyle* myPhDStyle = 0;
  if ( myPhDStyle==0 ) myPhDStyle = PhDStyle();
  gROOT->SetStyle("PhD");
  gROOT->ForceStyle();
}

TStyle* PhDStyle() 
{
  TStyle *myPhDStyle = new TStyle("PhD","PhD style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  myPhDStyle->SetFrameBorderMode(icol);
  myPhDStyle->SetFrameFillColor(icol);
  myPhDStyle->SetCanvasBorderMode(icol);
  myPhDStyle->SetCanvasColor(icol);
  myPhDStyle->SetPadBorderMode(icol);
  myPhDStyle->SetPadColor(icol);
  myPhDStyle->SetStatColor(icol);

  // set the paper & margin sizes
  myPhDStyle->SetPaperSize(20,26);

  // set margin sizes
  myPhDStyle->SetPadTopMargin(0.12);
  myPhDStyle->SetPadRightMargin(0.18);
  myPhDStyle->SetPadBottomMargin(0.11);
  myPhDStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  myPhDStyle->SetTitleOffset(1.3,"xyz"); //The offset has to be set for all axis first, because there is no SetTitleZOffset function
  myPhDStyle->SetTitleXOffset(1.05); //Then set the offset of the titles for the x- and y-axis separetly
  myPhDStyle->SetTitleYOffset(1.4);

  // use large fonts
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  myPhDStyle->SetTextFont(font);

  myPhDStyle->SetTextSize(tsize);
  myPhDStyle->SetLabelFont(font,"x");
  myPhDStyle->SetTitleFont(font,"x");
  myPhDStyle->SetLabelFont(font,"y");
  myPhDStyle->SetTitleFont(font,"y");
  myPhDStyle->SetLabelFont(font,"z");
  myPhDStyle->SetTitleFont(font,"z");
  
  myPhDStyle->SetLabelSize(tsize,"x");
  myPhDStyle->SetTitleSize(tsize,"x");
  myPhDStyle->SetLabelSize(tsize,"y");
  myPhDStyle->SetTitleSize(tsize,"y");
  myPhDStyle->SetLabelSize(tsize,"z");
  myPhDStyle->SetTitleSize(tsize,"z");

	//Force scientific notation with exponent
	TGaxis::SetMaxDigits(4);

  // use bold lines and markers
  myPhDStyle->SetMarkerStyle(20);
  myPhDStyle->SetMarkerSize(1.2);
  myPhDStyle->SetHistLineWidth(2.);
  myPhDStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //myPhDStyle->SetErrorX(0.001);
  // get rid of error bar caps
  myPhDStyle->SetEndErrorSize(0.);

  myPhDStyle->SetOptStat(1111);
  myPhDStyle->SetOptFit(1111);

  // put tick marks on top and RHS of plots
  myPhDStyle->SetPadTickX(1);
  myPhDStyle->SetPadTickY(1);

  myPhDStyle->SetPalette(53); //Black body color scale

  return myPhDStyle;

}
