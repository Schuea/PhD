#include <iostream>
#include <sstream>
#include <vector>

#include <TSystem.h>
#include <TROOT.h>
#include <TGaxis.h>
#include <TPaletteAxis.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TList.h>
#include <TH3.h>
#include <TF1.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TImage.h>
#include <TPaveLabel.h>

//set maximum of colour scale
double Max = 310.;
//Studied subdetector
std::string Subdetector = "SiVertexBarrel";
std::string specialname = "ILC250_DanielJeans"; //Change this name!

namespace
{
  //_______________________________________________________________________
  TH3D * clone_sizes(const TH3D * src)
  {
    TString name(src->GetName() + TString("_clone"));
    //
    const double xMin = src->GetXaxis()->GetXmin();
    const double xMax = src->GetXaxis()->GetXmax();
    const int xBins = src->GetNbinsX();
    //
    const double yMin = src->GetYaxis()->GetXmin();
    const double yMax = src->GetYaxis()->GetXmax();
    const int yBins = src->GetNbinsY();
    //
    const double zMin = src->GetZaxis()->GetXmin();
    const double zMax = src->GetZaxis()->GetXmax();
    const int zBins = src->GetNbinsZ();

    return new TH3D(name.Data(), name.Data(), xBins, xMin, xMax, yBins, yMin, yMax, zBins, zMin, zMax);
  }

  //_______________________________________________________________________
  TH3D * transform_hist(const TH3D * src, double (*map)(double))
  {
    TH3D * dst = clone_sizes(src);
    bool Bin_with_MaxValue_set = 0;

    for(int i = 1; i <= src->GetNbinsX(); ++i)
    {
      for(int j = 1; j <= src->GetNbinsY(); ++j)
      {
        for(int k = 1; k <= src->GetNbinsZ(); ++k)
        {
          const double val = src->GetBinContent(i, j, k);
          dst->SetBinContent(i, j, k, map(val));

          //set ONE bin of the original histogram to the maximum value, then make it invisible in the transfer function
          if(val < 1 && Bin_with_MaxValue_set == 0){
            dst->SetBinContent(i, j, k, Max+1.);
            Bin_with_MaxValue_set = 1;			
          }				
        }
      }
    }

    return dst;
  }

  //_______________________________________________________________________
  double log_scale(double val)
  {
    // if(val <= 0.) val = 10e-50;
    // return TMath::Log(val);
    return val;
  }

  //_______________________________________________________________________
  double transfer_function(const double * x, const Double_t *)
  {
    //make empty bins invisible
    if(x){
      if(*x > Max) return 0.;
      else if(*x > 3./4.*Max && *x < Max) return 1.;
      else if(*x > Max/2. && *x < 3./4.*Max) return 0.8;
      else if(*x > Max/4. && *x < Max/2.) return 0.4 ;
      else if(*x > 1. && *x < Max/4.) return 0.2 ;
      else return 0. ;
    }
    else return 0.;
  }

}

//_______________________________________________________________________
void th3()
{				
  //gROOT->SetBatch(kTRUE);
  gStyle->SetCanvasPreferGL(1);

  TFile* input = TFile::Open(("../output/OccupancyMap_"+Subdetector+"_"+specialname+".root").c_str(), "READ");
  TFile* output= new TFile(("../output/OccupancyMap_3D_"+Subdetector+"_"+specialname+".root").c_str(), "RECREATE");
  std::vector<TCanvas*> Canvas_;
  std::vector<TH3D*> Hits_3D_;
  std::vector<TH3D*> Hits_3D_Modified;
  std::vector<TList*> lf;

  TImage *img = 0;
  gSystem->Unlink(("../output/gl_3D_gifplot_"+ Subdetector+".gif").c_str());

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(100);
  gStyle->SetCanvasPreferGL(1);

  for(int i = 0; i <= 50; ++i){
    std::stringstream i_ss;
    i_ss << i;

    TH3D * temp = dynamic_cast<TH3D*>(input->Get(("Layer_" + i_ss.str() + "_3D").c_str()));
    Hits_3D_.push_back(temp);
    if(!Hits_3D_.at(i))
    {
      std::cout<<"Could not find the histogram or type is not TH3 derived\n";
      continue;
    }

    TH3D* temp2 = transform_hist(Hits_3D_.at(i), log_scale);
    Hits_3D_Modified.push_back(temp2);
    Hits_3D_Modified.at(i)->GetXaxis()->SetTitle("z (mm)");
    Hits_3D_Modified.at(i)->GetXaxis()->CenterTitle();
    Hits_3D_Modified.at(i)->GetYaxis()->SetTitle("x (mm)");
    Hits_3D_Modified.at(i)->GetYaxis()->CenterTitle();
    Hits_3D_Modified.at(i)->GetZaxis()->SetTitle("y (mm)");
    Hits_3D_Modified.at(i)->GetZaxis()->CenterTitle();



    //Transfer function.
    TList * templf = Hits_3D_Modified.at(i)->GetListOfFunctions();
    lf.push_back(templf);
    if(lf.at(i))
    {
      TF1 * tf = new TF1("TransferFunction", transfer_function);
      lf.at(i)->Add(tf);
    }

    TCanvas* tempcanvas = new TCanvas("glTH3",("gl_3D_plot"+i_ss.str()).c_str(), 200, 10, 600, 600);
    Canvas_.push_back(tempcanvas);

    Canvas_.at(i)->cd();
    Canvas_.at(i)->Update();

    TPaveLabel *title = new TPaveLabel(0.24, 0.85, 0.76, 0.91,("Layer "+i_ss.str()).c_str());
    //title->SetFillColor(32);
    title->Draw();
    Canvas_.at(i)->Update();
    TPad *th3Pad  = new TPad("box", "box", 0.01, 0.01, 0.99, 0.84);
    th3Pad->SetTopMargin(0);
    th3Pad->SetBottomMargin(0);
    th3Pad->Draw();
    Canvas_.at(i)->Update();
    th3Pad->cd();
    gStyle->SetPalette(1);
    Hits_3D_Modified.at(i)->SetContour(Max);
    Hits_3D_Modified.at(i)->SetDirectory(0);
    Hits_3D_Modified.at(i)->Draw("glcolz");
    Hits_3D_Modified.at(i)->Write();

    Canvas_.at(i)->Update();
    Canvas_.at(i)->Write();
    //Canvas_.at(i)->SaveAs(("../output/HitMap_3D_"+specialname+"_"+ Subdetector +"_Layer_"+i_ss.str()+".C").c_str());
    Canvas_.at(i)->SaveAs(("../output/HitMap_3D_"+specialname+"_"+ Subdetector +"_Layer_"+i_ss.str()+".png").c_str());
    //Canvas_.at(i)->Update();

    delete Hits_3D_Modified.at(i);
    delete title;
    delete th3Pad;
    delete Canvas_.at(i);

    delete img; // delete previous image
    // Read image data. Image can be in any format, e.g. png, gif, etc.
    img = TImage::Open(("../output/HitMap_3D_"+specialname+"_"+ Subdetector +"_Layer_"+i_ss.str()+".png").c_str());
    img->WriteImage(("../output/gl_3D_gifplot_"+ Subdetector +".gif+").c_str()); // 60 centiseconds delay
    //img->WriteImage(("../output/gl_3D_gifplot_"+ Subdetector +".gif+60++").c_str()); // 60 centiseconds delay
  }


  output->Close();
  input->Close();
}
