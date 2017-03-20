#include <iostream>
#include <memory>
#include <string>

#include "Style.h"

#include "TLegend.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

void PrintOptions(){
  cout << "-i or --input [INPUTFILENAME]" << endl;
  exit(0);
}

int main(int const argc, char const * const * const argv){
  UsePhDStyle();
  string input_file_name_1("");
  string input_file_name_2("");

  //Go through the options
  for(int i = 0; i < argc; ++i){
    if(string(argv[i]) == "-h" || string(argv[i]) == "--help"){
      PrintOptions();
    }
    if(string(argv[i]) == "-i1" || string(argv[i]) == "--input1"){
      input_file_name_1 = argv[i+1];
    }
    if(string(argv[i]) == "-i2" || string(argv[i]) == "--input2"){
      input_file_name_2 = argv[i+1];
    }
  }

  string file_stem_1 = input_file_name_1.substr(input_file_name_1.find('_')+1,input_file_name_1.find('.') - input_file_name_1.find('_')-1);
  string file_stem_2 = input_file_name_2.substr(input_file_name_2.find('_')+1,input_file_name_2.find('.') - input_file_name_2.find('_')-1);
  
  //Open the file
  shared_ptr< TFile > file_1(TFile::Open(input_file_name_1.c_str()));
  shared_ptr< TFile > file_2(TFile::Open(input_file_name_2.c_str()));
  
  //Get the TTrees
  shared_ptr< TTree > file1_sampler0(static_cast< TTree* >(file_1->Get("AnalysisUser_Sampler0")));
  shared_ptr< TTree > file1_sampler1(static_cast< TTree* >(file_1->Get("AnalysisUser_Sampler1")));
  shared_ptr< TTree > file1_sampler2(static_cast< TTree* >(file_1->Get("AnalysisUser_Sampler2")));

  shared_ptr< TTree > file2_sampler0(static_cast< TTree* >(file_2->Get("AnalysisUser_Sampler0")));
  shared_ptr< TTree > file2_sampler1(static_cast< TTree* >(file_2->Get("AnalysisUser_Sampler1")));
  shared_ptr< TTree > file2_sampler2(static_cast< TTree* >(file_2->Get("AnalysisUser_Sampler2")));
  
  //Vertex Models as they are for the 3 trees
  file1_sampler0->Draw("VertexModel >> file1vertexmodel0");
  shared_ptr< TH1F > file1_vertexhisto0(static_cast< TH1F* >(gPad->GetPrimitive("file1vertexmodel0")));
  file1_sampler1->Draw("VertexModel >> file1vertexmodel1");
  shared_ptr< TH1F > file1_vertexhisto1(static_cast< TH1F* >(gPad->GetPrimitive("file1vertexmodel1")));
  file1_sampler2->Draw("VertexModel >> file1vertexmodel2");
  shared_ptr< TH1F > file1_vertexhisto2(static_cast< TH1F* >(gPad->GetPrimitive("file1vertexmodel2")));

  file2_sampler0->Draw("VertexModel >> file2vertexmodel0");
  shared_ptr< TH1F > file2_vertexhisto0(static_cast< TH1F* >(gPad->GetPrimitive("file2vertexmodel0")));
  file2_sampler1->Draw("VertexModel >> file2vertexmodel1");
  shared_ptr< TH1F > file2_vertexhisto1(static_cast< TH1F* >(gPad->GetPrimitive("file2vertexmodel1")));
  file2_sampler2->Draw("VertexModel >> file2vertexmodel2");
  shared_ptr< TH1F > file2_vertexhisto2(static_cast< TH1F* >(gPad->GetPrimitive("file2vertexmodel2")));

  //Draw everything nicely
  TCanvas c("vertex_model_bar_chart","A bar chart of the vertex models for the three samplers", 800, 600);
  c.SetLogy();

  TLegend leg(0.45,0.6,0.75,0.8);
  leg.SetTextSize(0.03);
  leg.AddEntry(file1_vertexhisto0.get(),"Vertical Collimator","l");
  leg.AddEntry(file1_vertexhisto1.get(),"RHUL Detector 2","l");
  leg.AddEntry(file1_vertexhisto2.get(),"IP","l");

  file1_vertexhisto0->LabelsOption("v");
  file1_vertexhisto0->SetLineColor(kRed);
  file1_vertexhisto0->SetStats(kFALSE);
  file1_vertexhisto0->GetXaxis()->SetLabelSize(0.025);
  
  file1_vertexhisto1->LabelsOption("v");
  file1_vertexhisto1->SetLineColor(kBlue);
  file1_vertexhisto1->SetStats(kFALSE);
  file1_vertexhisto1->GetXaxis()->SetLabelSize(0.025);
  
  file1_vertexhisto2->LabelsOption("v");
  file1_vertexhisto2->SetLineColor(kGreen + 3);
  file1_vertexhisto2->SetStats(kFALSE);
  file1_vertexhisto2->GetXaxis()->SetLabelSize(0.025);

  file1_vertexhisto0->Draw();
  leg.Draw();
  c.Print((file_stem_1+"vertex_model_sampler_0.cxx").c_str());
  c.Print((file_stem_1+"vertex_model_sampler_0.pdf").c_str());
  file1_vertexhisto1->Draw();
  leg.Draw();
  c.Print((file_stem_1+"vertex_model_sampler_1.cxx").c_str());
  c.Print((file_stem_1+"vertex_model_sampler_1.pdf").c_str());
  file1_vertexhisto2->Draw();
  leg.Draw();
  c.Print((file_stem_1+"vertex_model_sampler_2.cxx").c_str());
  c.Print((file_stem_1+"vertex_model_sampler_2.pdf").c_str());


  file2_vertexhisto0->LabelsOption("v");
  file2_vertexhisto0->SetLineColor(kRed);
  file2_vertexhisto0->SetStats(kFALSE);
  file2_vertexhisto0->GetXaxis()->SetLabelSize(0.025);
  
  file2_vertexhisto1->LabelsOption("v");
  file2_vertexhisto1->SetLineColor(kBlue);
  file2_vertexhisto1->SetStats(kFALSE);
  file2_vertexhisto1->GetXaxis()->SetLabelSize(0.025);
  
  file2_vertexhisto2->LabelsOption("v");
  file2_vertexhisto2->SetLineColor(kGreen + 3);
  file2_vertexhisto2->SetStats(kFALSE);
  file2_vertexhisto2->GetXaxis()->SetLabelSize(0.025);

  file2_vertexhisto0->Draw();
  leg.Draw();
  c.Print((file_stem_2+"vertex_model_sampler_0.cxx").c_str());
  c.Print((file_stem_2+"vertex_model_sampler_0.pdf").c_str());
  file2_vertexhisto1->Draw();
  leg.Draw();
  c.Print((file_stem_2+"vertex_model_sampler_1.cxx").c_str());
  c.Print((file_stem_2+"vertex_model_sampler_1.pdf").c_str());
  file2_vertexhisto2->Draw();
  leg.Draw();
  c.Print((file_stem_2+"vertex_model_sampler_2.cxx").c_str());
  c.Print((file_stem_2+"vertex_model_sampler_2.pdf").c_str());

  //Overlay of the vertex models
  vector< string > file1_vertexhisto0_labels;
  for(int i = 1; i < file1_vertexhisto0->GetNbinsX(); ++i){
    file1_vertexhisto0_labels.push_back(file1_vertexhisto0->GetXaxis()->GetBinLabel(i));
  }
  vector< string > file1_vertexhisto1_labels;
  for(int i = 1; i < file1_vertexhisto1->GetNbinsX(); ++i){
    file1_vertexhisto1_labels.push_back(file1_vertexhisto1->GetXaxis()->GetBinLabel(i));
  }
  vector< string > file1_vertexhisto2_labels;
  for(int i = 1; i < file1_vertexhisto2->GetNbinsX(); ++i){
    file1_vertexhisto2_labels.push_back(file1_vertexhisto2->GetXaxis()->GetBinLabel(i));
  }
  
  vector< string > file2_vertexhisto0_labels;
  for(int i = 1; i < file2_vertexhisto0->GetNbinsX(); ++i){
    file2_vertexhisto0_labels.push_back(file2_vertexhisto0->GetXaxis()->GetBinLabel(i));
  }
  vector< string > file2_vertexhisto1_labels;
  for(int i = 1; i < file2_vertexhisto1->GetNbinsX(); ++i){
    file2_vertexhisto1_labels.push_back(file2_vertexhisto1->GetXaxis()->GetBinLabel(i));
  }
  vector< string > file2_vertexhisto2_labels;
  for(int i = 1; i < file2_vertexhisto2->GetNbinsX(); ++i){
    file2_vertexhisto2_labels.push_back(file2_vertexhisto2->GetXaxis()->GetBinLabel(i));
  }

  //Bins start 1 size big, but will grow as we add more labels
  shared_ptr< TH1F > file1_ordered_vertexhisto0(new TH1F("file1_ordered_vertexhisto0","",1,0,1));
  shared_ptr< TH1F > file2_ordered_vertexhisto0(new TH1F("file2_ordered_vertexhisto0","",1,0,1));
  //Go through each bin from file 2
  for(int i = 1; i < file2_vertexhisto0->GetNbinsX(); ++i){
    bool isFile2Alone(true);
    //Go through each bin from file 1
    for(int j = 1; j < file1_vertexhisto0->GetNbinsX(); ++j){
      //Check if the label for this bin is the same in file 1 and file 2
      if(string(file2_vertexhisto0->GetXaxis()->GetBinLabel(i)) != string(file1_vertexhisto0->GetXaxis()->GetBinLabel(j))) continue;
      isFile2Alone = false;
    }
    //If a bin label exists in file 2 that isn't in file 1, we have to create it in file 1 and fill it in file 2
    if(isFile2Alone){
      file1_ordered_vertexhisto0->Fill(file2_vertexhisto0->GetXaxis()->GetBinLabel(i),0); // Fill the other histogram with a zero-weighted event to force creation of the bin
      for(int j = 0; j < file2_vertexhisto0->GetBinContent(i); ++j){
        file2_ordered_vertexhisto0->Fill(file2_vertexhisto0->GetXaxis()->GetBinLabel(i),1);
      }
    }
  }
  for(int i = 1; i < file1_vertexhisto0->GetNbinsX(); ++i){
    bool isFile1Alone(true);
    for(int j = 1; j < file2_vertexhisto0->GetNbinsX(); ++j){
      if(string(file1_vertexhisto0->GetXaxis()->GetBinLabel(i)) != string(file2_vertexhisto0->GetXaxis()->GetBinLabel(j))) continue;
      isFile1Alone = false;
      for(int k = 0; k < file1_vertexhisto0->GetBinContent(i); ++k){
        file1_ordered_vertexhisto0->Fill(file1_vertexhisto0->GetXaxis()->GetBinLabel(i),1);
      }
      for(int k = 0; k < file2_vertexhisto0->GetBinContent(j); ++k){
        file2_ordered_vertexhisto0->Fill(file2_vertexhisto0->GetXaxis()->GetBinLabel(j),1);
      }
    }
    if(isFile1Alone){
      file2_ordered_vertexhisto0->Fill(file1_vertexhisto0->GetXaxis()->GetBinLabel(i),0); // Fill the other histogram with a zero-weighted event to force creation of the bin
      for(int j = 0; j < file1_vertexhisto0->GetBinContent(i); ++j){
        file1_ordered_vertexhisto0->Fill(file1_vertexhisto0->GetXaxis()->GetBinLabel(i),1);
      }
    }
  }
  file1_ordered_vertexhisto0->LabelsDeflate();
  file2_ordered_vertexhisto0->LabelsDeflate();

  file1_ordered_vertexhisto0->LabelsOption("v");
  file1_ordered_vertexhisto0->GetXaxis()->SetLabelSize(0.025);

  file1_ordered_vertexhisto0->SetLineColor(kRed);
  file2_ordered_vertexhisto0->SetLineColor(kBlue);
  file1_ordered_vertexhisto0->SetStats(kFALSE);
  file1_ordered_vertexhisto0->Draw("hist");
  file2_ordered_vertexhisto0->Draw("hist,same");


  TLegend legcomp(0.5,0.6,0.8,0.8);
  legcomp.SetTextSize(0.03);
  string leg1entry("");
  string leg2entry("");
  if(file_stem_1 == "20000_24_12_12_034"){
    leg1entry = "Open Collimator";
    leg2entry = "Closed Collimator";
  } else{
    leg1entry = "Closed Collimator";
    leg2entry = "Open Collimator";
  }
  //legcomp.SetHeader("Vertical Collimator","C");
  legcomp.SetHeader("Vertical Collimator");
  legcomp.AddEntry(file1_ordered_vertexhisto0.get(),leg1entry.c_str(),"l");
  legcomp.AddEntry(file2_ordered_vertexhisto0.get(),leg2entry.c_str(),"l");
  legcomp.Draw();

  c.Print(string(file_stem_1 + "-"+file_stem_2+"_Comparison_Sampler0.pdf").c_str());
  c.Print(string(file_stem_1 + "-"+file_stem_2+"_Comparison_Sampler0.cxx").c_str());
  
  //Bins start 1 size big, but will grow as we add more labels
  shared_ptr< TH1F > file1_ordered_vertexhisto1(new TH1F("file1_ordered_vertexhisto1","",1,0,1));
  shared_ptr< TH1F > file2_ordered_vertexhisto1(new TH1F("file2_ordered_vertexhisto1","",1,0,1));
  //Go through each bin from file 2
  for(int i = 1; i < file2_vertexhisto1->GetNbinsX(); ++i){
    bool isFile2Alone(true);
    //Go through each bin from file 1
    for(int j = 1; j < file1_vertexhisto1->GetNbinsX(); ++j){
      //Check if the label for this bin is the same in file 1 and file 2
      if(string(file2_vertexhisto1->GetXaxis()->GetBinLabel(i)) != string(file1_vertexhisto1->GetXaxis()->GetBinLabel(j))) continue;
      isFile2Alone = false;
    }
    //If a bin label exists in file 2 that isn't in file 1, we have to create it in file 1 and fill it in file 2
    if(isFile2Alone){
      file1_ordered_vertexhisto1->Fill(file2_vertexhisto1->GetXaxis()->GetBinLabel(i),0); // Fill the other histogram with a zero-weighted event to force creation of the bin
      for(int j = 0; j < file2_vertexhisto1->GetBinContent(i); ++j){
        file2_ordered_vertexhisto1->Fill(file2_vertexhisto1->GetXaxis()->GetBinLabel(i),1);
      }
    }
  }
  for(int i = 1; i < file1_vertexhisto1->GetNbinsX(); ++i){
    bool isFile1Alone(true);
    for(int j = 1; j < file2_vertexhisto1->GetNbinsX(); ++j){
      if(string(file1_vertexhisto1->GetXaxis()->GetBinLabel(i)) != string(file2_vertexhisto1->GetXaxis()->GetBinLabel(j))) continue;
      isFile1Alone = false;
      for(int k = 0; k < file1_vertexhisto1->GetBinContent(i); ++k){
        file1_ordered_vertexhisto1->Fill(file1_vertexhisto1->GetXaxis()->GetBinLabel(i),1);
      }
      for(int k = 0; k < file2_vertexhisto1->GetBinContent(j); ++k){
        file2_ordered_vertexhisto1->Fill(file2_vertexhisto1->GetXaxis()->GetBinLabel(j),1);
      }
    }
    if(isFile1Alone){
      file2_ordered_vertexhisto1->Fill(file1_vertexhisto1->GetXaxis()->GetBinLabel(i),0); // Fill the other histogram with a zero-weighted event to force creation of the bin
      for(int j = 0; j < file1_vertexhisto1->GetBinContent(i); ++j){
        file1_ordered_vertexhisto1->Fill(file1_vertexhisto1->GetXaxis()->GetBinLabel(i),1);
      }
    }
  }
  file1_ordered_vertexhisto1->LabelsDeflate();
  file2_ordered_vertexhisto1->LabelsDeflate();

  file1_ordered_vertexhisto1->LabelsOption("v");
  file1_ordered_vertexhisto1->GetXaxis()->SetLabelSize(0.025);

  file1_ordered_vertexhisto1->SetLineColor(kRed);
  file2_ordered_vertexhisto1->SetLineColor(kBlue);
  file1_ordered_vertexhisto1->SetStats(kFALSE);
  file1_ordered_vertexhisto1->Draw("hist");
  file2_ordered_vertexhisto1->Draw("hist,same");

  legcomp.SetHeader("RHUL Detector 2");
  //legcomp.SetHeader("RHUL Detector 2","C");
  legcomp.Draw();

  c.Print(string(file_stem_1 + "-"+file_stem_2+"_Comparison_Sampler1.pdf").c_str());
  c.Print(string(file_stem_1 + "-"+file_stem_2+"_Comparison_Sampler1.cxx").c_str());
  
  //Bins start 1 size big, but will grow as we add more labels
  shared_ptr< TH1F > file1_ordered_vertexhisto2(new TH1F("file1_ordered_vertexhisto2","",1,0,1));
  shared_ptr< TH1F > file2_ordered_vertexhisto2(new TH1F("file2_ordered_vertexhisto2","",1,0,1));
  //Go through each bin from file 2
  for(int i = 1; i < file2_vertexhisto2->GetNbinsX(); ++i){
    bool isFile2Alone(true);
    //Go through each bin from file 1
    for(int j = 1; j < file1_vertexhisto2->GetNbinsX(); ++j){
      //Check if the label for this bin is the same in file 1 and file 2
      if(string(file2_vertexhisto2->GetXaxis()->GetBinLabel(i)) != string(file1_vertexhisto2->GetXaxis()->GetBinLabel(j))) continue;
      isFile2Alone = false;
    }
    //If a bin label exists in file 2 that isn't in file 1, we have to create it in file 1 and fill it in file 2
    if(isFile2Alone){
      file1_ordered_vertexhisto2->Fill(file2_vertexhisto2->GetXaxis()->GetBinLabel(i),0); // Fill the other histogram with a zero-weighted event to force creation of the bin
      for(int j = 0; j < file2_vertexhisto2->GetBinContent(i); ++j){
        file2_ordered_vertexhisto2->Fill(file2_vertexhisto2->GetXaxis()->GetBinLabel(i),1);
      }
    }
  }
  for(int i = 1; i < file1_vertexhisto2->GetNbinsX(); ++i){
    bool isFile1Alone(true);
    for(int j = 1; j < file2_vertexhisto2->GetNbinsX(); ++j){
      if(string(file1_vertexhisto2->GetXaxis()->GetBinLabel(i)) != string(file2_vertexhisto2->GetXaxis()->GetBinLabel(j))) continue;
      isFile1Alone = false;
      for(int k = 0; k < file1_vertexhisto2->GetBinContent(i); ++k){
        file1_ordered_vertexhisto2->Fill(file1_vertexhisto2->GetXaxis()->GetBinLabel(i),1);
      }
      for(int k = 0; k < file2_vertexhisto2->GetBinContent(j); ++k){
        file2_ordered_vertexhisto2->Fill(file2_vertexhisto2->GetXaxis()->GetBinLabel(j),1);
      }
    }
    if(isFile1Alone){
      file2_ordered_vertexhisto2->Fill(file1_vertexhisto2->GetXaxis()->GetBinLabel(i),0); // Fill the other histogram with a zero-weighted event to force creation of the bin
      for(int j = 0; j < file1_vertexhisto2->GetBinContent(i); ++j){
        file1_ordered_vertexhisto2->Fill(file1_vertexhisto2->GetXaxis()->GetBinLabel(i),1);
      }
    }
  }
  file1_ordered_vertexhisto2->LabelsDeflate();
  file2_ordered_vertexhisto2->LabelsDeflate();

  file1_ordered_vertexhisto2->LabelsOption("v");
  file1_ordered_vertexhisto2->GetXaxis()->SetLabelSize(0.025);

  file1_ordered_vertexhisto2->SetLineColor(kRed);
  file2_ordered_vertexhisto2->SetLineColor(kBlue);
  file1_ordered_vertexhisto2->SetStats(kFALSE);
  file1_ordered_vertexhisto2->Draw("hist");
  file2_ordered_vertexhisto2->Draw("hist,same");

  legcomp.SetHeader("IP");
  //legcomp.SetHeader("IP","C");

  legcomp.Draw();

  c.Print(string(file_stem_1 + "-"+file_stem_2+"_Comparison_Sampler2.pdf").c_str());
  c.Print(string(file_stem_1 + "-"+file_stem_2+"_Comparison_Sampler2.cxx").c_str());
  
  //Vertex models for particular processes

}

  //for(int i = 1; i < vertexhisto0->GetNbinsX(); ++i){
  //  cout << vertexhisto0->GetXaxis()->GetBinLabel(i) << endl;
  //}
  //cout << endl;
  //for(int i = 1; i < vertexhisto1->GetNbinsX(); ++i){
  //  cout << vertexhisto1->GetXaxis()->GetBinLabel(i) << endl;
  //}
  //cout << endl;
  //for(int i = 1; i < vertexhisto2->GetNbinsX(); ++i){
  //  cout << vertexhisto2->GetXaxis()->GetBinLabel(i) << endl;
  //}
  //cout << endl;

/*
Attaching file AnalysisUserOutput_Halomodel_20000_16_8_8_008_event.root as _file0...
(TFile *) 0x7f852aa046f0
root [1] .ls
TFile**   AnalysisUserOutput_Halomodel_20000_16_8_8_008_event.root  Output ROOT file from BDSIM AnalysisUser
TFile*   AnalysisUserOutput_Halomodel_20000_16_8_8_008_event.root  Output ROOT file from BDSIM AnalysisUser
KEY: TTree AnalysisUser_Sampler0;1 AnalysisUser TTree containing vertex information of the particles recorded in Sampler0
KEY: TTree AnalysisUser_Sampler1;1 AnalysisUser TTree containing vertex information of the particles recorded in Sampler1
KEY: TTree AnalysisUser_Sampler2;1 AnalysisUser TTree containing vertex information of the particles recorded in Sampler2
root [2] AnalysisUser_Sampler0->Print("all")
 ******************************************************************************
 *Tree    :AnalysisUser_Sampler0: AnalysisUser TTree containing vertex information of the particles recorded in Sampler0 *
 *Entries :    48816 : Total =         2308052 bytes  File  Size =     406173 *
 *        :          : Tree compression factor =   5.70                       *
 ******************************************************************************
 *Br    0 :VertexModel : VertexModel[11]/C                                    *
 *Entries :    48816 : Total  Size=     736256 bytes  File Size  =     101745 *
 *Baskets :       30 : Basket Size=      32000 bytes  Compression=   7.23     *
 *............................................................................*
 *Br    1 :VertexX   : VertexX/F                                              *
 *Entries :    48816 : Total  Size=     196383 bytes  File Size  =      50116 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=   3.91     *
 *............................................................................*
 *Br    2 :VertexY   : VertexY/F                                              *
 *Entries :    48816 : Total  Size=     196383 bytes  File Size  =      52427 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=   3.74     *
 *............................................................................*
 *Br    3 :VertexZ   : VertexZ/F                                              *
 *Entries :    48816 : Total  Size=     196383 bytes  File Size  =      47801 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=   4.10     *
 *............................................................................*
 *Br    4 :VertexProcess : VertexProcess/I                                    *
 *Entries :    48816 : Total  Size=     196449 bytes  File Size  =       8996 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=  21.78     *
 *............................................................................*
 *Br    5 :VertexSubProcess : VertexSubProcess/I                              *
 *Entries :    48816 : Total  Size=     196482 bytes  File Size  =       9580 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=  20.46     *
 *............................................................................*
 *Br    6 :TrackID   : TrackID/I                                              *
 *Entries :    48816 : Total  Size=     196383 bytes  File Size  =      63877 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=   3.07     *
 *............................................................................*
 *Br    7 :PDGID     : PDGID/I                                                *
 *Entries :    48816 : Total  Size=     196361 bytes  File Size  =      12861 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=  15.23     *
 *............................................................................*
 *Br    8 :ParentID  : ParentID/I                                             *
 *Entries :    48816 : Total  Size=     196394 bytes  File Size  =      56858 *
 *Baskets :        7 : Basket Size=      32000 bytes  Compression=   3.45     *
 *............................................................................*
 */
