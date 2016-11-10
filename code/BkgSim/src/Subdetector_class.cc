#include "Subdetector_class.h"

#include <iostream>

void MCP::Initialize(){
  Name_ = "MCP";
}
void EcalBarrel::Initialize(){
  Name_ = "EcalBarrel";
  LengthLayerBin_ = 6;
  StartLayerBin_ = 13;
  NumberOfLayers_ = 31;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {40., 0., 0.01};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {10, 0, 10};                                           
  ROOTHisto_binning2D = {70, -1500, 1500, 70, -1500, 1500}; //xbins, xlow, xup, ybins, ylow, yup
  ROOTHisto_binning3D = {100, -2000, 2000, 100, -1500, 1500, 100, -1500, 1500}; //zbins, zlow, zup, xbins, xlow, xup, ybins, ylow, yup
  ROOTHisto_time = {500, 0, 5000};
}								
void EcalEndcap::Initialize(){
  Name_ = "EcalEndcap";
  LengthLayerBin_ = 6;
  StartLayerBin_ = 13;
  NumberOfLayers_ = 31;
  TotCellNumber_ = 102400;

  ROOTEnergyHisto_binning = {40., 0., 0.01};
  ROOTHisto_binning_occupancy = {30, 0, 30};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {35, 0, 70};
  ROOTHisto_binning2D = {70, -1500, 1500, 70, -1500, 1500};                    
  ROOTHisto_binning3D = {70, -2000, 2000, 70, -1500, 1500, 70, -1500, 1500};
  ROOTHisto_time = {250, 0, 1000};
}								
void HcalBarrel::Initialize(){
  Name_ = "HcalBarrel";
  LengthLayerBin_ = 6;
  StartLayerBin_ = 13;
  NumberOfLayers_ = 40;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {40., 0., 0.00002};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {6, 0, 6};                                           
  ROOTHisto_binning2D = {70, -2500, 2500, 70, -2500, 2500};                    
  ROOTHisto_binning3D = {100, -3100, 3100, 70, -2500, 2500, 70, -2500, 2500};
  ROOTHisto_time = {500, 0, 5000};
}								
void HcalEndcap::Initialize(){
  Name_ = "HcalEndcap";
  LengthLayerBin_ = 8;
  StartLayerBin_ = 13;
  NumberOfLayers_ = 45;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {250., 0., 0.00045};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {50, 0, 50};                                           
  ROOTHisto_binning2D = {70, -1500, 1500, 70, -1500, 1500};                    
  ROOTHisto_binning3D = {100, -3100, 3100, 70, -1500, 1500, 70, -1500, 1500};
  ROOTHisto_time = {500, 0, 13000};
}								
void MuonBarrel::Initialize(){
  Name_ = "MuonBarrel";
  LengthLayerBin_ = 8;
  StartLayerBin_ = 15;
  NumberOfLayers_ = 11;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {25., 0., 0.00003};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {10, 0, 10};                                           
  ROOTHisto_binning2D = {150, -6300, 6300, 150, -6300, 6300};                    
  ROOTHisto_binning3D = {100, -3100, 3100, 150, -6300, 6300, 150, -6300, 6300};
  ROOTHisto_time = {500, 0, 10000};
}								
void MuonEndcap::Initialize(){
  Name_ = "MuonEndcap";
  LengthLayerBin_ = 8;
  StartLayerBin_ = 15;
  NumberOfLayers_ = 11;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {35., 0., 0.0075};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {50, 0, 15000};                                           
  ROOTHisto_binning2D = {150, -6300, 6300, 150, -6300, 6300};                    
  ROOTHisto_binning3D = {150, -5700, 5700, 150, -6300, 6300, 150, -6300, 6300};
  ROOTHisto_time = {500, 0, 12000};
}								
void BeamCal::Initialize(){
  Name_ = "BeamCal";
  LengthLayerBin_ = 8;
  StartLayerBin_ = 8;
  NumberOfLayers_ = 50;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {45., 0., 2.5};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {150, 0, 2000000};                                           
  ROOTHisto_binning2D = {100, -140, 140, 100, -140, 140};                    
  ROOTHisto_binning3D = {150, -3200, 3200, 100, -140, 140, 100, -140, 140};
  ROOTHisto_time = {500, 0, 15000};
}								
void LumiCal::Initialize(){
  Name_ = "LumiCal";
  LengthLayerBin_ = 8;
  StartLayerBin_ = 8;
  NumberOfLayers_ = 30;
  TotCellNumber_ = 1024;

  ROOTEnergyHisto_binning = {40., 0., 0.045};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {100, 0, 210};                                           
  ROOTHisto_binning2D = {100, -230, 230, 100, -230, 230};                    
  ROOTHisto_binning3D = {150, -1900, 1900, 100, -230, 230, 100, -230, 230};
  ROOTHisto_time = {200, 0, 3500};
}								
void SiVertexEndcap::Initialize(){
  Name_ = "SiVertexEndcap";
  LengthLayerBin_ = 4;
  StartLayerBin_ = 9;
  NumberOfLayers_ = 4;
  TotCellNumber_ = 1024;

  //TO BE ADJUSTED:
  ROOTEnergyHisto_binning = {15., 0., 0.005};
  ROOTHisto_binning_occupancy = {500, 0, 20000};
  ROOTHisto_binning_occupancy_r = {9, 0, 90, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {100, 0, 210};                                           
  ROOTHisto_binning2D = {50, -100, 100, 50, -100, 100};
  ROOTHisto_binning3D = {150, -1900, 1900, 50, -100, 100, 50, -100, 100};
  ROOTHisto_time = {200, 0, 2000};
}
void SiVertexBarrel::Initialize(){
  Name_ = "SiVertexBarrel";
  LengthLayerBin_ = 4;
  StartLayerBin_ = 9;
  NumberOfLayers_ = 5;
  TotCellNumber_ = 1024;

  //TO BE ADJUSTED:
  ROOTEnergyHisto_binning = {40., 0., 0.045};
  ROOTHisto_binning_occupancy = {500, 0, 20000};
  ROOTHisto_binning_occupancy_r = {9, 0, 90, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {100, 0, 210};                                           
  ROOTHisto_binning2D = {50, -100, 100, 50, -100, 100};
  ROOTHisto_binning3D = {150, -1900, 1900, 50, -100, 100, 50, -100, 100};
  ROOTHisto_time = {200, 0, 2000};
}
void SiTrackerEndcap::Initialize(){
  Name_ = "SiTrackerEndcap";
  LengthLayerBin_ = 4;
  StartLayerBin_ = 9;
  NumberOfLayers_ = 4;
  TotCellNumber_ = 1024;

  //TO BE ADJUSTED:
  ROOTEnergyHisto_binning = {40., 0., 0.045};
  ROOTHisto_binning_occupancy = {45, 0, 90};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {100, 0, 210};                                           
  ROOTHisto_binning2D = {150, -1260, 1260, 150, -1260, 1260};
  ROOTHisto_binning3D = {150, -1900, 1900, 150, -1260, 1260, 150, -1260, 1260};
  ROOTHisto_time = {200, 0, 2000};
}
void SiTrackerBarrel::Initialize(){
  Name_ = "SiTrackerBarrel";
  LengthLayerBin_ = 4;
  StartLayerBin_ = 9;
  NumberOfLayers_ = 5;
  TotCellNumber_ = 1024;

  //TO BE ADJUSTED:
  ROOTEnergyHisto_binning = {40., 0., 0.045};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {100, 0, 210};                                           
  ROOTHisto_binning2D = {100, -230, 230, 100, -230, 230};                    
  ROOTHisto_binning3D = {150, -1900, 1900, 100, -230, 230, 100, -230, 230};
  ROOTHisto_time = {200, 0, 3500};
}
void SiTrackerForward::Initialize(){
  Name_ = "SiTrackerForward";
  LengthLayerBin_ = 4;
  StartLayerBin_ = 9;
  NumberOfLayers_ = 3;
  TotCellNumber_ = 1024;

  //TO BE ADJUSTED:
  ROOTEnergyHisto_binning = {40., 0., 0.045};
  ROOTHisto_binning_occupancy = {12, 0, 12};
  ROOTHisto_binning_occupancy_r = {10, 0, 2000, 12, 0, 12};
  ROOTHisto_binning_occupancy_phi = {16, -3.2, 3.2, 12, 0, 12};
  ROOTHisto_binning1D = {100, 0, 210};                                           
  ROOTHisto_binning2D = {100, -230, 230, 100, -230, 230};                    
  ROOTHisto_binning3D = {150, -1900, 1900, 100, -230, 230, 100, -230, 230};
  ROOTHisto_time = {200, 0, 2000};
}
int Subdetector::GetLayer(long long const cellID){
  LayerInfo_ = new LayerCodeInCellID();
  return LayerInfo_->GetLayer(cellID,StartLayerBin_,LengthLayerBin_);
}
int Subdetector::GetLayer(std::string const cellID){
  LayerInfo_ = new LayerCodeInCellID();
  return LayerInfo_->GetLayer(cellID,StartLayerBin_,LengthLayerBin_);
}

Subdetector::Subdetector(){
}
void Subdetector::Initialize(){
}
std::string Subdetector::GetName() const{
  return Name_;
}
int Subdetector::GetTotCellNumber() const {
  return TotCellNumber_;
}
int Subdetector::GetStartLayerBin() const {
  return StartLayerBin_;
}				
int Subdetector::GetLengthLayerBin() const {
  return LengthLayerBin_;
}				
int Subdetector::GetNumberOfLayers() const {
  return NumberOfLayers_;
}	
std::vector<float> Subdetector::GetROOTEnergyHisto_binning() const{
  return ROOTEnergyHisto_binning;
}
std::vector<float> Subdetector::GetROOTHisto_binning_occupancy() const{
  return ROOTHisto_binning_occupancy;
}
std::vector<float> Subdetector::GetROOTHisto_binning_occupancy_r() const{
  return ROOTHisto_binning_occupancy_r;
}
std::vector<float> Subdetector::GetROOTHisto_binning_occupancy_phi() const{
  return ROOTHisto_binning_occupancy_phi;
}
std::vector<float> Subdetector::GetROOTHisto_binning1D() const{
  return ROOTHisto_binning1D;
}				
std::vector<float> Subdetector::GetROOTHisto_binning2D() const{
  return ROOTHisto_binning2D;
}				
std::vector<float> Subdetector::GetROOTHisto_binning3D() const{
  return ROOTHisto_binning3D;
}				
std::vector<float> Subdetector::GetROOTHisto_time() const{
  return ROOTHisto_time;
}

MCP::MCP(){
}
MCP::~MCP(){
}

EcalBarrel::EcalBarrel(){
}
EcalBarrel::~EcalBarrel(){
}

EcalEndcap::EcalEndcap() {
}
EcalEndcap::~EcalEndcap() {
}

HcalBarrel::HcalBarrel(){
}
HcalBarrel::~HcalBarrel(){
}

HcalEndcap::HcalEndcap(){
}
HcalEndcap::~HcalEndcap(){
}

MuonBarrel::MuonBarrel(){
}
MuonBarrel::~MuonBarrel(){
}

MuonEndcap::MuonEndcap(){
}
MuonEndcap::~MuonEndcap(){
}

BeamCal::BeamCal(){
}
BeamCal::~BeamCal(){
}

LumiCal::LumiCal(){
}
LumiCal::~LumiCal(){
}

SiVertexEndcap::SiVertexEndcap(){
}
SiVertexEndcap::~SiVertexEndcap(){
}

SiVertexBarrel::SiVertexBarrel(){
}
SiVertexBarrel::~SiVertexBarrel(){
}

SiTrackerEndcap::SiTrackerEndcap(){
}
SiTrackerEndcap::~SiTrackerEndcap(){
}

SiTrackerBarrel::SiTrackerBarrel(){
}
SiTrackerBarrel::~SiTrackerBarrel(){
}

SiTrackerForward::SiTrackerForward(){
}
SiTrackerForward::~SiTrackerForward(){
}
