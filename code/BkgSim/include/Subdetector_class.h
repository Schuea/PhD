#ifndef SUBDETECTORCLASS_H_
#define SUBDETECTORCLASS_H_

#include "LayerCodeInCellID_class.h"
#include <vector>

class Subdetector{
  public:
    Subdetector();

    virtual void Initialize();

    std::string GetName() const;
    int GetStartLayerBin() const;
    int GetLengthLayerBin() const;				
    int GetNumberOfLayers() const;	
    int GetLayer(std::string const cellID);
    int GetLayer(long long const cellID);
    int GetTotCellNumber() const;

	std::vector<float> GetROOTEnergyHisto_binning() const;
	std::vector<float> GetROOTHisto_binning_occupancy() const;
	std::vector<float> GetROOTHisto_binning_occupancy_r() const;
	std::vector<float> GetROOTHisto_binning_occupancy_phi() const;
	std::vector<float> GetROOTHisto_binning1D() const;
    std::vector<float> GetROOTHisto_binning2D() const;
    std::vector<float> GetROOTHisto_binning3D() const;
    std::vector<float> GetROOTHisto_time() const;

  protected:
    LayerCodeInCellID *LayerInfo_;
    std::string Name_;
    int LengthLayerBin_;
    int StartLayerBin_;
    int NumberOfLayers_;
    int TotCellNumber_;
    std::vector<float> ROOTEnergyHisto_binning;
    std::vector<float> ROOTHisto_binning_occupancy;
    std::vector<float> ROOTHisto_binning_occupancy_r;
    std::vector<float> ROOTHisto_binning_occupancy_phi;
    std::vector<float> ROOTHisto_binning1D;
    std::vector<float> ROOTHisto_binning2D;
    std::vector<float> ROOTHisto_binning3D;
    std::vector<float> ROOTHisto_time;
};
class MCP : public Subdetector{
  public:
    MCP();
    ~MCP();
    void Initialize();
};
class EcalBarrel : public Subdetector{
  public:
    EcalBarrel();
    ~EcalBarrel();
    void Initialize();
};
class EcalEndcap : public Subdetector{
  public:
    EcalEndcap();
    ~EcalEndcap();
    void Initialize();
};
class HcalBarrel : public Subdetector{
  public:
    HcalBarrel();
    ~HcalBarrel();
    void Initialize();
};
class HcalEndcap : public Subdetector{
  public:
    HcalEndcap();
    ~HcalEndcap();
    void Initialize();
};
class MuonBarrel : public Subdetector{
  public:
    MuonBarrel();
    ~MuonBarrel();
    void Initialize();
};
class MuonEndcap : public Subdetector{
  public:
    MuonEndcap();
    ~MuonEndcap();
    void Initialize();
};
class BeamCal : public Subdetector{
  public:
    BeamCal();
    ~BeamCal();
    void Initialize();
};
class LumiCal : public Subdetector{
  public:
    LumiCal();
    ~LumiCal();
    void Initialize();
};
class SiVertexEndcap : public Subdetector{
  public:
    SiVertexEndcap();
    ~SiVertexEndcap();
    void Initialize();
};
class SiVertexBarrel : public Subdetector{
  public:
    SiVertexBarrel();
    ~SiVertexBarrel();
    void Initialize();
};
class SiTrackerEndcap : public Subdetector{
  public:
    SiTrackerEndcap();
    ~SiTrackerEndcap();
    void Initialize();
};
class SiTrackerBarrel : public Subdetector{
  public:
    SiTrackerBarrel();
    ~SiTrackerBarrel();
    void Initialize();
};
class SiTrackerForward : public Subdetector{
  public:
    SiTrackerForward();
    ~SiTrackerForward();
    void Initialize();
};
#endif /*SUBDETECTORCLASS_H_*/
