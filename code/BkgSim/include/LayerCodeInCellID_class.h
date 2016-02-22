#ifndef LAYERCODEINCELLIDCLASS_H_
#define LAYERCODEINCELLIDCLASS_H_

#include <bitset>

class LayerCodeInCellID {

  public:
    LayerCodeInCellID() {}

    virtual ~LayerCodeInCellID() {}

    int GetLayer(std::string const CellID, int const StartBin_layers, int const LengthBin_layers);
    int GetLayer(long long const cellid, int const StartBin_layers, int const LengthBin_layers);

  private:
    //unsigned long long Layer;
    int Layer;

    int FindLayer(std::string const CellID_, int const StartBin_layers, int const LengthBin_layers);
};
#endif /*LAYERCODEINCELLIDCLASS_H_*/
