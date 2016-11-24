#ifndef SUBDETECTOR
#define SUBDETECTOR

#include <string>
#include <vector>

class Subdetector{
  private:
    int _numberOfLayers;
    int _startBitLayer;
    int _lengthBitLayer;
    std::vector< double > _rMin;
    std::vector< double > _rMax;
    std::vector< double > _zHalf;
    double _cellSizeX;
    double _cellSizeY;
    double _cellSizeArea;
    double _length;//maybe layer dependent->vector
    std::vector< double > _rLayer;
    std::vector< double > _area;
    std::vector< double > _numberOfCells;
  public:
    Subdetector();
    Subdetector(std::string const subdetector_config_file);

    int getNumberOfLayers() const;
    int getStartBitLayer() const;
    int getLengthBitLayer() const;
    std::vector< double > getRMin() const;
    std::vector< double > getRMax() const;
    std::vector< double > getZHalf() const;
    double getCellSizeX() const;
    double getCellSizeY() const;
    double getCellSizeArea() const;
    double getLength() const;
    std::vector< double > getRLayer() const;
    std::vector< double > getArea() const;
    std::vector< double > getNumberOfCells() const;

    void setNumberOfLayers(int const number_of_layers);
    void setStartBitLayer(int const start_bit_layer);
    void setLengthBitLayer(int const length_bit_layer);
    void setRMin(std::vector< double > const r_min);
    void setRMax(std::vector< double > const r_max);
    void setZHalf(std::vector< double > const z_half);
    void setCellSizeX(double const cell_size_x);
    void setCellSizeY(double const cell_size_y);
    void setCellSizeArea(double const cell_size_area);
    void setLength(double const length);
    void setRLayer(std::vector< double > const r_layer);
    void setArea(std::vector< double > const area);
    void setNumberOfCells(std::vector< double > const number_of_cells);
};

#endif
