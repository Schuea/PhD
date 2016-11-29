#ifndef SUBDETECTOR_CLASS_NEW_H
#define SUBDETECTOR_CLASS_NEW_H

#include <string>
#include <vector>

class Subdetector{
  private:
    std::string _name;
    int _numberOfLayers;
    int _startBitLayer;
    int _lengthBitLayer;
    std::vector< double > _rMin;
    std::vector< double > _rMax;
    std::vector< double > _zHalf;
    std::vector< double > _length;
    double _cellSizeX;
    double _cellSizeY;
    double _cellSizeArea;
    std::vector< double > _rLayer;
    std::vector< double > _area;
    std::vector< double > _numberOfCells;
  public:
    Subdetector();
    Subdetector(std::string const subdetector_config_file);

    std::string getName() const;
    int getNumberOfLayers() const;
    int getStartBitLayer() const;
    int getLengthBitLayer() const;
    std::vector< double > getRMin() const;
    std::vector< double > getRMax() const;
    std::vector< double > getZHalf() const;
    double getCellSizeX() const;
    double getCellSizeY() const;
    double getCellSizeArea() const;
    std::vector< double > getLength() const;
    std::vector< double > getRLayer() const;
    std::vector< double > getArea() const;
    std::vector< double > getNumberOfCells() const;

    void setName(std::string const name);
    void setNumberOfLayers(int const number_of_layers);
    void setStartBitLayer(int const start_bit_layer);
    void setLengthBitLayer(int const length_bit_layer);
    void setRMin(std::vector< double > const r_min);
    void setRMax(std::vector< double > const r_max);
    void setZHalf(std::vector< double > const z_half);
    void setCellSizeX(double const cell_size_x);
    void setCellSizeY(double const cell_size_y);
    void setCellSizeArea(double const cell_size_area);
    void setLength(std::vector< double > const length);
    void setRLayer(std::vector< double > const r_layer);
    void setArea(std::vector< double > const area);
    void setNumberOfCells(std::vector< double > const number_of_cells);

    std::vector< double > ConvertCSVToVectorDoubles(std::string const csv);
    void FillUpVector( std::vector< double > & vec);
};

#endif
