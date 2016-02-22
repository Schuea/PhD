#ifndef USEFULFUNCTIONS
#define USEFULFUNCTIONS

#include "TH1.h"

#include <algorithm>
#include <utility>
#include <vector>
#include <sstream>

float FindMax(float const value, float max);
float FindMax(int const value, float max);
float FindMax(float const value, int max);
int FindMax(int const value, int max);
float FindMin(float const value, float min);
float FindMin(float const value, int min);
float FindMin(int const value, float min);
int FindMin(int const value, int min);

std::string Convert_FloatToString (float number);
void NormalizeHistogram(TH1* histo, float size = 1.0);

template <class T, class Q>
std::vector <T> operator* (const Q c, std::vector <T> A)
{
    std::transform (A.begin (), A.end (), A.begin (),
                 std::bind1st (std::multiplies <T> () , c)) ;
    return A ;
}

struct CompareFirst
{
  CompareFirst(int val) : val_(val) {}

  bool operator()(const std::pair< int,int >& elem) const {
    return val_ == elem.first;
  }
  private:
    int val_;
};
#endif
