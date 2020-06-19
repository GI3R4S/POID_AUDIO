#ifndef WINDOW_APPLICATOR_H
#define WINDOW_APPLICATOR_H
#include "matplotlibcpp.h"
#include <complex>

namespace POID_DGMK
{

enum class WindowFunctionType
{
  Rectangle,
  Hamming,
  Hann
};

class WindowApplicator
{
public:
  using TComplex = std::complex<double>;
  using CArray = std::valarray<TComplex>;
  using TComplexRepresentation = std::vector<CArray>;

  WindowApplicator(int aWindowSize, int aHopSize, WindowFunctionType aWindowFunctionType);

  std::vector<double> MergeDoubleSegments(const std::vector<std::vector<double>>& aSegments,
                                          int aMaxSize);
  std::vector<double> MergeComplexSegments(const std::vector<CArray>& aSegments, int aMaxSize);

  void ApplyWindowFunction(std::vector<double>& aData);
  static void ApplyWindowFunction(std::vector<double>& aData,
                                  WindowFunctionType aWindowFunctionType);

private:
  int mWindowSize;
  int mHopSize;
  double mFactor = 1;
  WindowFunctionType mWindowFunctionType;
};

} // namespace POID_DGMK
#endif // WINDOW_APPLICATOR_H