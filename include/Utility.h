#ifndef UTILITY_H
#define UTILITY_H

#include "AudioFile.h"
#include <complex>
#include <string>
#include <valarray>
#include <vector>

namespace POID_DGMK
{

enum class WindowFunctionType
{
  Rectangle,
  Hamming,
  Hanning
};

class Utility
{
public:
  static std::vector<std::vector<double>> GetSegmentedSignal(const std::vector<double>& aSamples,
                                                             int aWindowSize,
                                                             int aHopSize);

  static bool LoadSound(std::string& aFileName, AudioFile<double>& aSoundToUpdate);

  static void LoadSoundUntilSuccessful(std::string& aFileName, AudioFile<double>& aSoundToUpdate);

  static void SaveSignalToFile(const std::string aFileName, const std::vector<double>& aData);

  static void ApplyWindowFunction(std::vector<double>& aData,
                                  WindowFunctionType aWindowFunctionType);

  using TComplex = std::complex<double>;
  using CArray = std::valarray<TComplex>;
  using TComplexRepresentation = std::vector<CArray>;

  static void FFT(CArray& x);
  static void IFFT(CArray& x);
};

} // namespace POID_DGMK
#endif // UTILITY_H