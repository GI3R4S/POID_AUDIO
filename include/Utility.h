#ifndef UTILITY_H
#define UTILITY_H

#include "AudioFile.h"
#include <complex>
#include <string>
#include <valarray>
#include <vector>

namespace POID_DGMK {

class Utility {
public:
  static void ShowPlot(const std::vector<double> &aData,
                       const std::string &aName);

  static std::vector<std::vector<double>>
  GetSegmentedSamples(const std::vector<double> &aSamples, int aWindowSize);

  static bool LoadSound(std::string &aFileName,
                        AudioFile<double> &aSoundToUpdate);

  static void LoadSoundUntilSuccessful(std::string &aFileName,
                                       AudioFile<double> &aSoundToUpdate);

  static std::vector<double>
  GetBaseFreqPlotData(const AudioFile<double> &aAudioSource,
                      const std::vector<double> &aPitches, int aWindowSize);

  using TComplex = std::complex<double>;
  using CArray = std::valarray<TComplex>;
  using TComplexRepresentation = std::vector<CArray>;

  static void FFT(CArray &x);
  static void IFFT(CArray &x);
};

} // namespace POID_DGMK
#endif // UTILITY_H