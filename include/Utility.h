#ifndef UTILITY_H
#define UTILITY_H

#include "AudioFile.h"
#include <string>
#include <vector>

namespace POID_DGMK {

class Utility {
public:
  // static HistogramData GenerateHistogramData(const Image &aImage);

  // static void DumpChannelData(const std::vector<uint8_t> &aXValues,
  //                             const std::vector<float> &aYValues,
  //                             const std::string &aFileName);

  // static void GenerateFiles(const HistogramData &aHistogramData);

  static bool LoadSound(std::string &aFileName,
                        AudioFile<double> &aSoundToUpdate);

  static void LoadSoundUntilSuccessful(std::string &aFileName,
                                       AudioFile<double> &aSoundToUpdate);

  static AudioFile<double> &SelectSource(AudioFile<double> &aBaseSound,
                                         AudioFile<double> &aModifiedSound,
                                         AudioFile<double> &aWorkSound);

  // static void ViewHistogram();

  static void ViewMenu();

  static void HandleUserChoice(std::string &aInputFile,
                               AudioFile<double> &aBaseSound,
                               AudioFile<double> &aModifiedSound,
                               AudioFile<double> &aWorkSound);
};

} // namespace POID_DGMK
#endif // UTILITY_H