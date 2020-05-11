#include "Utility.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <stdio.h>

namespace POID_DGMK {

// HistogramData Utility::GenerateHistogramData(const Image &aImage) {
//   HistogramData histogramData;
//   const auto &imgData = aImage.GetImageData();

//   std::map<uint8_t, std::size_t> redChannelNumbers;
//   std::map<uint8_t, std::size_t> greenChannelNumbers;
//   std::map<uint8_t, std::size_t> blueChannelNumbers;
//   std::map<uint8_t, std::size_t> averageChannelNumbers;

//   for (uint16_t i = 0; i <= 255; i++) {
//     redChannelNumbers.insert({i, 0});
//     greenChannelNumbers.insert({i, 0});
//     blueChannelNumbers.insert({i, 0});
//     averageChannelNumbers.insert({i, 0});
//   }

//   for (int i = 0; i < imgData.width(); ++i) {
//     for (int j = 0; j < imgData.height(); ++j) {
//       uint8_t redChannelValue = imgData(i, j, 0);
//       uint8_t greenChannelValue = imgData(i, j, 1);
//       uint8_t blueChannelValue = imgData(i, j, 2);
//       uint8_t averageValue = static_cast<uint8_t>(
//           std::round((0.299 * redChannelValue + 0.587 * greenChannelValue +
//                       0.114 * blueChannelValue)));

//       const auto &itRed = redChannelNumbers.find(redChannelValue);
//       const auto &itGreen = greenChannelNumbers.find(greenChannelValue);
//       const auto &itBlue = blueChannelNumbers.find(blueChannelValue);
//       const auto &itAvg = averageChannelNumbers.find(averageValue);

//       ++(redChannelNumbers[redChannelValue]);
//       ++(greenChannelNumbers[greenChannelValue]);
//       ++(blueChannelNumbers[blueChannelValue]);
//       ++(averageChannelNumbers[averageValue]);
//     }

//     const auto numberOfPixels = imgData.width() * imgData.height();
//     for (uint16_t i = 0; i <= 255; ++i) {
//       histogramData.mColorValues.push_back(i);

//       histogramData.mRedChannelOccurrences.push_back(redChannelNumbers[i]);
//       histogramData.mGreenChannelOccurrences.push_back(greenChannelNumbers[i]);
//       histogramData.mBlueChannelOccurrences.push_back(blueChannelNumbers[i]);
//       histogramData.nAverageOccurrences.push_back(averageChannelNumbers[i]);

//       histogramData.mRedChannelProbabilities.push_back(
//           1.0f * redChannelNumbers[i] / numberOfPixels);
//       histogramData.mGreenChannelProbabilities.push_back(
//           1.0f * greenChannelNumbers[i] / numberOfPixels);
//       histogramData.mBlueChannelProbabilities.push_back(
//           1.0f * blueChannelNumbers[i] / numberOfPixels);
//       histogramData.nAverageProbabilities.push_back(
//           1.0f * averageChannelNumbers[i] / numberOfPixels);
//     }
//   }

//   return histogramData;
// }

// void Utility::DumpChannelData(const std::vector<uint8_t> &aXValues,
//                               const std::vector<float> &aYValues,
//                               const std::string &aFileName) {
//   std::ofstream outputFile;
//   std::string outputFilePath;
//   outputFilePath.append(RESOURCES_DIR);
//   outputFilePath.append(aFileName);
//   outputFile.open(outputFilePath, std::ios::out | std::ios::trunc);
//   assert((aXValues.size() == aYValues.size()) && aFileName.c_str());
//   for (size_t i = 0; i < aXValues.size(); ++i) {
//     outputFile << std::to_string(aXValues[i]) << '\t'
//                << std::to_string(aYValues[i]) << '\n';
//   }

//   outputFile.close();
// }

// void Utility::GenerateFiles(const HistogramData &aHistogramData) {
//   DumpChannelData(aHistogramData.mColorValues,
//                   aHistogramData.mRedChannelProbabilities, "redChannel.dat");
//   DumpChannelData(aHistogramData.mColorValues,
//                   aHistogramData.mGreenChannelProbabilities,
//                   "greenChannel.dat");
//   DumpChannelData(aHistogramData.mColorValues,
//                   aHistogramData.mBlueChannelProbabilities, "blueChannel.dat");
//   DumpChannelData(aHistogramData.mColorValues,
//                   aHistogramData.nAverageProbabilities, "averageRGB.dat");

//   float redChannelMax =
//       *std::max_element(aHistogramData.mRedChannelProbabilities.begin(),
//                         aHistogramData.mRedChannelProbabilities.end());
//   float greenChannelMax =
//       *std::max_element(aHistogramData.mGreenChannelProbabilities.begin(),
//                         aHistogramData.mGreenChannelProbabilities.end());
//   float blueChannelMax =
//       *std::max_element(aHistogramData.mBlueChannelProbabilities.begin(),
//                         aHistogramData.mBlueChannelProbabilities.end());
//   float avgMax = *std::max_element(aHistogramData.nAverageProbabilities.begin(),
//                                    aHistogramData.nAverageProbabilities.end());
//   float yMax =
//       std::max({redChannelMax, greenChannelMax, blueChannelMax, avgMax});
//   yMax *= 1.1f;

//   std::string filePath;
//   filePath.append(RESOURCES_DIR);
//   filePath.append("drawHistogramRGB.sh");
//   filePath.append(" ");
//   filePath.append(RESOURCES_DIR);
//   filePath.append(" -e \"yMax='");
//   filePath.append(std::to_string(yMax));
//   filePath.append("'\"");

//   system(filePath.c_str());
// }

bool Utility::LoadSound(std::string &aFileName,
                        AudioFile<double> &aSoundToUpdate) {
  std::cout << "Insert file name from directory 'resources': \n";
  std::cin >> aFileName;

  std::string filePath;
  filePath.append(RESOURCES_DIR);
  filePath.append(aFileName);

  return aSoundToUpdate.load(filePath);
}

void Utility::LoadSoundUntilSuccessful(std::string &aFileName,
                                       AudioFile<double> &aSoundToUpdate) {
  if (!LoadSound(aFileName, aSoundToUpdate)) {
    LoadSoundUntilSuccessful(aFileName, aSoundToUpdate);
  }
}

AudioFile<double> &Utility::SelectSource(AudioFile<double> &aBaseSound,
                                         AudioFile<double> &aModifiedSound,
                                         AudioFile<double> &aWorkSound) {
  char choice;
  do {
    std::cout
        << "Select source picture: '1' - base, '2' - modified, '3' - work\n";
    std::cin >> choice;
  } while (choice != '1' && choice != '2' && choice != '3');

  if (choice == '1') {
    return aBaseSound;
  }
  if (choice == '2') {
    return aModifiedSound;
  }
  if (choice == '3') {
    return aWorkSound;
  }
  return aBaseSound;
}

void Utility::ViewMenu() {
  std::cout << "==============================================\n"
            << "[1] Deduce pitch using autocorelation\n"
            << "[L] Load another file\n"
            << "[R] Reset work sound\n"
            << "[S] Save modified sound\n"
            << "[Q] Exit program\n"
            << "==============================================\n";
}

// void Utility::ViewHistogram() {
//   char input;
//   do {
//     std::cout << "View histogram for picture '1' - base, '2' - modified, '3' "
//                  "- work\n";
//     std::cin >> input;
//   } while (input != '1' && input != '2' && input != '3');

//   Image chosenImage;
// }

void Utility::HandleUserChoice(std::string &aInputFile, AudioFile<double> &aBaseSound,
                               AudioFile<double> &aModifiedSound, AudioFile<double> &aWorkSound) {
  char userChoice;
  std::cout << "\n Select action: \n";
  std::cin >> userChoice;
  switch (userChoice) {
  case '1': {
    break;
  }
  case 'q':
  case 'Q': {
    exit(0);
    break;
  }

  case 's':
  case 'S': {
    std::string fileName;
    std::string fileNameTmp;
    fileName.append(RESOURCES_DIR);
    
    std::cout << "Insert file name: " << std::endl;
    std::cin >> fileNameTmp;
    fileName.append("/");
    fileName.append(fileNameTmp);

    aModifiedSound.save(fileName);
    break;
  }
  case 'u':
  case 'U': {
    aModifiedSound = aWorkSound;
    break;
  }
  case 'r':
  case 'R': {
    aWorkSound = aBaseSound;
    break;
  }
  default: {
    std::cout << "Unhandled option chosen\n";
  }
  }
}
} // namespace POID_DGMK