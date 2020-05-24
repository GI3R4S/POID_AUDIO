#include "UserChoiceHandler.h"
#include "matplotlibcpp.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <stdio.h>

#include <algorithm>
#include <assert.h>
#include <complex>
#include <valarray>
#include <vector>

namespace {

} // namespace

namespace POID_DGMK {

void UserChoiceHandler::ViewMenu() {
  std::cout << "==============================================\n"
            << "[1] Detect pitch using autocorelation\n"
            << "[2] Detect pitch using cepstral analysis\n"
            << "[L] Load another file\n"
            << "[S] Save modified sound\n"
            << "[V] Display sample info\n"
            << "[Q] Exit program\n"
            << "==============================================\n";
}

void UserChoiceHandler::HandleUserChoice(std::string &aInputFile,
                                         AudioFile<double> &aBaseSound) {
  char userChoice;
  std::cout << "\n Select action: \n";
  std::cin >> userChoice;
  switch (userChoice) {
  case '1': {
    auto &audioSource = aBaseSound;
    int windowSize = 4096;
    auto batches =
        Utility::GetSegmentedSamples(audioSource.samples[0], windowSize);
    std::vector<double> pitches;

    for (const auto &batch : batches) {

      std::vector<double> transformedSignal;

      for (int m = 1; m < batch.size(); ++m) {
        double sum = 0;
        for (int i = 0; i < batch.size(); ++i) {
          double value = batch[i];

          value *= m + i < batch.size() ? batch[m + i] : 0;
          sum += value;
        }
        transformedSignal.push_back(sum);
      }

      // Utility::ShowPlot(batch, "Discrete time signal");
      // Utility::ShowPlot(transformedSignal, "Autocorelation of discrete
      // signal");
      std::vector<double> reducedVector;
      bool isIncreaseReached = false;
      int i = 0;
      for (; i < transformedSignal.size(); ++i) {
        if (transformedSignal[i + 1] > transformedSignal[i]) {
          isIncreaseReached = true;
        }

        if (isIncreaseReached) {
          reducedVector.push_back(transformedSignal[i]);
        }
      }

      auto baseFrequency = 0;

      if (!reducedVector.empty()) {
        auto maxElement =
            std::max_element(reducedVector.begin(), reducedVector.end());
        auto maxIndex = maxElement - reducedVector.begin();
        auto maxVal = *maxElement;

        const auto sampleRate = audioSource.getSampleRate();
        baseFrequency =
            static_cast<double>(sampleRate) /
            ((transformedSignal.size() - reducedVector.size()) + maxIndex + 1);

      } else {
        baseFrequency = 0;
      }

      // std::cout << "Max val: " << maxVal << std::endl;
      // std::cout << "Max index: " << maxIndex << std::endl;
      // std::cout << "SampleRate: " << audioSource.getSampleRate() <<
      // std::endl; std::cout << "Base frequency: " << baseFrequency <<
      // std::endl; ShowPlot(batch, "base"); ShowPlot(transformedSignal,
      // "auto");

      pitches.push_back(baseFrequency);
    }

    std::vector<double> pitchData =
        Utility::GetBaseFreqPlotData(audioSource, pitches, windowSize);

    assert(audioSource.samples[0].size() == pitchData.size());

    Utility::ShowPlot(audioSource.samples[0], "Discrete time domain signal");
    Utility::ShowPlot(pitchData, "Pitch detected using autocorrelation");

    break;
  }
  case '2': {
    auto &audioSource = aBaseSound;

    int windowSize = 4096;
    auto batches =
        Utility::GetSegmentedSamples(audioSource.samples[0], windowSize);
    std::vector<double> pitches;

    for (auto &batch : batches) {

      CArray freqDomainData(windowSize);

      // Utility::ShowPlot(batch, "Basic signal before windowing");
      // okienkowanie przez funkcje Hamminga
      for (int i = 0; i < batch.size(); ++i) {
        batch[i] *= (0.53836 -
                     0.46164 * cos(2 * M_PI * i / (1.0 * (batch.size() - 1))));
      }

      // Utility::ShowPlot(batch, "Basic signal windowed by Hamming function");

      for (int i = 0; i < batch.size(); ++i) {
        freqDomainData[i] = batch[i];
      }

      Utility::FFT(freqDomainData);

      CArray halfOfFreqData(windowSize / 2);

      for (int i = 0; i < windowSize / 2; ++i) {
        halfOfFreqData[i] = freqDomainData[i];
      }

      CArray firstSpectrum(windowSize / 2);
      std::vector<double> firstSpectrumVec(windowSize / 2);

      for (int i = 0; i < windowSize / 2; ++i) {
        double val = log(sqrt(pow(halfOfFreqData[i].real(), 2) +
                              pow(halfOfFreqData[i].imag(), 2)) +
                         M_E);
        firstSpectrum[i] = val;
        firstSpectrumVec[i] = val;
      }

      // ShowPlot(firstSpectrumVec, "Freq domain - log(mag) signal");

      Utility::FFT(firstSpectrum);

      CArray halfOfFirstSpectrum(windowSize / 4);
      for (int i = 0; i < windowSize / 4; ++i) {
        halfOfFirstSpectrum[i] = firstSpectrum[i];
      }

      std::vector<double> secondSpectrum(windowSize / 4);
      for (int i = 0; i < windowSize / 4; ++i) {
        double val = 0;
        if (i >= 10) {
          val = sqrt(pow(halfOfFirstSpectrum[i].real(), 2) +
                     pow(halfOfFirstSpectrum[i].imag(), 2));
        }
        secondSpectrum[i] = val;
      }

      // Utility::ShowPlot(secondSpectrum, "Bla");
      auto maxElement =
          std::max_element(secondSpectrum.begin(), secondSpectrum.end());
      auto maxIndex = maxElement - secondSpectrum.begin();

      const auto sampleRate = audioSource.getSampleRate();

      double baseFrequency = static_cast<double>(sampleRate) / (maxIndex * 2);

      pitches.push_back(baseFrequency);
    }

    std::vector<double> pitchData =
        Utility::GetBaseFreqPlotData(audioSource, pitches, windowSize);

    assert(audioSource.samples[0].size() == pitchData.size());

    Utility::ShowPlot(audioSource.samples[0], "Discrete time domain signal");
    Utility::ShowPlot(pitchData, "Pitch detected using cepstrum");

    break;
  }
  case 'l':
  case 'L': {

    POID_DGMK::Utility::LoadSoundUntilSuccessful(aInputFile, aBaseSound);
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
    break;
  }
  case 'v':
  case 'V': {
    auto &audioSource = aBaseSound;
    audioSource.printSummary();

    std::cin.ignore();
    std::cin.get();
    break;
  }
  default: {
    std::cout << "Unhandled option chosen\n";
    break;
  }
  }
}

} // namespace POID_DGMK