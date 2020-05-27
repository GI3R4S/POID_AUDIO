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

namespace POID_DGMK {

void UserChoiceHandler::ViewMenu() {
  std::cout << "==============================================\n"
            << "[1] Compare pitch detection between autocorelation and "
               "cepstrum analysis\n"
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

    char windowSizeChoice;
    int windowSize = 0;
    const std::vector<char> kPossibleChoices{'1', '2', '3'};

    do {
      std::cout
          << "\nSelect window size:\n'1' = 1024\n'2' = 2048\n'3' = 4096\n";
      std::cin >> windowSizeChoice;
    } while (std::find(kPossibleChoices.begin(), kPossibleChoices.end(),
                       windowSizeChoice) == kPossibleChoices.end());

    switch (windowSizeChoice) {
    case '1': {
      windowSize = 1024;
      break;
    }
    case '2': {
      windowSize = 2048;
      break;
    }
    case '3': {
      windowSize = 4096;
      break;
    }
    }

    std::cout << "Loaded window size: " << windowSize << std::endl;
    auto plotData =
        Utility::GetSegmentedSamples(audioSource.samples[0], windowSize);
    std::vector<std::vector<double>> batches = plotData.segmentedPlotData;
    std::vector<int> mask = plotData.mask;

    std::vector<double> pitchesAutoCorrelation;

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
      // Utility::ShowPlot(transformedSignal, "Autocorelation of discrete  signal");

      int maxIndexAutocorr = Utility::FindIndexOfMaximum(transformedSignal);
      double baseFrequency =
          maxIndexAutocorr == -1
              ? 0
              : static_cast<double>(audioSource.getSampleRate()) /
                    (maxIndexAutocorr + 1);

      // std::endl; ShowPlot(batch, "base"); ShowPlot(transformedSignal,
      // "auto");

      pitchesAutoCorrelation.push_back(baseFrequency);
    }

    std::vector<double> pitchDataAutocorrelation = Utility::GetBaseFreqPlotData(
        audioSource, pitchesAutoCorrelation, windowSize);

    std::vector<double> pitchesCepstrum;

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
        double val = val = sqrt(pow(halfOfFirstSpectrum[i].real(), 2) +
                                pow(halfOfFirstSpectrum[i].imag(), 2));

        secondSpectrum[i] = val;
      }

      int maxIndexCep = Utility::FindIndexOfMaximum(secondSpectrum);

      double baseFrequencyCep =
          maxIndexCep == -1 ? 0
                            : static_cast<double>(audioSource.getSampleRate()) /
                                  (maxIndexCep * 2);

      pitchesCepstrum.push_back(baseFrequencyCep);
    }

    std::vector<double> pitchDataCepstrum =
        Utility::GetBaseFreqPlotData(audioSource, pitchesCepstrum, windowSize);

    assert(audioSource.samples[0].size() == pitchDataAutocorrelation.size());
    assert(audioSource.samples[0].size() == pitchDataCepstrum.size());

    for (int i = 0; i < audioSource.samples[0].size(); ++i) {
      pitchDataAutocorrelation[i] *= mask[i];
      pitchDataCepstrum[i] *= mask[i];
    }

    // Utility::ShowPlot(audioSource.samples[0], "Discrete time domain signal");

    matplotlibcpp::subplot(2, 1, 1);
    matplotlibcpp::title("Discrete time domain signal");
    matplotlibcpp::named_plot("Signal", audioSource.samples[0]);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(2, 1, 2);
    matplotlibcpp::named_plot("Pitch detected by autocorrelation",
                              pitchDataAutocorrelation);
    matplotlibcpp::named_plot("Pitch detected by cepstrum", pitchDataCepstrum);
    matplotlibcpp::title("Pitch detection comparision");
    matplotlibcpp::legend();
    matplotlibcpp::show();

    auto autocorrPitchData =
        Utility::GeneratePitchSignal(audioSource, pitchDataAutocorrelation);
    auto cepstrumPitchData =
        Utility::GeneratePitchSignal(audioSource, pitchDataCepstrum);

    matplotlibcpp::subplot(3, 1, 1);
    matplotlibcpp::title("Discrete time domain signal");
    matplotlibcpp::named_plot("Signal", audioSource.samples[0]);
    matplotlibcpp::legend();

    matplotlibcpp::subplot(3, 1, 2);
    matplotlibcpp::title("Sine wave generated from autocorrelation pitch");
    matplotlibcpp::named_plot("Signal", autocorrPitchData);
    matplotlibcpp::legend();

    matplotlibcpp::subplot(3, 1, 3);
    matplotlibcpp::title("Sine wave generated from cepstrum pitch");
    matplotlibcpp::named_plot("Signal", cepstrumPitchData);
    matplotlibcpp::legend();

    matplotlibcpp::show();

    AudioFile<double>::AudioBuffer bufferAuto;
    AudioFile<double>::AudioBuffer bufferCep;

    bufferAuto.resize(1);
    bufferCep.resize(1);

    // 3. Set number of samples per channel
    bufferAuto[0].resize(autocorrPitchData.size());
    bufferCep[0].resize(cepstrumPitchData.size());

    bufferAuto[0] = autocorrPitchData;
    bufferCep[0] = cepstrumPitchData;

    AudioFile<double> autoCorrFile;
    AudioFile<double> cepstrumFile;

    autoCorrFile.setAudioBuffer(bufferAuto);
    cepstrumFile.setAudioBuffer(bufferCep);

    autoCorrFile.setBitDepth(24);
    cepstrumFile.setBitDepth(24);
    autoCorrFile.setSampleRate(44100);
    cepstrumFile.setSampleRate(44100);

    autoCorrFile.save("./autoCorrPitch.wav", AudioFileFormat::Wave);
    cepstrumFile.save("./cepstrumPitch.wav", AudioFileFormat::Wave);
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
} // namespace POID_DGMK

} // namespace POID_DGMK