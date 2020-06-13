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

namespace POID_DGMK
{

void UserChoiceHandler::ViewMenu()
{
  std::cout << "==============================================\n"
            << "[1] Compare pitch detection between autocorelation and "
               "cepstrum analysis\n"
            << "[L] Load another file\n"
            << "[S] Save modified sound\n"
            << "[V] Display sample info\n"
            << "[Q] Exit program\n"
            << "==============================================\n";
}

void UserChoiceHandler::HandleUserChoice(std::string& aInputFile, AudioFile<double>& aBaseSound)
{
  char userChoice;
  std::cout << "\n Select action: \n";
  std::cin >> userChoice;

  switch (userChoice)
  {
  case '1':
  {
    auto& audioSource = aBaseSound;

    char windowSizeChoice;
    char windowFunctionTypeChoice;
    int windowSize = 2049;
    int hopSize = 1028;
    int filterLength = 399;
    int cutoff_freq = 500;
    WindowFunctionType windowFunctionType;
    const std::vector<char> kPossibleChoices{'1', '2', '3'};

    std::vector<double> signal;
    for (int i = 0; i < 44100; ++i)
    {
      signal.push_back(1);
    }

    auto batchesWithHops =
      Utility::GetSegmentedSignal(audioSource.samples[0], windowSize, hopSize);

    // for (auto& batch : batchesWithHops)
    // {
    //   Utility::ApplyWindowFunction(batch, WindowFunctionType::Hamming);
    // }

    // obliczenie coefficentów
    std::vector<double> coefficients;
    for (int i = 0; i < filterLength; ++i)
    {
      if (i == (filterLength - 1) / 2)
      {
        coefficients.push_back(2.0 * cutoff_freq / 44100.0);
      }
      else
      {
        double value =
          sin((2 * M_PI * cutoff_freq / (44100.0)) * (i - ((filterLength - 1) / 2.0))) /
          (M_PI * (i - ((filterLength - 1) / 2)));
        coefficients.push_back(value);
      }
    }

    // zastosowanie na nich okna
    auto coeffCopy = coefficients;
    Utility::ApplyWindowFunction(coefficients, WindowFunctionType::Hamming);

    // matplotlibcpp::subplot(2, 1, 1);
    // matplotlibcpp::title("Plots");
    // matplotlibcpp::named_plot("Coeff Before", coeffCopy);
    // matplotlibcpp::legend();
    // matplotlibcpp::subplot(2, 1, 2);
    // matplotlibcpp::named_plot("Coeff After", coefficients);
    // matplotlibcpp::legend();
    // matplotlibcpp::show();

    auto batchesWithHopsCopy = batchesWithHops;

    // Uzupełnianie zerami do rozwiązania w dziedzinie czasu
    for (int i = 0; i < batchesWithHops.size(); ++i)
    {
      std::vector<double> tmp;

      for (int j = 0; j < filterLength - 1; ++j)
      {
        tmp.push_back(0);
      }

      auto tmp2 = batchesWithHops[i];

      tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());

      for (int j = 0; j < filterLength - 1; ++j)
      {
        tmp.push_back(0);
      }

      batchesWithHops[i] = tmp;
    }

    // Nakładanie filtra w dziedzinie czasu
    std::vector<std::vector<double>> timeDomainOutput;

    for (int i = 0; i < batchesWithHops.size(); ++i)
    {
      timeDomainOutput.push_back(std::vector<double>());

      for (int j = filterLength; j < batchesWithHops[i].size() - filterLength; ++j)
      {
        double sum = 0;
        for (int k = 0; k < filterLength; ++k)
        {
          sum += batchesWithHops[i][j - k] * coefficients[k];
        }
        timeDomainOutput[i].push_back(sum);
      }
    }

    // Mergowanie nachodzących okien
    std::vector<double> mergedTimeDomainSignal(audioSource.samples[0].size(), 0);

    for (int i = 0; i < timeDomainOutput.size(); ++i)
    {
      if (i == 0)
      {
        for (int j = 0; j < timeDomainOutput[i].size(); ++j)
        {
          mergedTimeDomainSignal[j] += timeDomainOutput[i][j];
        }
      }
      else
      {
        for (int j = 0; j < timeDomainOutput[i].size(); ++j)
        {
          if (i * hopSize + j >= audioSource.samples[0].size())
          {
            break;
          }
          mergedTimeDomainSignal[i * hopSize + j] += timeDomainOutput[i][j];
        }
      }
    }

    // matplotlibcpp::subplot(2, 1, 1);
    // matplotlibcpp::title("Plots");
    // matplotlibcpp::named_plot("Before", audioSource.samples[0]);
    // matplotlibcpp::legend();
    // matplotlibcpp::subplot(2, 1, 2);
    // matplotlibcpp::named_plot("After", mergedTimeDomainSignal);
    // matplotlibcpp::legend();
    // matplotlibcpp::show();

    // UZUPEŁNIANIE ZERAMI DO FFT

    auto batchesWithHopsWithAddedZeros = batchesWithHopsCopy;
    auto coefficientsWithAddedZeros = coeffCopy;

    int desiredSize = filterLength + windowSize - 1;

    std::cout << "DUPA1" << std::endl;
    while (coefficientsWithAddedZeros.size() != desiredSize)
    {
      coefficientsWithAddedZeros.push_back(0);
    }
    std::cout << "DUPA2" << std::endl;
    for (auto& batch : batchesWithHopsWithAddedZeros)
    {
      while (batch.size() != desiredSize)
      {
        batch.push_back(0);
      }
      assert(coefficientsWithAddedZeros.size() == batch.size());
    }
    std::cout << "DUPA3" << std::endl;
    CArray coefficientsFFT(coefficientsWithAddedZeros.size());
    std::vector<CArray> batchesFFT(batchesWithHopsWithAddedZeros.size(),
                                   CArray(coefficientsWithAddedZeros.size()));

    for (int i = 0; i < coefficientsWithAddedZeros.size(); ++i)
    {
      coefficientsFFT[i] = coefficientsWithAddedZeros[i];
    }

    std::cout << "DUPA4" << std::endl;
    Utility::FFT(coefficientsFFT);
    std::cout << "DUPA5" << std::endl;
    for (int i = 0; i < batchesFFT.size(); ++i)
    {
      for (int j = 0; j < batchesFFT[i].size(); ++j)
      {
        batchesFFT[i][j] = batchesWithHopsWithAddedZeros[i][j];
      }
      Utility::FFT(batchesFFT[i]);
      batchesFFT[i] = batchesFFT[i] * coefficientsFFT;
    }

    for (auto& batch : batchesFFT)
    {
      Utility::IFFT(batch);
    }

    std::vector<double> freqDomainResult(audioSource.samples[0].size());

    for (int i = 0; i < batchesFFT.size(); ++i)
    {
      if (i == 0)
      {
        for (int j = 0; j < batchesFFT[i].size(); ++j)
        {
          freqDomainResult[j] += batchesFFT[i][j].real();
        }
      }
      else
      {
        for (int j = 0; j < batchesFFT[i].size(); ++j)
        {
          if (i * hopSize + j >= audioSource.samples[0].size())
          {
            break;
          }
          freqDomainResult[i * hopSize + j] += batchesFFT[i][j].real();
        }
      }
    }

    matplotlibcpp::subplot(2, 1, 1);
    matplotlibcpp::title("Plots");
    matplotlibcpp::named_plot("Before", audioSource.samples[0]);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(2, 1, 2);
    matplotlibcpp::named_plot("After", freqDomainResult);
    matplotlibcpp::legend();
    matplotlibcpp::show();

    std::cout << "DUPA6" << std::endl;

    // std::vector<double> plotData;
    // matplotlibcpp::named_plot("Data", batchesWithHops[0]);
    // matplotlibcpp::show();

    // do {
    //   std::cout
    //       << "\nSelect window size:\n'1' = 1024\n'2' = 2048\n'3' = 4096\n";
    //   std::cin >> windowSizeChoice;
    // } while (std::find(kPossibleChoices.begin(), kPossibleChoices.end(),
    //                    windowSizeChoice) == kPossibleChoices.end());

    // do {
    //   std::cout << "\nSelect window type:\n'1': Rectangle\n'2': Hamming\n'3':
    //   "
    //                "Hanning\n";
    //   std::cin >> windowFunctionTypeChoice;
    // } while (std::find(kPossibleChoices.begin(), kPossibleChoices.end(),
    //                    windowFunctionTypeChoice) == kPossibleChoices.end());

    // switch (windowSizeChoice) {
    // case '1': {
    //   windowSize = 1024;
    //   break;
    // }
    // case '2': {
    //   windowSize = 2048;
    //   break;
    // }
    // case '3': {
    //   windowSize = 4096;
    //   break;
    // }
    // }

    // switch (windowFunctionTypeChoice) {
    // case '1': {
    //   windowFunctionType = WindowFunctionType::Rectangle;
    //   break;
    // }
    // case '2': {
    //   windowFunctionType = WindowFunctionType::Hamming;
    //   break;
    // }
    // case '3': {
    //   windowFunctionType = WindowFunctionType::Hanning;
    //   break;
    // }
    // }

    // std::cout << "Insert hop size:\n" << std::endl;
    // std::cin >> hopSize;

    // std::cout << "Loaded window size: " << windowSize << std::endl;
    // auto plotData =
    //     Utility::GetSegmentedSamples(audioSource.samples[0], windowSize);
    // std::vector<std::vector<double>> batches = plotData.segmentedPlotData;
    // std::vector<int> mask = plotData.mask;

    // std::vector<double> pitchesAutoCorrelation;

    // for (const auto &batch : batches) {

    //   std::vector<double> transformedSignal;

    //   for (int m = 1; m < batch.size(); ++m) {
    //     double sum = 0;
    //     for (int i = 0; i < batch.size(); ++i) {
    //       double value = batch[i];

    //       value *= m + i < batch.size() ? batch[m + i] : 0;
    //       sum += value;
    //     }
    //     transformedSignal.push_back(sum);
    //   }

    //   // Utility::ShowPlot(batch, "Discrete time signal");
    //   // Utility::ShowPlot(transformedSignal, "Autocorelation of discrete
    //   // signal");

    //   int maxIndexAutocorr = Utility::FindIndexOfMaximum(transformedSignal);
    //   double baseFrequency =
    //       maxIndexAutocorr == -1
    //           ? 0
    //           : static_cast<double>(audioSource.getSampleRate()) /
    //                 (maxIndexAutocorr + 1);

    //   // std::endl; ShowPlot(batch, "base"); ShowPlot(transformedSignal,
    //   // "auto");

    //   pitchesAutoCorrelation.push_back(baseFrequency);
    // }

    // std::vector<double> pitchDataAutocorrelation =
    // Utility::GetBaseFreqPlotData(
    //     audioSource, pitchesAutoCorrelation, windowSize);

    // std::vector<double> pitchesCepstrum;

    // for (auto batch : batches) {

    //   CArray freqDomainData(windowSize);

    //   // Utility::ShowPlot(batch, "Basic signal before windowing");
    //   // okienkowanie przez funkcje Hamminga

    //   Utility::ApplyWindowFunction(batch, windowFunctionType);

    //   // Utility::ShowPlot(batch, "Basic signal windowed by Hamming
    //   function");

    //   for (int i = 0; i < batch.size(); ++i) {
    //     freqDomainData[i] = batch[i];
    //   }

    //   Utility::FFT(freqDomainData);
    //   CArray halfOfFreqData(windowSize / 2);

    //   for (int i = 0; i < windowSize / 2; ++i) {
    //     halfOfFreqData[i] = freqDomainData[i];
    //   }

    //   CArray firstSpectrum(windowSize / 2);
    //   std::vector<double> firstSpectrumVec(windowSize / 2);

    //   for (int i = 0; i < windowSize / 2; ++i) {
    //     double val = log(sqrt(pow(halfOfFreqData[i].real(), 2) +
    //                           pow(halfOfFreqData[i].imag(), 2)) +
    //                      M_E);
    //     firstSpectrum[i] = val;
    //     firstSpectrumVec[i] = val;
    //   }

    //   // ShowPlot(firstSpectrumVec, "Freq domain - log(mag) signal");
    //   Utility::FFT(firstSpectrum);

    //   CArray halfOfFirstSpectrum(windowSize / 4);
    //   for (int i = 0; i < windowSize / 4; ++i) {
    //     halfOfFirstSpectrum[i] = firstSpectrum[i];
    //   }

    //   std::vector<double> secondSpectrum(windowSize / 4);
    //   for (int i = 0; i < windowSize / 4; ++i) {
    //     double val = val = sqrt(pow(halfOfFirstSpectrum[i].real(), 2) +
    //                             pow(halfOfFirstSpectrum[i].imag(), 2));

    //     secondSpectrum[i] = val;
    //   }

    //   int maxIndexCep = Utility::FindIndexOfMaximum(secondSpectrum);

    //   double baseFrequencyCep =
    //       maxIndexCep == -1 ? 0
    //                         :
    //                         static_cast<double>(audioSource.getSampleRate())
    //                         /
    //                               (maxIndexCep * 2);

    //   pitchesCepstrum.push_back(baseFrequencyCep);
    // }

    // std::vector<double> pitchDataCepstrum =
    //     Utility::GetBaseFreqPlotData(audioSource, pitchesCepstrum,
    //     windowSize);

    // assert(audioSource.samples[0].size() == pitchDataAutocorrelation.size());
    // assert(audioSource.samples[0].size() == pitchDataCepstrum.size());

    // for (int i = 0; i < audioSource.samples[0].size(); ++i) {
    //   pitchDataAutocorrelation[i] *= mask[i];
    //   pitchDataCepstrum[i] *= mask[i];
    // }

    // // Utility::ShowPlot(audioSource.samples[0], "Discrete time domain
    // signal");

    // matplotlibcpp::subplot(2, 1, 1);
    // matplotlibcpp::title("Discrete time domain signal");
    // matplotlibcpp::named_plot("Signal", audioSource.samples[0]);
    // matplotlibcpp::legend();
    // matplotlibcpp::subplot(2, 1, 2);
    // matplotlibcpp::named_plot("Pitch detected by autocorrelation",
    //                           pitchDataAutocorrelation);
    // matplotlibcpp::named_plot("Pitch detected by cepstrum",
    // pitchDataCepstrum); matplotlibcpp::title("Pitch detection comparision");
    // matplotlibcpp::legend();
    // matplotlibcpp::show();

    // auto autocorrPitchData =
    //     Utility::GeneratePitchSignal(audioSource, pitchDataAutocorrelation);
    // auto cepstrumPitchData =
    //     Utility::GeneratePitchSignal(audioSource, pitchDataCepstrum);

    // matplotlibcpp::subplot(3, 1, 1);
    // matplotlibcpp::title("Discrete time domain signal");
    // matplotlibcpp::named_plot("Signal", audioSource.samples[0]);
    // matplotlibcpp::legend();

    // matplotlibcpp::subplot(3, 1, 2);
    // matplotlibcpp::title("Sine wave generated from autocorrelation pitch");
    // matplotlibcpp::named_plot("Signal", autocorrPitchData);
    // matplotlibcpp::legend();

    // matplotlibcpp::subplot(3, 1, 3);
    // matplotlibcpp::title("Sine wave generated from cepstrum pitch");
    // matplotlibcpp::named_plot("Signal", cepstrumPitchData);
    // matplotlibcpp::legend();

    // matplotlibcpp::show();

    // AudioFile<double>::AudioBuffer bufferAuto;
    // AudioFile<double>::AudioBuffer bufferCep;

    // bufferAuto.resize(1);
    // bufferCep.resize(1);

    // // 3. Set number of samples per channel
    // bufferAuto[0].resize(autocorrPitchData.size());
    // bufferCep[0].resize(cepstrumPitchData.size());

    // bufferAuto[0] = autocorrPitchData;
    // bufferCep[0] = cepstrumPitchData;

    // AudioFile<double> autoCorrFile;
    // AudioFile<double> cepstrumFile;

    // autoCorrFile.setAudioBuffer(bufferAuto);
    // cepstrumFile.setAudioBuffer(bufferCep);

    // autoCorrFile.setBitDepth(24);
    // cepstrumFile.setBitDepth(24);
    // autoCorrFile.setSampleRate(44100);
    // cepstrumFile.setSampleRate(44100);

    // autoCorrFile.save("./autoCorrPitch.wav", AudioFileFormat::Wave);
    // cepstrumFile.save("./cepstrumPitch.wav", AudioFileFormat::Wave);
    break;
  }
  case 'l':
  case 'L':
  {

    POID_DGMK::Utility::LoadSoundUntilSuccessful(aInputFile, aBaseSound);
    break;
  }

  case 'q':
  case 'Q':
  {
    exit(0);
    break;
  }
  case 's':
  case 'S':
  {
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
  case 'V':
  {
    auto& audioSource = aBaseSound;
    audioSource.printSummary();

    std::cin.ignore();
    std::cin.get();
    break;
  }
  default:
  {
    std::cout << "Unhandled option chosen\n";
    break;
  }
  }
} // namespace POID_DGMK

} // namespace POID_DGMK