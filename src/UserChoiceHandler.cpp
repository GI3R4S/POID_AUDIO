#include "UserChoiceHandler.h"
#include "WindowApplicator.h"
#include "matplotlibcpp.h"
#include <algorithm>
#include <assert.h>
#include <chrono>
#include <complex>
#include <fstream>
#include <map>
#include <numeric>
#include <stdio.h>
#include <valarray>
#include <vector>

namespace POID_DGMK
{

void UserChoiceHandler::ViewMenu()
{
  std::cout << "==============================================\n"
            << "[1] Compare time and frequency domain low pass filter\n"
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
    // configuration variables - start
    auto& audioSource = aBaseSound;
    char windowFunctionTypeChoice;
    int windowSize = 0;
    int hopSize = 0;
    bool isCasuative = 0;
    int filterLength = 0;
    int cutoff_freq = 0;
    int nValue = 0;
    WindowFunctionType windowFunctionType;
    // configuration variables - end

    // input of params - start
    std::cout << "\nInsert size of window: \n";
    std::cin >> windowSize;

    const std::vector<char> kPossibleChoices{'1', '2', '3'};

    do
    {
      std::cout << "\nSelect window type:\n'1': Rectangle\n'2': Hamming\n'3': "
                   "Hann\n";
      std::cin >> windowFunctionTypeChoice;
    } while (std::find(kPossibleChoices.begin(), kPossibleChoices.end(),
                       windowFunctionTypeChoice) == kPossibleChoices.end());

    switch (windowFunctionTypeChoice)
    {
    case '1':
    {
      windowFunctionType = WindowFunctionType::Rectangle;
      break;
    }
    case '2':
    {
      windowFunctionType = WindowFunctionType::Hamming;
      break;
    }
    case '3':
    {
      windowFunctionType = WindowFunctionType::Hann;
      break;
    }
    }

    std::cout << "Insert size of hop: \n";
    std::cin >> hopSize;

    std::cout << "Insert number of coefficients - should be odd number: \n";
    std::cin >> filterLength;

    std::cout << "Insert border frequency value for lowpass filter - should be "
                 "between 1 - 1000: \n";
    std::cin >> cutoff_freq;

    std::cout << "Insert value of N - 0 for default value is first 2^n GT "
                 "M(WINDOW SIZE) + L(NUMBER OF COEFFICIENTS) - 1 = "
              << Utility::GetFirstPowerOf2GT(windowSize + filterLength - 1)
              << std::endl;

    std::cin >> nValue;

    std::cout
      << "Is window casuative ? Type 0 if no, any other number otherwise: \n";
    std::cin >> isCasuative;

    const double kSignalMin = *std::min_element(audioSource.samples[0].begin(),
                                                audioSource.samples[0].end());
    const double kSignalMax = *std::max_element(audioSource.samples[0].begin(),
                                                audioSource.samples[0].end());
    // input of params - end

    // segmentation of signal - common for both variants

    auto segments =
      Utility::GetSegmentedSignal(audioSource.samples[0], windowSize, hopSize);

    WindowApplicator windowApplicator(windowSize, hopSize, windowFunctionType);

    for (int i = 0; i < segments.size(); ++i)
    {
      windowApplicator.ApplyWindowFunction(segments[i]);
    }

    // computation of coefficients - start
    std::vector<double> coefficients;
    int half = (filterLength - 1) / 2.0;
    for (int i = 0; i < filterLength; ++i)
    {
      if (i == half)
      {
        coefficients.push_back(2.0 * cutoff_freq / 44100.0);
      }
      else
      {
        double value = sin(2.0 * M_PI * cutoff_freq / 44100.0 * (i - half)) /
                       (M_PI * (i - half));
        coefficients.push_back(value);
      }
    }

    // computation of coefficients - end

    // applying window function on coefficients
    auto coefficientsBeforeWindowing = coefficients;
    WindowApplicator::ApplyWindowFunction(coefficients, windowFunctionType);
    auto coefficientsAfterWindowing = coefficients;

    auto segmentedSignalTimeDomain = segments;
    auto segmentedSignalFreqDomain = segments;

    ////////////////////////////////////////////////////////////////////////////
    // Time domain solution - start
    auto timeDomainStart = std::chrono::steady_clock::now();
    auto segmentsCopy = segments;

    // Convolution of batches and coefficients - start
    std::vector<std::vector<double>> timeDomainOutput;

    for (int i = 0; i < segmentedSignalTimeDomain.size(); ++i)
    {
      auto convolution =
        Utility::Convolution(segmentedSignalTimeDomain[i], coefficientsAfterWindowing);
      timeDomainOutput.push_back(std::move(convolution));
    }
    // Convolution of batches and coefficients - end

    // Merging modified segments, addition of values on connetction - start
    std::vector<double> timeDomainResult = windowApplicator.MergeDoubleSegments(
      timeDomainOutput, audioSource.samples[0].size());

    // Utility::Normalize(timeDomainResult, kSignalMin, kSignalMax);

    // Merging modified segments, addition of values on connetction - end

    auto timeDomainEnd = std::chrono::steady_clock::now();
    // Time domain solution - end
    ///////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////
    // Freq domain solution - start
    auto freqDomainStart = std::chrono::steady_clock::now();
    auto coefficientsWithAddedZeros = coefficientsAfterWindowing;

    int desiredSize =
      nValue == 0 ? Utility::GetFirstPowerOf2GT(windowSize + filterLength - 1) : nValue;

    // Adding zeroes to coefficients - start
    if (isCasuative)
    {
      while (coefficientsWithAddedZeros.size() != desiredSize)
      {
        coefficientsWithAddedZeros.push_back(0);
      }
    }
    else
    {
      // auto coefficientsCpy = coefficients;
      // Adding zeroes to coefficients - start
      std::vector<double> firstHalf;
      std::vector<double> secondHalf;
      std::vector<double> zeroes(desiredSize - coefficients.size(), 0);

      for (int i = 0; i < coefficientsWithAddedZeros.size(); ++i)
      {
        if (i < coefficientsWithAddedZeros.size() / 2.0)
        {
          secondHalf.push_back(coefficientsWithAddedZeros[i]);
        }
        else
        {
          firstHalf.push_back(coefficientsWithAddedZeros[i]);
        }
      }

      coefficientsWithAddedZeros.clear();
      coefficientsWithAddedZeros.insert(coefficientsWithAddedZeros.end(),
                                        firstHalf.begin(), firstHalf.end());
      coefficientsWithAddedZeros.insert(coefficientsWithAddedZeros.end(),
                                        zeroes.begin(), zeroes.end());
      coefficientsWithAddedZeros.insert(coefficientsWithAddedZeros.end(),
                                        secondHalf.begin(), secondHalf.end());
    }

    // Adding zeroes to segments - start
    for (auto& batch : segmentedSignalFreqDomain)
    {
      while (batch.size() != desiredSize)
      {
        batch.push_back(0);
      }
      assert(coefficientsWithAddedZeros.size() == batch.size());
    }

    // Adding zeroes to segments - end
    // Performing FFT transforms, multiplication of batches and coefficiennts- start
    CArray coefficientsFFT(coefficientsWithAddedZeros.size());
    std::vector<CArray> segmentsFFT(segmentedSignalFreqDomain.size(),
                                    CArray(coefficientsWithAddedZeros.size()));

    for (int i = 0; i < coefficientsWithAddedZeros.size(); ++i)
    {
      coefficientsFFT[i] = coefficientsWithAddedZeros[i];
    }

    Utility::FFT(coefficientsFFT);
    for (int i = 0; i < segmentsFFT.size(); ++i)
    {
      for (int j = 0; j < segmentsFFT[i].size(); ++j)
      {
        segmentsFFT[i][j] = segmentedSignalFreqDomain[i][j];
      }
      Utility::FFT(segmentsFFT[i]);
      segmentsFFT[i] = segmentsFFT[i] * coefficientsFFT;
    }
    // Performing FFT transforms, multiplication of batches and coefficiennts- end

    // Performing IFFT transform to get result - start
    for (auto& batch : segmentsFFT)
    {
      Utility::IFFT(batch);
    }
    // Performing IFFT transform to get result - end

    // Retrieving result - merging segments - start
    std::vector<double> freqDomainResult = windowApplicator.MergeComplexSegments(
      segmentsFFT, audioSource.samples[0].size());

    // Utility::Normalize(freqDomainResult, kSignalMin, kSignalMax);
    // Retrieving result - merging segments - end

    auto freqDomainEnd = std::chrono::steady_clock::now();
    // Freq domain result - end
    ///////////////////////////////////////////////////////////////////////////

    std::stringstream timeDomainDurationString;
    timeDomainDurationString
      << "Time domain duration: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(timeDomainEnd - timeDomainStart)
           .count()
      << " ms\n";

    std::stringstream freqDomainDurationString;
    timeDomainDurationString
      << "Freq domain duration: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(freqDomainEnd - freqDomainStart)
           .count()
      << " ms\n";

    // plotting - start
    std::cout << timeDomainDurationString.str() + freqDomainDurationString.str()
              << std::endl;

    matplotlibcpp::subplot(6, 1, 1);
    matplotlibcpp::title(timeDomainDurationString.str() +
                         freqDomainDurationString.str());
    matplotlibcpp::named_plot("Basic signal", audioSource.samples[0]);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(6, 1, 2);
    matplotlibcpp::named_plot(
      "Coefficients before application of window function", coefficientsBeforeWindowing);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(6, 1, 3);
    matplotlibcpp::named_plot(
      "Coefficients after application of window function", coefficientsAfterWindowing);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(6, 1, 4);
    matplotlibcpp::named_plot("Coefficients after adding zeroes", coefficientsWithAddedZeros);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(6, 1, 5);
    matplotlibcpp::named_plot("Time domain low pass filter", timeDomainResult);
    matplotlibcpp::legend();
    matplotlibcpp::subplot(6, 1, 6);
    matplotlibcpp::named_plot("Frequency domain low pass filter", freqDomainResult);
    matplotlibcpp::legend();
    matplotlibcpp::show();
    // plotting - end

    Utility::SaveSignalToFile("time", timeDomainResult);
    Utility::SaveSignalToFile("freq", freqDomainResult);

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
}

} // namespace POID_DGMK