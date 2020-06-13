#include "Utility.h"
#include "matplotlibcpp.h"
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdio.h>

#include <algorithm>
#include <assert.h>
#include <tuple>

namespace POID_DGMK
{
namespace plt = matplotlibcpp;

void Utility::FFT(CArray& x)
{
  const size_t N = x.size();
  if (N <= 1)
    return;

  CArray even = x[std::slice(0, N / 2, 2)];
  CArray odd = x[std::slice(1, N / 2, 2)];

  FFT(even);
  FFT(odd);

  for (size_t k = 0; k < N / 2; ++k)
  {
    TComplex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
    x[k] = even[k] + t;
    x[k + N / 2] = even[k] - t;
  }
}

void Utility::IFFT(CArray& x)
{
  x = x.apply(std::conj);
  FFT(x);
  x = x.apply(std::conj);
  x /= x.size();
}

std::vector<double> Utility::GetBaseFreqPlotData(const AudioFile<double>& aAudioSource,
                                                 const std::vector<double>& aPitches,
                                                 int aWindowSize)
{

  constexpr double kEpsilon = 5;

  const double kMin = *std::min(aPitches.begin(), aPitches.end());
  const double kMax = *std::max(aPitches.begin(), aPitches.end());
  const double kAvg =
    std::accumulate(aPitches.begin(), aPitches.end(), 0) / aPitches.size();

  std::vector<std::vector<double>> pitchesGroups;
  std::vector<double> pitchesGroup;

  for (int i = 0; i < aPitches.size(); ++i)
  {
    if (pitchesGroup.empty())
    {
      pitchesGroup.push_back(aPitches[i]);
    }
    else
    {
      auto tmpVec = pitchesGroup;
      tmpVec.push_back(aPitches[i]);
      auto min = *std::min_element(tmpVec.begin(), tmpVec.end());
      auto max = *std::max_element(tmpVec.begin(), tmpVec.end());
      if (max - min <= kEpsilon)
      {
        pitchesGroup.push_back(aPitches[i]);
      }
      else
      {
        pitchesGroups.push_back(pitchesGroup);
        pitchesGroup.clear();
        pitchesGroup.push_back(aPitches[i]);
      }
    }
  }
  if (!pitchesGroup.empty())
  {
    pitchesGroups.push_back(pitchesGroup);
  }

  std::vector<double> pitchData;

  auto lastGroupSize = aAudioSource.samples[0].size() % aWindowSize;

  for (int i = 0; i < pitchesGroups.size(); ++i)
  {
    const auto avg =
      std::accumulate(pitchesGroups[i].begin(), pitchesGroups[i].end(), 0) /
      pitchesGroups[i].size();
    auto groupSize = 0;

    if (i == pitchesGroups.size() - 1)
    {
      groupSize = (pitchesGroups[i].size() - 1) * aWindowSize + lastGroupSize;
    }
    else
    {
      groupSize = pitchesGroups[i].size() * aWindowSize;
    }

    for (int j = 0; j < groupSize; ++j)
    {
      pitchData.push_back(avg);
    }
  }

  return pitchData;
}

bool Utility::LoadSound(std::string& aFileName, AudioFile<double>& aSoundToUpdate)
{
  std::cout << "Insert file name from directory 'resources': \n";
  std::cin >> aFileName;

  std::string filePath;
  filePath.append(RESOURCES_DIR);
  filePath.append(aFileName);

  return aSoundToUpdate.load(filePath);
}

void Utility::LoadSoundUntilSuccessful(std::string& aFileName, AudioFile<double>& aSoundToUpdate)
{
  if (!LoadSound(aFileName, aSoundToUpdate))
  {
    LoadSoundUntilSuccessful(aFileName, aSoundToUpdate);
  }
}

void Utility::ShowPlot(const std::vector<double>& aData, const std::string& aName)
{
  plt::plot(aData);
  plt::title(aName);
  plt::show();
}

PlotData Utility::GetSegmentedSamples(const std::vector<double>& aSamples, int aWindowSize)
{

  std::vector<double> samplesCopy(aSamples.begin(), aSamples.end());

  double average = 0;
  for (int i = 0; i < samplesCopy.size(); ++i)
  {
    average += std::abs(samplesCopy[i]);
  }

  average /= samplesCopy.size();
  average *= 1.5;

  std::cout << "average: " << average << std::endl;
  std::vector<int> mask;
  std::vector<std::vector<double>> batches;

  int numberOfFullBatches = aSamples.size() / aWindowSize;
  int numberpOfSamplesInLastBatch = aSamples.size() % aWindowSize;

  bool wasPreviousPowerNotSufficent = false;

  for (int i = 1; i <= numberOfFullBatches; ++i)
  {
    std::vector<double> batch;
    for (int j = (i - 1) * aWindowSize; j < i * aWindowSize; ++j)
    {
      batch.push_back(aSamples[j]);
    }

    bool isMinimalValueNotExceeded =
      *std::max_element(
        batch.begin(),
        batch.end())<average&& * std::min_element(batch.begin(), batch.end())> -
      1 * average;

    int valToInsert = isMinimalValueNotExceeded || wasPreviousPowerNotSufficent ? 0 : 1;

    wasPreviousPowerNotSufficent = isMinimalValueNotExceeded;

    for (int j = 0; j < batch.size(); ++j)
    {
      mask.push_back(valToInsert);
    }

    batches.push_back(batch);
  }

  std::vector<double> lastBatch;
  for (int i = numberOfFullBatches * aWindowSize;
       i < numberOfFullBatches * aWindowSize + numberpOfSamplesInLastBatch; ++i)
  {
    lastBatch.push_back(aSamples[i]);
  }

  bool isMinimalValueNotExceeded =
    *std::max_element(
      lastBatch.begin(),
      lastBatch.end())<average&& * std::min_element(lastBatch.begin(), lastBatch.end())> -
    1 * average;

  int valToInsert = isMinimalValueNotExceeded ? 0 : 1;

  for (int j = 0; j < lastBatch.size(); ++j)
  {
    mask.push_back(valToInsert);
  }

  batches.push_back(lastBatch);

  PlotData toReturn;
  toReturn.segmentedPlotData = batches;
  toReturn.mask = mask;

  return toReturn;
}

std::vector<double> Utility::GeneratePitchSignal(const AudioFile<double>& aSource,
                                                 const std::vector<double>& aPitchData)
{
  assert(aPitchData.size() > 1);
  std::vector<double> pitchSignal;

  double prevPitchValue = std::numeric_limits<double>::lowest();
  int numberOfSamplesPerPeriod = 0;
  for (int i = 0; i < aPitchData.size(); ++i)
  {
    if (prevPitchValue != aPitchData[i])
    {
      prevPitchValue = aPitchData[i];
      numberOfSamplesPerPeriod = aSource.getSampleRate() / aPitchData[i];
    }

    if (aPitchData[i] != 0)
    {
      int modulo = i % numberOfSamplesPerPeriod;
      double factor = 1.0 * modulo / numberOfSamplesPerPeriod;
      double arg = (2 * M_PI) * factor;
      double sinVal = sin(arg);
      pitchSignal.push_back(sinVal);
    }
    else
    {
      pitchSignal.push_back(0);
    }
  }

  return pitchSignal;
}

int Utility::FindIndexOfMaximum(const std::vector<double> aData)
{
  std::vector<double> reducedVector;
  bool isIncreaseReached = false;
  int i = 0;
  for (; i < aData.size(); ++i)
  {
    if (aData[i + 1] > aData[i])
    {
      isIncreaseReached = true;
    }

    if (isIncreaseReached)
    {
      reducedVector.push_back(aData[i]);
    }
  }

  if (!reducedVector.empty())
  {
    auto maxElement = std::max_element(reducedVector.begin(), reducedVector.end());
    auto maxIndex =
      (aData.size() - reducedVector.size()) + (maxElement - reducedVector.begin());

    return maxIndex;
  }
  else
  {
    return -1;
  }
}

void Utility::ApplyWindowFunction(std::vector<double>& aData, WindowFunctionType aWindowFunctionType)
{
  switch (aWindowFunctionType)
  {
  case WindowFunctionType::Rectangle:
  {
    break;
  }
  case WindowFunctionType::Hamming:
  {
    for (int i = 0; i < aData.size(); ++i)
    {
      if (i == 0 || i == aData.size() - 1)
      {
        aData[i] *=
          0.5 * (0.53836 - 0.46164 * cos(2 * M_PI * i / (1.0 * (aData.size() - 1))));
      }
      else
      {
        aData[i] *=
          (0.53836 - 0.46164 * cos(2 * M_PI * i / (1.0 * (aData.size() - 1))));
      }
    }

    break;
  }
  case WindowFunctionType::Hanning:
  {
    for (int i = 0; i < aData.size(); ++i)
    {
      if (i == 0 || i == aData.size() - 1)
      {
        aData[i] *= 0.5 * 0.5 * (1 - cos((2 * M_PI * i) / (aData.size() - 1)));
      }
      else
      {
        aData[i] *= 0.5 * (1 - cos((2 * M_PI * i) / (aData.size() - 1)));
      }
    }
    break;
  }
  }
}

std::vector<std::vector<double>> Utility::GetSegmentedSignal(const std::vector<double>& aSamples,
                                                             int aWindowSize,
                                                             int aHopSize)

{
  std::vector<std::vector<double>> batches;
  int numberOfHops = aSamples.size() / aHopSize;

  for (int m = 0; m <= numberOfHops; ++m)
  {
    std::cout << "=======================================" << std::endl;
    int firstIndexOfBatch = m * aHopSize;
    int maxIndexOfBatch = m * aHopSize + aWindowSize;
    std::cout << "First index: " << firstIndexOfBatch << std::endl;
    std::cout << "Max index: " << maxIndexOfBatch << std::endl;

    std::vector<double> batch;

    if (maxIndexOfBatch > aSamples.size())
    {
      for (int i = firstIndexOfBatch; i < aSamples.size(); ++i)
      {
        batch.push_back(aSamples[i]);
      }
      int samplesToFill = aWindowSize - batch.size();
      for (int i = 0; i < samplesToFill; ++i)
      {
        batch.push_back(0);
      }
      std::cout << "Ups, requires filling: " << std::endl;
    }
    else
    {
      for (int i = firstIndexOfBatch; i < maxIndexOfBatch; ++i)
      {
        batch.push_back(aSamples[i]);
      }
    }
    std::cout << "Batch size: " << batch.size() << std::endl;
    std::cout << "=======================================" << std::endl;
    batches.push_back(batch);
  }

  return batches;
}

} // namespace POID_DGMK