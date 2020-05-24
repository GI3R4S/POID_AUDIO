#include "Utility.h"
#include "matplotlibcpp.h"
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdio.h>

#include <algorithm>
#include <assert.h>

namespace POID_DGMK {
namespace plt = matplotlibcpp;

void Utility::FFT(CArray &x) {
  const size_t N = x.size();
  if (N <= 1)
    return;

  CArray even = x[std::slice(0, N / 2, 2)];
  CArray odd = x[std::slice(1, N / 2, 2)];

  FFT(even);
  FFT(odd);

  for (size_t k = 0; k < N / 2; ++k) {
    TComplex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
    x[k] = even[k] + t;
    x[k + N / 2] = even[k] - t;
  }
}

void Utility::IFFT(CArray &x) {
  x = x.apply(std::conj);
  FFT(x);
  x = x.apply(std::conj);
  x /= x.size();
}

std::vector<double>
Utility::GetBaseFreqPlotData(const AudioFile<double> &aAudioSource,
                             const std::vector<double> &aPitches,
                             int aWindowSize) {

  constexpr double kEpsilon = 5;

  const double kMin = *std::min(aPitches.begin(), aPitches.end());
  const double kMax = *std::max(aPitches.begin(), aPitches.end());
  const double kAvg =
      std::accumulate(aPitches.begin(), aPitches.end(), 0) / aPitches.size();

  std::vector<std::vector<double>> pitchesGroups;
  std::vector<double> pitchesGroup;

  for (int i = 0; i < aPitches.size(); ++i) {
    if (pitchesGroup.empty()) {
      pitchesGroup.push_back(aPitches[i]);
    } else {
      auto tmpVec = pitchesGroup;
      tmpVec.push_back(aPitches[i]);
      auto min = *std::min_element(tmpVec.begin(), tmpVec.end());
      auto max = *std::max_element(tmpVec.begin(), tmpVec.end());
      if (max - min <= kEpsilon) {
        pitchesGroup.push_back(aPitches[i]);
      } else {
        pitchesGroups.push_back(pitchesGroup);
        pitchesGroup.clear();
        pitchesGroup.push_back(aPitches[i]);
      }
    }
  }
  if (!pitchesGroup.empty()) {
    pitchesGroups.push_back(pitchesGroup);
  }

  std::vector<double> pitchData;

  auto lastGroupSize = aAudioSource.samples[0].size() % aWindowSize;

  for (int i = 0; i < pitchesGroups.size(); ++i) {
    const auto avg =
        std::accumulate(pitchesGroups[i].begin(), pitchesGroups[i].end(), 0) /
        pitchesGroups[i].size();
    auto groupSize = 0;

    if (i == pitchesGroups.size() - 1) {
      groupSize = (pitchesGroups[i].size() - 1) * aWindowSize + lastGroupSize;
    } else {
      groupSize = pitchesGroups[i].size() * aWindowSize;
    }

    for (int j = 0; j < groupSize; ++j) {
      pitchData.push_back(avg);
    }
  }

  return pitchData;
}

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

void Utility::ShowPlot(const std::vector<double> &aData,
                       const std::string &aName) {
  plt::plot(aData);
  plt::title(aName);
  plt::show();
}

std::vector<std::vector<double>>
Utility::GetSegmentedSamples(const std::vector<double> &aSamples,
                             int aWindowSize) {
  std::vector<std::vector<double>> batches;

  int numberOfFullBatches = aSamples.size() / aWindowSize;
  int numberpOfSamplesInLastBatch = aSamples.size() % aWindowSize;

  for (int i = 1; i <= numberOfFullBatches; ++i) {
    std::vector<double> batch;
    for (int j = (i - 1) * aWindowSize; j < i * aWindowSize; ++j) {
      batch.push_back(aSamples[j]);
    }

    batches.push_back(batch);
  }

  std::vector<double> lastBatch;
  for (int i = numberOfFullBatches * aWindowSize;
       i < numberOfFullBatches * aWindowSize + numberpOfSamplesInLastBatch;
       ++i) {
    lastBatch.push_back(aSamples[i]);
  }
  batches.push_back(lastBatch);

  return batches;
}

} // namespace POID_DGMK