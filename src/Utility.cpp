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

bool Utility::LoadSound(std::string& aFileName, AudioFile<double>& aSoundToUpdate)
{
  std::cout << "Insert '1' for 'lin.wav', '2' for 'log.wav', , '3' for "
               "'piano.wav', '4' for 'violin.wav': \n";
  std::cin >> aFileName;

  std::string filePath;
  if (aFileName == "1")
  {
    filePath.append(RESOURCES_DIR);
    filePath.append("lin.wav");
  }
  else if (aFileName == "2")
  {
    filePath.append(RESOURCES_DIR);
    filePath.append("log.wav");
  }
  else if (aFileName == "3")
  {
    filePath.append(RESOURCES_DIR);
    filePath.append("piano.wav");
  }
  else if (aFileName == "4")
  {
    filePath.append(RESOURCES_DIR);
    filePath.append("violin.wav");
  }
  else
  {
    return false;
  }

  return aSoundToUpdate.load(filePath);
}

void Utility::LoadSoundUntilSuccessful(std::string& aFileName, AudioFile<double>& aSoundToUpdate)
{
  if (!LoadSound(aFileName, aSoundToUpdate))
  {
    LoadSoundUntilSuccessful(aFileName, aSoundToUpdate);
  }
}

void Utility::SaveSignalToFile(std::string aFileName, const std::vector<double>& aData)
{
  AudioFile<double>::AudioBuffer buffer;
  buffer.resize(1);

  buffer[0].resize(aData.size());
  buffer[0] = aData;

  AudioFile<double> bufferFile;

  bufferFile.setAudioBuffer(buffer);

  bufferFile.setBitDepth(24);
  bufferFile.setSampleRate(44100);

  bufferFile.save("./" + aFileName + ".wav", AudioFileFormat::Wave);
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
    int firstIndexOfBatch = m * aHopSize;
    int maxIndexOfBatch = m * aHopSize + aWindowSize;

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
    }
    else
    {
      for (int i = firstIndexOfBatch; i < maxIndexOfBatch; ++i)
      {
        batch.push_back(aSamples[i]);
      }
    }
    batches.push_back(batch);
  }

  return batches;
}

} // namespace POID_DGMK