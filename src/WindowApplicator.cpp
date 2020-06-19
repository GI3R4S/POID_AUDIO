#include "WindowApplicator.h"
#include "Utility.h"

namespace POID_DGMK
{

WindowApplicator::WindowApplicator(int aWindowSize, int aHopSize, WindowFunctionType aWindowFunctionType)
  : mWindowSize(aWindowSize), mHopSize(aHopSize), mWindowFunctionType(aWindowFunctionType)
{
  std::vector<double> signal(44100, 1);

  auto segments = Utility::GetSegmentedSignal(signal, aWindowSize, aHopSize);

  for (int i = 0; i < segments.size(); ++i)
  {
    ApplyWindowFunction(segments[i], aWindowFunctionType);
  }

  auto mergedSignal = MergeDoubleSegments(segments, signal.size());

  mFactor = *std::max_element(mergedSignal.begin(), mergedSignal.end());

  for (auto& value : mergedSignal)
  {
    value /= mFactor;
  }

  matplotlibcpp::named_plot("Shape of window function", segments[0]);
  matplotlibcpp::legend();
  matplotlibcpp::show();

  matplotlibcpp::named_plot("Windows COLA check for chosen configuration", mergedSignal);
  matplotlibcpp::legend();
  matplotlibcpp::show();
}

std::vector<double> WindowApplicator::MergeDoubleSegments(
  const std::vector<std::vector<double>>& aSegments, int aMaxSize)
{
  std::vector<double> result(aMaxSize, 0);

  for (int i = 0; i < aSegments.size(); ++i)
  {
    if (i == 0)
    {
      for (int j = 0; j < aSegments[i].size(); ++j)
      {
        result[j] += aSegments[i][j];
      }
    }
    else
    {
      for (int j = 0; j < aSegments[i].size(); ++j)
      {
        if (i * mHopSize + j >= aMaxSize)
        {
          break;
        }
        result[i * mHopSize + j] += aSegments[i][j];
      }
    }
  }

  return result;
}

std::vector<double> WindowApplicator::MergeComplexSegments(const std::vector<CArray>& aSegments,
                                                           int aMaxSize)
{
  std::vector<double> result(aMaxSize, 0);

  for (int i = 0; i < aSegments.size(); ++i)
  {
    if (i == 0)
    {
      for (int j = 0; j < aSegments[i].size(); ++j)
      {
        result[j] += aSegments[i][j].real();
      }
    }
    else
    {
      for (int j = 0; j < aSegments[i].size(); ++j)
      {
        if (i * mHopSize + j >= result.size())
        {
          break;
        }
        result[i * mHopSize + j] += aSegments[i][j].real();
      }
    }
  }

  return result;
}

void WindowApplicator::ApplyWindowFunction(std::vector<double>& aData)
{
  ApplyWindowFunction(aData, mWindowFunctionType);
  for (auto& value : aData)
  {
    value /= mFactor;
  }
}

void WindowApplicator::ApplyWindowFunction(std::vector<double>& aData,
                                           WindowFunctionType aWindowFunctionType)
{
  switch (aWindowFunctionType)
  {
  case WindowFunctionType::Rectangle:
  {
    break;
  }
  case WindowFunctionType::Hamming:
  {
    std::vector<double> mockSignal(aData.size(), 1);

    for (int i = 0; i < aData.size(); ++i)
    {
      if (i == 0 || i == aData.size() - 1)
      {
        aData[i] *=
          0.5 * (0.54 - 0.46 * cos(2.0 * M_PI * i / (1.0 * (aData.size() - 1.0))));
      }
      else
      {
        aData[i] *= (0.54 - 0.46 * cos(2.0 * M_PI * i / (1.0 * (aData.size() - 1.0))));
      }
    }

    break;
  }
  case WindowFunctionType::Hann:
  {
    for (int i = 0; i < aData.size(); ++i)
    {
      aData[i] *= 0.5 * (1 - cos((2 * M_PI * i) / (aData.size() - 1)));
    }
    break;
  }
  }
}
} // namespace POID_DGMK
