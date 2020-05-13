#include "Utility.h"
#include <complex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <valarray>
#include <vector>

int main()
{
  std::string inputFileName;

  AudioFile<double> baseAudioFile;
  AudioFile<double> modifiedAudioFile;
  AudioFile<double> workAudioFile;

  POID_DGMK::Utility::LoadSoundUntilSuccessful(inputFileName, baseAudioFile);

  modifiedAudioFile = baseAudioFile;
  workAudioFile = baseAudioFile;

  while (true)
  {
    std::system("clear");

    POID_DGMK::Utility::ViewMenu();
    POID_DGMK::Utility::HandleUserChoice(inputFileName, baseAudioFile,
                                                   modifiedAudioFile, workAudioFile);
  }
}
