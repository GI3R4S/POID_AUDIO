#include "UserChoiceHandler.h"
#include "Utility.h"
#include <complex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <valarray>
#include <vector>

int main() {
  std::string inputFileName;

  AudioFile<double> baseAudioFile;

  POID_DGMK::Utility::LoadSoundUntilSuccessful(inputFileName, baseAudioFile);

  while (true) {
    std::system("clear");

    POID_DGMK::UserChoiceHandler::ViewMenu();
    POID_DGMK::UserChoiceHandler::HandleUserChoice(
        inputFileName, baseAudioFile);
  }
}
