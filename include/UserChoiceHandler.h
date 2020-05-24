#ifndef USER_CHOICE_HANDLER_H
#define USER_CHOICE_HANDLER_H

#include "AudioFile.h"
#include "Utility.h"
#include <complex>
#include <string>
#include <valarray>
#include <vector>

namespace POID_DGMK {

class UserChoiceHandler {
public:
  using TComplex = std::complex<double>;
  using CArray = std::valarray<TComplex>;
  using TComplexRepresentation = std::vector<CArray>;

  static void ViewMenu();

  static void HandleUserChoice(std::string &aInputFile,
                               AudioFile<double> &aBaseSound);
};

} // namespace POID_DGMK
#endif // USER_CHOICE_HANDLER_H