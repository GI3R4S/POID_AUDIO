#ifndef I_IMAGE_PROCESSOR_H
#define I_IMAGE_PROCESSOR_H

#include "CImg.h"
#include "FFTUtilities.h"
#include "Image.h"
#include "ImageDisplayManager.h"

namespace POID_DGMK
{

class IImageProcessor
{
public:
  virtual Image Process(const Image& aImg) = 0;
};
} // namespace POID_DGMK

#endif // I_IMAGE_PROCESSOR_H