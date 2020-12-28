#ifndef __INCLUDED_IMAGE_H__
#define __INCLUDED_IMAGE_H__

#include <png.h>
#include <omp.h>
#include <stdlib.h>

template <typename P> 
class ImageClass {
 public:
  P* pixel;
  int width;
  int height;

  ImageClass (int const _width, int const _height);
  ImageClass(char const * file_name);
  ~ImageClass();

  void WriteToFile(char const * file_name);
};

#endif
