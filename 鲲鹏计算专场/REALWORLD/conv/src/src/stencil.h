#ifndef __INCLUDED_STENCIL_H__
#define __INCLUDED_STENCIL_H__

#include "image.h"
#include <omp.h>

template<typename P>
void ApplyStencil(ImageClass<P> & img_in, ImageClass<P> & img_out);

#endif
