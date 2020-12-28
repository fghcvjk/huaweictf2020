#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include "image.h"
#include "stencil.h"

#define P float

int main(int argc, char** argv) {

  ImageClass<P> img_in ("input.png");
  ImageClass<P> img_out(img_in.width, img_in.height);
  
  for (int i = 1; i <= 4; i++) {
    ApplyStencil(img_in, img_out);
    ApplyStencil(img_out, img_in);
  }

  ApplyStencil(img_in, img_out);

  FILE * fp;

  fp = fopen("data.dat", "wb+");
  P * out = img_out.pixel;
  for (int j = 1; j < img_in.width-1; j++) {
      for (int i = 1; i < img_in.height-1; i++) {
          fwrite(&out[i*img_in.width + j], sizeof(P), 1, fp);
      }
  }
  fclose(fp);
}
