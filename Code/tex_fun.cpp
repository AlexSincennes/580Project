/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include    <time.h>
#include    <stdlib.h>
#include    "rend.h"
#include    <math.h>

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

// gradient vectors at each grid point X,Y
float gradient[50][50][2];

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd); 
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */

  if (u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0) {
      fprintf(stderr, "u or v parameter out of bounds\n");
      return GZ_FAILURE;
  }

  // bilinear interpolation
  // Color(p) = s t C + (1-s) t D + s (1-t) B + (1-s) (1-t) A 
  
  float posx = u * xs;
  float posy = v * ys;

  // x-dir fractional dist
  float s = posx - floor(posx);
  // y-dir fractional dist
  float t = posy - floor(posy);

  GzColor a, b, c, d;
  memcpy(a, image[(int)floor(posx) + xs*(int)floor(posy)], sizeof(GzColor));
  memcpy(b, image[(int)ceil(posx) + xs*(int)floor(posy)], sizeof(GzColor));
  memcpy(c, image[(int)floor(posx) + xs*(int)ceil(posy)], sizeof(GzColor));
  memcpy(d, image[(int)ceil(posx) + xs*(int)ceil(posy)], sizeof(GzColor));


  // check that a,b,c,d not out of bounds
  // if they are, give it the value of a valid texel
  if ((int)ceil(posx) >= 1.0f) {
      memcpy(b, a, sizeof(GzColor));
      memcpy(d, c, sizeof(GzColor));
  }

  if ((int)ceil(posy) >= 1.0f) {
      memcpy(c, a, sizeof(GzColor));
      memcpy(d, b, sizeof(GzColor));
  }

  color[RED] = s * t * a[RED] + (1.0f - s) * t * d[RED] + s*(1.0f - t) * b[RED] + (1.0f - s) * (1.0f - t) * a[RED];
  color[GREEN] = s * t * a[GREEN] + (1.0f - s) * t * d[GREEN] + s*(1.0f - t) * b[GREEN] + (1.0f - s) * (1.0f - t) * a[GREEN];
  color[BLUE] = s * t * a[BLUE] + (1.0f - s) * t * d[BLUE] + s*(1.0f - t) * b[BLUE] + (1.0f - s) * (1.0f - t) * a[BLUE];

  return GZ_SUCCESS;

}

// 2D Perlin Noise & Helpers as per Wikipedia's
// pseudocode section on Perlin Noise

// Function to linearly interpolate between a0 and a1
// Weight w should be in the range [0.0, 1.0]
float lerp(float a0, float a1, float w) {
    return (1.0 - w)*a0 + w*a1;
}

// Computes the dot product of the distance and gradient vectors.
float dotGridGradient(int ix, int iy, float x, float y) {

    // Compute the distance vector
    float dx = x - (double)ix;
    float dy = y - (double)iy;

    // Compute the dot-product
    return (dx*gradient[iy][ix][0] + dy*gradient[iy][ix][1]);
}

// Compute Perlin noise at coordinates x, y
float perlin(float x, float y) {

    // Determine grid cell coordinates
    int x0 = (x > 0.0 ? (int)x : (int)x - 1);
    int x1 = x0 + 1;
    int y0 = (y > 0.0 ? (int)y : (int)y - 1);
    int y1 = y0 + 1;

    // Determine interpolation weights
    // Could also use higher order polynomial/s-curve here
    float sx = x - (double)x0;
    float sy = y - (double)y0;

    // Interpolate between grid point gradients
    float n0, n1, ix0, ix1, value;
    n0 = dotGridGradient(x0, y0, x, y);
    n1 = dotGridGradient(x1, y0, x, y);
    ix0 = lerp(n0, n1, sx);
    n0 = dotGridGradient(x0, y1, x, y);
    n1 = dotGridGradient(x1, y1, x, y);
    ix1 = lerp(n0, n1, sx);
    value = lerp(ix0, ix1, sy);

    return value;
}

// turbulence function using perlin noise
float turbulence(float u, float v) {
    float t = 0;
    
    const int K = 6;

    for (int i = 0; i < K; i++) {
        int pow2i = pow(2.0f, i);
        t += 1.0f / pow2i * abs(perlin(pow2i*u, pow2i*v));
    }

    return t;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
    const int SIZE = 50;

    if (reset) {          /* generate gradient */
        srand(time(NULL));
        for (int x = 0; x < SIZE; x++) {
            for (int y = 0; y < SIZE; y++) {
                GzCoord grad_vec;
                grad_vec[X] = (float)(rand() % 100) / 100.0f;
                grad_vec[Y] = (float)(rand() % 100) / 100.0f;
                grad_vec[Z] = 0;
                normalizeGzCoord(grad_vec);
                gradient[y][x][X] = grad_vec[X];
                gradient[y][x][Y] = grad_vec[Y];
            }
        }
        reset = 0;
    }

    float turbulence_value = turbulence(u,v);
    // show it in red because scene lighting favours red
    color[RED] = turbulence_value * 2.0f;
    color[BLUE] = turbulence_value * 0.25f;
    color[GREEN] = turbulence_value * 0.25f;

    return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

