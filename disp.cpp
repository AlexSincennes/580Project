/*   CS580 HW1 display functions to be completed   */

#include   "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"
#include	<limits>


int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- pass back pointer 
 */
	*framebuffer = (char*)malloc(3 * sizeof(char) * width * height);
	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* HW1.2 create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	
	*display = (GzDisplay*)malloc(sizeof(GzDisplay));
	(*display)->xres = xRes;
	(*display)->yres = yRes;
	(*display)->fbuf = (GzPixel*)malloc(sizeof(GzPixel) * (*display)->xres * (*display)->yres);
	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* HW1.3 clean up, free memory */
	free(display->fbuf);
	free(display);
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* HW1.4 pass back values for a display */
	*xRes = display->xres;
	*yRes = display->yres;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* HW1.5 set everything to some default values - start a new frame */
	for (int i = 0; i < display->xres; i++) {
		for (int j = 0; j < display->yres; j++) {
			display->fbuf[i + display->xres*j] = GzPixel{ GREY, GREY, GREY, 1, INT_MAX };
		}
	}
	return GZ_SUCCESS;
}
 

int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.6 write pixel values into the display */
	// bounds check
	if (i < 0 || i >= display->xres || j < 0 || j >= display->yres)
		return GZ_SUCCESS;

	if (r < 0)
		r = 0;
	else if (r > GZINTESITY_MAX)
		r = GZINTESITY_MAX;

	if (b < 0)
		b = 0;
	else if (b > GZINTESITY_MAX)
		b = GZINTESITY_MAX;

	if (g < 0)
		g = 0;
	else if (g > GZINTESITY_MAX)
		g = GZINTESITY_MAX;

	display->fbuf[i + display->xres*j] = GzPixel{ r, g, b, a, z };
	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.7 pass back a pixel value to the display */

	if (i < 0 || i >= display->xres || j < 0 || j >= display->yres)
		return GZ_SUCCESS;

	GzPixel tmp = display->fbuf[i + display->xres*j];
	if (r != 0) *r = tmp.red;
	if (g != 0) *g = tmp.green;
	if (b != 0) *b = tmp.blue;
	if (a != 0) *a = tmp.alpha;
	if (z != 0) *z = tmp.z;
	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{

/* HW1.8 write pixels to ppm file -- "P6 %d %d 255\r" */
	fprintf(outfile, "P6 %d %d 255\r", display->xres, display->yres);
	for (int j = 0; j < display->yres; ++j) {
		for (int i = 0; i < display->xres; ++i) {
			unsigned char pixel[3];
			pixel[0] = display->fbuf[i + display->xres*j].red >> 4;
			pixel[1] = display->fbuf[i + display->xres*j].green >> 4;
			pixel[2] = display->fbuf[i + display->xres*j].blue >> 4;
			fwrite(pixel, 1, 3, outfile);
		}
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

/* HW1.9 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	for (int j = 0; j < display->yres; ++j) {
		for (int i = 0; i < display->xres; ++i) {
			framebuffer[3 * (i + display->xres*j)] = display->fbuf[i + display->xres*j].blue >> 4;
			framebuffer[3 * (i + display->xres*j) + 1] = display->fbuf[i + display->xres*j].green >> 4;
			framebuffer[3 * (i + display->xres*j)+2] = display->fbuf[i + display->xres*j].red >> 4;
		}
	}
	return GZ_SUCCESS;
}