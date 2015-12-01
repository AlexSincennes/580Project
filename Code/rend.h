#include "disp.h" /* include your own disp.h file (e.g. hw1)*/

#include <vector>

/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(200)//(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(0.0)//(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(0.0)//(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#ifndef GZRENDER
#define GZRENDER
typedef struct {			/* define a renderer */
  GzDisplay		*display;
  GzCamera		camera;
  short		    matlevel;	        /* top of stack - current xform */
  GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
  GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
  GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
  GzColor		flatcolor;          /* color state for flat shaded triangles */
  int			interp_mode;
  int			numlights;
  GzLight		lights[MAX_LIGHTS];
  GzLight		ambientlight;
  GzColor		Ka, Kd, Ks;
  float		    spec;		/* specular power */
  GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
}  GzRender;
#endif

#ifndef GZPROJECT
#define GZPROJECT
typedef struct {
    GzCoord pos;
    GzCoord normal;
    GzColor Ka, Kd, Ks;
} GzVertex;

typedef struct {
    GzVertex* vertices[3];
} GzTriangle;

typedef struct {
    std::vector<GzTriangle*> tris;
} GzWorldSpaceTriangles;
#endif

// Function declaration
// HW2
int GzNewRender(GzRender **render, GzDisplay *display);
int GzFreeRender(GzRender *render);
int GzBeginRender(GzRender	*render);
int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList);
int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList);

// HW3
int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix	matrix);
int GzPopMatrix(GzRender *render);

// HW5
int GzFreeTexture();

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);

//general line equation helper
// where out is {A,B,C}
int findGeneralLineEq(const GzCoord* tail, const GzCoord* head, GzCoord &out);

// sort vertices CW
int sortVerticesCW(const GzCoord* (&tri)[3], GzCoord* (&normals)[3]);

// linear expression evaluation algorithm
int lee(GzRender *render, const GzCoord* (&tri)[3], GzCoord* (&normals)[3], GzTextureIndex* (&tex)[3]);

// Raycast Render:
int raycast_render(GzRender *render, GzWorldSpaceTriangles *tris);

// general line equation helper
// where out is {A,B,C}
int findGeneralPlaneEq(const GzCoord* (&tri)[3], GzCoord &out, float &d);

#define A	0
#define B	1
#define C	2

void normalizeGzCoord(GzCoord &v);

void TracePath(GzRay ray, GzWorldSpaceTriangles *tris, int depth, GzRender* render,GzColor &retColor); 
