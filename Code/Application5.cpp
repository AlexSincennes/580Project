// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"
#include "time.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "ppot.asc"
#define OUTFILE "output.ppm"

extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */

void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 
 	m_nWidth = 256;		// frame buffer and display width
	m_nHeight = 256;    // frame buffer and display height

	status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

	status |= GzNewDisplay(&m_pDisplay, m_nWidth, m_nHeight);

	status |= GzGetDisplayParams(m_pDisplay, &xRes, &yRes); 
 
	status |= GzNewRender(&m_pRender, m_pDisplay); 

#if 0 	/* set up app-defined camera if desired, else use camera defaults */
    camera.position[X] = -3;
    camera.position[Y] = -25;
    camera.position[Z] = -4;

    camera.lookat[X] = 7.8;
    camera.lookat[Y] = 0.7;
    camera.lookat[Z] = 6.5;

    camera.worldup[X] = -0.2;
    camera.worldup[Y] = 1.0;
    camera.worldup[Z] = 0.0;

    camera.FOV = 63.7;              /* degrees *              /* degrees */

	status |= GzPutCamera(m_pRender, &camera); 
#endif 

	/* Start Renderer */
	status |= GzBeginRender(m_pRender);

	/* Light */
	GzLight	light1 = { { -0.7071, 0.7071, 0 }, { 0.9, 0.9, 0.9 } };
	GzLight	light2 = { { 0, -0.7071, -0.7071 }, { 0.9, 0.9, 0.9 } };
	GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.9, 0.9, 0.9} };
	GzLight	ambientlight = { {0, 0, 0}, {0.9, 0.9, 0.9} };
	/* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
	GzColor ambientCoefficient = { 0.1, 0.1, 0.1 };
	GzColor diffuseCoefficient = {0.7, 0.7, 0.7};

/* 
  renderer is ready for frame --- define lights and shader at start of frame 
*/

        /*
         * Tokens associated with light parameters
         */
        nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[0] = (GzPointer)&light1;
        nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[1] = (GzPointer)&light2;
        nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[2] = (GzPointer)&light3;
        //status |= GzPutAttribute(m_pRender, 1, nameListLights, valueListLights);
		status |= GzPutAttribute(m_pRender, 3, nameListLights, valueListLights);


        nameListLights[0] = GZ_AMBIENT_LIGHT;
        valueListLights[0] = (GzPointer)&ambientlight;
        status |= GzPutAttribute(m_pRender, 1, nameListLights, valueListLights);

        /*
         * Tokens associated with shading 
         */
        nameListShader[0]  = GZ_DIFFUSE_COEFFICIENT;
        valueListShader[0] = (GzPointer)diffuseCoefficient;

	/* 
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode  
	*/
        nameListShader[1]  = GZ_INTERPOLATE;
        interpStyle = GZ_COLOR;         /* Phong shading */
        valueListShader[1] = (GzPointer)&interpStyle;

        nameListShader[2]  = GZ_AMBIENT_COEFFICIENT;
        valueListShader[2] = (GzPointer)ambientCoefficient;
        nameListShader[3]  = GZ_SPECULAR_COEFFICIENT;
        valueListShader[3] = (GzPointer)specularCoefficient;
        nameListShader[4]  = GZ_DISTRIBUTION_COEFFICIENT;
        specpower = 32;
        valueListShader[4] = (GzPointer)&specpower;


        status |= GzPutAttribute(m_pRender, 5, nameListShader, valueListShader);


	if (status) exit(GZ_FAILURE); 

	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Render() 
{
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */ 
	char		dummy[256]; 
	int			status; 

	time_t start, stop;
	time(&start);

	/* Initialize Display */
	status |= GzInitDisplay(m_pDisplay); 
	
	/* 
	* Tokens associated with triangle vertex values 
	*/ 
	//nameListTriangle[2] = GZ_TEXTURE_INDEX;  

	// I/O File open
	FILE *infile;
	if( (infile  = fopen( INFILE , "r" )) == NULL )
	{
         AfxMessageBox( "The input file was not opened\n" );
		 return GZ_FAILURE;
	}

	FILE *outfile;
	if( (outfile  = fopen( OUTFILE , "wb" )) == NULL )
	{
         AfxMessageBox( "The output file was not opened\n" );
		 return GZ_FAILURE;
	}

	/* 
	* Walk through the list of triangles, set color 
	* and render each triangle 
	*/ 

    GzWorldSpaceTriangles* wst = new GzWorldSpaceTriangles();

    // set kd,ks,ka (copied from above for now)
    /* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 }; 
    GzColor ambientCoefficient = { 0.8, 0.1, 0.1 };
    GzColor diffuseCoefficient = { 0.7, 0.7, 0.7 }; 
	GzColor emissionCoefficient = { 0.7, 0.7, 0.7 };

	int left_wall = 2;
	int right_wall = 5;
	int back_wall = 8;
	int floor_wall = 11;
	int ceil_wall = 14;
	int light_plane = 17;
	int front_wall = 20;



	int counter = 0;
	while( fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[0][0]), &(vertexList[0][1]),  
		&(vertexList[0][2]), 
		&(normalList[0][0]), &(normalList[0][1]), 	
		&(normalList[0][2]), 
		&(uvList[0][0]), &(uvList[0][1]) ); 
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[1][0]), &(vertexList[1][1]), 	
		&(vertexList[1][2]), 
		&(normalList[1][0]), &(normalList[1][1]), 	
		&(normalList[1][2]), 
		&(uvList[1][0]), &(uvList[1][1]) ); 
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[2][0]), &(vertexList[2][1]), 	
		&(vertexList[2][2]), 
		&(normalList[2][0]), &(normalList[2][1]), 	
		&(normalList[2][2]), 
		&(uvList[2][0]), &(uvList[2][1]) ); 

		if (counter <= left_wall)
		{
			
			ambientCoefficient[0] = 0;
			ambientCoefficient[1] = 1;
			ambientCoefficient[2] = 0;

			diffuseCoefficient[0] = .25;
			diffuseCoefficient[1] = .75;
			diffuseCoefficient[2] = .25;

			emissionCoefficient[0] = 0.0f;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			
		}
		else if (counter <= right_wall)
		{
			ambientCoefficient[0] = 0;
			ambientCoefficient[1] = 0;
			ambientCoefficient[2] = 1;

			emissionCoefficient[0] = 0.0;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			diffuseCoefficient[0] = .25;
			diffuseCoefficient[1] = .25;
			diffuseCoefficient[2] = .75;
		}
		else if (counter <= back_wall)
		{
			ambientCoefficient[0] = .75;
			ambientCoefficient[1] = .75;
			ambientCoefficient[2] = .75;

			emissionCoefficient[0] = 0.0f;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			diffuseCoefficient[0] = .75;
			diffuseCoefficient[1] = .75;
			diffuseCoefficient[2] = .75;
		}
		else if (counter <= floor_wall)
		{
			ambientCoefficient[0] = .75;
			ambientCoefficient[1] = .75;
			ambientCoefficient[2] = .75;

			emissionCoefficient[0] = 0.0f;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			diffuseCoefficient[0] = .75;
			diffuseCoefficient[1] = .75;
			diffuseCoefficient[2] = .75;
		}
		else if (counter <= ceil_wall)
		{
			ambientCoefficient[0] = .75;
			ambientCoefficient[1] = .75;
			ambientCoefficient[2] = .75;

			emissionCoefficient[0] = 0.0f;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			diffuseCoefficient[0] = .75;
			diffuseCoefficient[1] = .75;
			diffuseCoefficient[2] = .75;
		}
		else if (counter <= light_plane)
		{
			ambientCoefficient[0] = 1;
			ambientCoefficient[1] = 1;
			ambientCoefficient[2] = 1;

			emissionCoefficient[0] = 2.0f;
			emissionCoefficient[1] = 2.0f;
			emissionCoefficient[2] = 2.0f;

			diffuseCoefficient[0] = 1;
			diffuseCoefficient[1] = 1;
			diffuseCoefficient[2] = 1;
		}
		else if (counter <= front_wall)
		{
			ambientCoefficient[0] = .75;
			ambientCoefficient[1] = .75;
			ambientCoefficient[2] = .75;

			emissionCoefficient[0] = 0.0f;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			diffuseCoefficient[0] = .75;
			diffuseCoefficient[1] = .75;
			diffuseCoefficient[2] = .75;
		}
		else //Cube
		{
			ambientCoefficient[0] = .4;
			ambientCoefficient[1] = .3;
			ambientCoefficient[2] = .3;

			emissionCoefficient[0] = 0.0f;
			emissionCoefficient[1] = 0.0f;
			emissionCoefficient[2] = 0.0f;

			diffuseCoefficient[0] = .4;
			diffuseCoefficient[1] = .3;
			diffuseCoefficient[2] = .3;
		}

        GzTriangle* tri = new GzTriangle();
        GzVertex* v0 = new GzVertex();
        memcpy(v0->Ka, ambientCoefficient, sizeof(GzColor));
        memcpy(v0->Kd, diffuseCoefficient, sizeof(GzColor));
        memcpy(v0->Ks, specularCoefficient, sizeof(GzColor));
		memcpy(v0->Le, emissionCoefficient, sizeof(GzColor));
        memcpy(v0->pos, vertexList[0], sizeof(GzCoord));
        memcpy(v0->normal, normalList[0], sizeof(GzCoord));

        GzVertex* v1 = new GzVertex();
        memcpy(v1->Ka, ambientCoefficient, sizeof(GzColor));
        memcpy(v1->Kd, diffuseCoefficient, sizeof(GzColor));
        memcpy(v1->Ks, specularCoefficient, sizeof(GzColor));
		memcpy(v1->Le, emissionCoefficient, sizeof(GzColor));
        memcpy(v1->pos, vertexList[1], sizeof(GzCoord));
        memcpy(v1->normal, normalList[1], sizeof(GzCoord));

        GzVertex* v2 = new GzVertex();
        memcpy(v2->Ka, ambientCoefficient, sizeof(GzColor));
        memcpy(v2->Kd, diffuseCoefficient, sizeof(GzColor));
        memcpy(v2->Ks, specularCoefficient, sizeof(GzColor));
		memcpy(v2->Le, emissionCoefficient, sizeof(GzColor));
        memcpy(v2->pos, vertexList[2], sizeof(GzCoord));
        memcpy(v2->normal, normalList[2], sizeof(GzCoord));

        tri->vertices[0] = v0;
        tri->vertices[1] = v1;
        tri->vertices[2] = v2;

        wst->tris.push_back(tri);
		counter++;
	} 

    // DO RAYCAST HERE

    // e.g.
	raycast_render(m_pRender, wst);

	GzFlushDisplay2File(outfile, m_pDisplay); 	/* write out or update display to file*/
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, m_pDisplay);	// write out or update display to frame buffer

	time(&stop);
	double diff = difftime(stop,start);
	int hrs = (int)diff / 3600;
	int mins = ((int)diff / 60) - (hrs * 60);
	int secs = (int)diff - (hrs * 3600) - (mins*60);
	CString str;
	str.Format("\r\nTime taken: %i hrs, %i mins, %i secs \n\n", hrs, mins, secs);
	/* 
	 * Close file
	 */ 

	if( fclose( infile ) )
      AfxMessageBox( "The input file was not closed\n" );

	if( fclose( outfile ) )
      AfxMessageBox( "The output file was not closed\n" );


	AfxMessageBox(str);
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	status |= GzFreeRender(m_pRender); 
	status |= GzFreeDisplay(m_pDisplay);
	status |= GzFreeTexture();
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}



