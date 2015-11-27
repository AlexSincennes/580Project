/* CS580 Homework 4 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"gz.h"
#include	"rend.h"
#include	<algorithm>
#include	<limits>

#define PI 3.141592653589793f
#define MAXREFLECTANCE 10

// HW3 helpers
float degToRad(float deg) {
    return deg * PI / 180.0f;
}

void normalizeGzCoord(GzCoord &v) {

    float vlength = sqrt(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
    if (vlength != 0) {
        v[X] = v[X] / vlength;
        v[Y] = v[Y] / vlength;
        v[Z] = v[Z] / vlength;
    }
}

float dotProduct(GzCoord v1, GzCoord v2) {
    return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}

void crossProduct(GzCoord &v1, GzCoord &v2, GzCoord* out) {

    GzCoord ret = { v1[Y] * v2[Z] - v1[Z] * v2[Y], v1[Z] * v2[X] - v1[X] * v2[Z], v1[X] * v2[Y] - v1[Y] * v2[X] };
    (*out)[0] = ret[0];
    (*out)[1] = ret[1];
    (*out)[2] = ret[2];
}

void MVmultiply(const GzMatrix &mat, const GzCoord &v, GzCoord *out) {

    float vh[4] = { v[0], v[1], v[2], 1 };

    float result[4] = { 0, 0, 0, 0 };

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i] += mat[i][j] * vh[j];
        }
    }

    GzCoord result_3d = { result[0] / result[3], result[1] / result[3], result[2] / result[3] };

    (*out)[0] = result_3d[0];
    (*out)[1] = result_3d[1];
    (*out)[2] = result_3d[2];
}

void MMmultiply(const GzMatrix& m1, const GzMatrix& m2, GzMatrix* out) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            float sum = 0;
            for (int k = 0; k < 4; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            (*out)[i][j] = sum;
        }
    }
}

//HW4 helpers
float length(GzCoord &c) {
    float ret = std::sqrt(c[X] * c[X] + c[Y] * c[Y] + c[Z] * c[Z]);
    return ret;
}

// triangle verts argument order is the same sorting order as LEE (clockwise, topmost first)
void barycentricCoords(const GzCoord &vert0, const GzCoord &vert1, const GzCoord &vert2, GzCoord &point, GzCoord* out) {

    //using Cramer's rule (from Christer Ericson's Real-Time Collision Detection)
    GzCoord vert0z, vert1z, vert2z, pointz;
    memcpy(vert0z, vert0, sizeof(GzCoord));
    memcpy(vert1z, vert1, sizeof(GzCoord));
    memcpy(vert2z, vert2, sizeof(GzCoord));
    memcpy(pointz, point, sizeof(GzCoord));

    //reduce z range to prevent overflow
    vert0z[Z] /= 100000;
    vert1z[Z] /= 100000;
    vert2z[Z] /= 100000;
    pointz[Z] /= 100000;

    GzCoord v0, v1, v2;
    v0[X] = vert0z[X] - vert2z[X];
    v0[Y] = vert0z[Y] - vert2z[Y];
    v0[Z] = vert0z[Z] - vert2z[Z];
    
    
    v1[X] = vert1z[X] - vert2z[X];
    v1[Y] = vert1z[Y] - vert2z[Y];
    v1[Z] = vert1z[Z] - vert2z[Z];
    
    
    v2[X] = pointz[X] - vert2z[X];
    v2[Y] = pointz[Y] - vert2z[Y];
    v2[Z] = pointz[Z] - vert2z[Z];

    float v0DOTv0 = dotProduct(v0, v0);
    float v0DOTv1 = dotProduct(v0, v1);
    float v1DOTv1 = dotProduct(v1, v1);
    float v2DOTv0 = dotProduct(v2, v0);
    float v2DOTv1 = dotProduct(v2, v1);

    float d     =  v0DOTv0 * v1DOTv1 - v0DOTv1 * v0DOTv1;
    (*out)[X]   = (v1DOTv1 * v2DOTv0 - v0DOTv1 * v2DOTv1) / d;
    (*out)[Y]   = (v0DOTv0 * v2DOTv1 - v0DOTv1 * v2DOTv0) / d;
    (*out)[Z]   = 1.0f - (*out)[X] - (*out)[Y];
}



int GzRotXMat(float degree, GzMatrix mat) {
    // Create rotate matrix : rotate along x axis
    // Pass back the matrix using mat value

    float rad = degToRad(degree);

    float cosD = cos(rad);
    float sinD = sin(rad);

    GzMatrix m =
    {
        1.0f, 0, 0, 0,
        0, cosD, -sinD, 0,
        0, sinD, cosD, 0,
        0, 0, 0, 1.0f
    };

    memcpy(mat, m, sizeof(GzMatrix));

    return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat) {
    // Create rotate matrix : rotate along y axis
    // Pass back the matrix using mat value

    float rad = degToRad(degree);

    float cosD = cos(rad);
    float sinD = sin(rad);

    GzMatrix m =
    {
        cosD, 0, sinD, 0,
        0, 1, 0, 0,
        -sinD, 0, cosD, 0,
        0, 0, 0, 1.0f
    };

    memcpy(mat, m, sizeof(GzMatrix));

    return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat) {
    // Create rotate matrix : rotate along z axis
    // Pass back the matrix using mat value

    float rad = degToRad(degree);

    float cosD = cos(rad);
    float sinD = sin(rad);

    GzMatrix m =
    {
        cosD, -sinD, 0, 0,
        sinD, cosD, 0, 0,
        0, 0, 1.0f, 0,
        0, 0, 0, 1.0f
    };

    memcpy(mat, m, sizeof(GzMatrix));

    return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat) {
    // Create translation matrix
    // Pass back the matrix using mat value

    GzMatrix m =
    {
        1.0f, 0, 0, translate[X],
        0, 1.0f, 0, translate[Y],
        0, 0, 1.0f, translate[Z],
        0, 0, 0, 1.0f
    };

    memcpy(mat, m, sizeof(GzMatrix));

    return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat) {
    // Create scaling matrix
    // Pass back the matrix using mat value

    GzMatrix m =
    {
        scale[X], 0, 0, 0,
        0, scale[Y], 0, 0,
        0, 0, scale[Z], 0,
        0, 0, 0, 1.0f
    };

    memcpy(mat, m, sizeof(GzMatrix));

    return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay	*display) {
    /*
    - malloc a renderer struct
    - setup Xsp and anything only done once
    - save the pointer to display
    - init default camera
    */

    *render = (GzRender*)malloc(sizeof(GzRender));
    (*render)->display = display;

    (*render)->camera = GzCamera();
    (*render)->camera.FOV = DEFAULT_FOV;
    GzCoord origin = { 0, 0, 0 };
    memcpy((*render)->camera.lookat, origin, sizeof(GzCoord));
    GzCoord cam_pos = { DEFAULT_IM_X, DEFAULT_IM_Y, DEFAULT_IM_Z };

    memcpy((*render)->camera.position, cam_pos, sizeof(GzCoord));
    GzCoord up = { 0, 1.0f, 0 };
    memcpy((*render)->camera.worldup, up, sizeof(GzCoord));

    (*render)->matlevel = 0;

    GzMatrix m_sp =
    {
        display->xres / 2.0f, 0, 0, display->xres / 2.0f,
        0, -display->yres / 2.0f, 0, display->yres / 2.0f,
        0, 0, float(INT_MAX), 0,
        0, 0, 0, 1.0f
    };

    memcpy((*render)->Xsp, m_sp, sizeof(GzMatrix));

    (*render)->numlights = 0;

    return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render) {
    /*
    -free all renderer resources
    */
    free(render);
    return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render) {
    /*
    - setup for start of each frame - init frame buffer color,alpha,z
    - compute Xiw and projection xform Xpi from camera definition
    - init Ximage - put Xsp at base of stack, push on Xpi and Xiw
    - now stack contains Xsw and app can push model Xforms when needed
    */

    float d_inverse = tan(degToRad(render->camera.FOV) / 2.0f);

    GzCoord cam_Z = {
        render->camera.lookat[X] - render->camera.position[X],
        render->camera.lookat[Y] - render->camera.position[Y],
        render->camera.lookat[Z] - render->camera.position[Z],
    };
    normalizeGzCoord(cam_Z);


    float updotz = dotProduct(render->camera.worldup, cam_Z);
    GzCoord cam_Y = {
        render->camera.worldup[X] - updotz * cam_Z[X],
        render->camera.worldup[Y] - updotz * cam_Z[Y],
        render->camera.worldup[Z] - updotz * cam_Z[Z],
    };
    normalizeGzCoord(cam_Y);

    GzCoord* cam_X = (GzCoord*)malloc(sizeof(GzCoord));
    crossProduct(cam_Y, cam_Z, cam_X);

    float minus_cx = -dotProduct(*cam_X, render->camera.position);
    float minus_cy = -dotProduct(cam_Y, render->camera.position);
    float minus_cz = -dotProduct(cam_Z, render->camera.position);

    GzMatrix xiw =
    {
        (*cam_X)[X], (*cam_X)[Y], (*cam_X)[Z], minus_cx,
        cam_Y[X], cam_Y[Y], cam_Y[Z], minus_cy,
        cam_Z[X], cam_Z[Y], cam_Z[Z], minus_cz,
        0, 0, 0, 1.0f
    };


    GzMatrix xpi =
    {
        1.0f, 0, 0, 0,
        0, 1.0f, 0, 0,
        0, 0, d_inverse, 0,
        0, 0, d_inverse, 1.0f
    };

    free(cam_X);

    memcpy(render->camera.Xiw, xiw, sizeof(GzMatrix));
    memcpy(render->camera.Xpi, xpi, sizeof(GzMatrix));

    GzInitDisplay(render->display);

    GzMatrix i = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
    memcpy(render->Ximage[render->matlevel], i, sizeof(GzMatrix));
    memcpy(render->Xnorm[render->matlevel], i, sizeof(GzMatrix));

    if (!GzPushMatrix(render, render->Xsp))
        if (!GzPushMatrix(render, render->camera.Xpi))
            if (!GzPushMatrix(render, render->camera.Xiw))
                return GZ_SUCCESS;


    return GZ_FAILURE;
}

int GzPutCamera(GzRender *render, GzCamera *camera) {
    /*
    - overwrite renderer camera structure with new camera definition
    */
    render->camera = *camera;
    // make sure new data being fed is correct
    normalizeGzCoord(render->camera.worldup);
    return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix) {
    /*
    - push a matrix onto the Ximage stack
    - check for stack overflow
    */
    if (render->matlevel < MATLEVELS) {

        GzMatrix m = {};
        memcpy(m, matrix, sizeof(GzMatrix));

        GzMatrix* top = (GzMatrix*)malloc(sizeof(GzMatrix));

        MMmultiply(render->Ximage[render->matlevel], m, top);
        memcpy(render->Ximage[render->matlevel+1], top, sizeof(GzMatrix));


        // skip Xsp && Xpi for normals (i.e. push identity)
        if (memcmp(m, render->Xsp, sizeof(GzMatrix)) == 0 || memcmp(m, render->camera.Xpi, sizeof(GzMatrix)) == 0) {
            GzMatrix i = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
            memcpy(render->Xnorm[render->matlevel+1], i, sizeof(GzMatrix));
            render->matlevel++;
            return GZ_SUCCESS;
        }

        GzMatrix* norm_top = (GzMatrix*)malloc(sizeof(GzMatrix));

        //make "m" a unitary rotation matrix

        // scale factor
        float k = std::sqrt(m[0][0] * m[0][0] + m[0][1] * m[0][1] + m[0][2] * m[0][2]);


        GzMatrix r = {
            m[0][0] / k, m[0][1] / k, m[0][2] / k, 0,
            m[1][0] / k, m[1][1] / k, m[1][2] / k, 0,
            m[2][0] / k, m[2][1] / k, m[2][2] / k, 0,
            0, 0, 0, 1,
        };

        MMmultiply(render->Xnorm[render->matlevel], r, norm_top);
        memcpy(render->Xnorm[render->matlevel+1], norm_top, sizeof(GzMatrix));

        render->matlevel++;

        free(norm_top);
        free(top);

        return GZ_SUCCESS;
    }
    else
        return GZ_FAILURE;
}

int GzPopMatrix(GzRender *render) {
    /*
    - pop a matrix off the Ximage stack
    - check for stack underflow
    */
    if (render->matlevel > 0) {
        render->matlevel--;
        return GZ_SUCCESS;
    }
    else
        return GZ_FAILURE;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList,
    GzPointer	*valueList) /* void** valuelist */
{
    /*
    - set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
    - later set shaders, interpolaters, texture maps, and lights
    */

    if (numAttributes == 0)
        return GZ_FAILURE;
    for (int i = 0; i < numAttributes; i++) {
        GzToken token = nameList[i];
        GzPointer val = valueList[i];
        switch (token) {
        case GZ_NULL_TOKEN:
            break;
        case GZ_RGB_COLOR: {
            GzColor* c = static_cast<GzColor*>(val);
            render->flatcolor[0] = (*c)[0];
            render->flatcolor[1] = (*c)[1];
            render->flatcolor[2] = (*c)[2];
            break;
        }
        case GZ_DIRECTIONAL_LIGHT:
            render->lights[render->numlights++] = *(GzLight*)val;
            // normalize in case given dir isn't
            normalizeGzCoord(render->lights[render->numlights-1].direction);
            break;
        case GZ_AMBIENT_LIGHT:
            render->ambientlight = *(GzLight*)val;
            break;
        case GZ_DIFFUSE_COEFFICIENT:
            render->Kd[0] = (*(GzColor*)val)[0];
            render->Kd[1] = (*(GzColor*)val)[1];
            render->Kd[2] = (*(GzColor*)val)[2];
            break;
        case GZ_AMBIENT_COEFFICIENT:
            render->Ka[0] = (*(GzColor*)val)[0];
            render->Ka[1] = (*(GzColor*)val)[1];
            render->Ka[2] = (*(GzColor*)val)[2];
            break;
        case GZ_SPECULAR_COEFFICIENT:
            render->Ks[0] = (*(GzColor*)val)[0];
            render->Ks[1] = (*(GzColor*)val)[1];
            render->Ks[2] = (*(GzColor*)val)[2];
            break;
        case GZ_DISTRIBUTION_COEFFICIENT:
            render->spec = *(float*)val;
            break;
        case GZ_INTERPOLATE:
            render->interp_mode = *(int*)val;
            break;
        case GZ_TEXTURE_MAP:
            render->tex_fun = (GzTexture)val;
            break;
        default:
            return GZ_FAILURE;
        }
    }

    return GZ_SUCCESS;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
/* numParts : how many names and values */
{
    /*
    - pass in a triangle description with tokens and values corresponding to
    GZ_POSITION:3 vert positions in model space
    - Xform positions of verts using matrix on top of stack
    - Clip - just discard any triangle with any vert(s) behind view plane
    - optional: test for triangles with all three verts off-screen (trivial frustum cull)
    - invoke triangle rasterizer
    */

    if (numParts == 0)
        return GZ_FAILURE;

    GzCoord *tri0 = (GzCoord *)malloc(sizeof(GzCoord));
    GzCoord *tri1 = (GzCoord *)malloc(sizeof(GzCoord));
    GzCoord *tri2 = (GzCoord *)malloc(sizeof(GzCoord));

    GzCoord *n0 = (GzCoord *)malloc(sizeof(GzCoord));
    GzCoord *n1 = (GzCoord *)malloc(sizeof(GzCoord));
    GzCoord *n2 = (GzCoord *)malloc(sizeof(GzCoord));

    for (int i = 0; i < numParts; i++) {
        GzToken token = nameList[i];
        GzPointer val = valueList[i];
        switch (token) {
        case GZ_NULL_TOKEN:
            break;
        case GZ_POSITION: {

            const GzCoord* tri[3];
            tri[0] = &((GzCoord*)val)[0];
            tri[1] = &((GzCoord*)val)[1];
            tri[2] = &((GzCoord*)val)[2];

            MVmultiply(render->Ximage[render->matlevel], *tri[0], tri0);
            MVmultiply(render->Ximage[render->matlevel], *tri[1], tri1);
            MVmultiply(render->Ximage[render->matlevel], *tri[2], tri2);

            break;
        }
        case GZ_NORMAL: {
            // assume normals are normalized when input
            const GzCoord* normals[3];
            normals[0] = &((GzCoord*)val)[0];
            normals[1] = &((GzCoord*)val)[1];
            normals[2] = &((GzCoord*)val)[2];

            MVmultiply(render->Xnorm[render->matlevel], *normals[0], n0);
            MVmultiply(render->Xnorm[render->matlevel], *normals[1], n1);
            MVmultiply(render->Xnorm[render->matlevel], *normals[2], n2);

            break;
        }
        default:
            return GZ_FAILURE;
        }
    }

    const GzCoord* tri_w[3] = {};
    tri_w[0] = tri0;
    tri_w[1] = tri1;
    tri_w[2] = tri2;

    // cull negative z
    if (!((*tri_w[0])[2] < 0 && (*tri_w[1])[2] < 0 && (*tri_w[2])[2] < 0)) {
        GzCoord* n_w[3] = {};
        n_w[0] = n0;
        n_w[1] = n1;
        n_w[2] = n2;

        sortVerticesCW(tri_w, n_w);
        

    }
    return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
    return(short)((int)(color * ((1 << 12) - 1)));
}

/* HW 2 */

int findGeneralLineEq(const GzCoord* tail, const GzCoord* head, GzCoord &out) {
    // solve for A, B, C in Ax+By+C=0

    out[A] = (*head)[Y] - (*tail)[Y];
    out[B] = (*tail)[X] - (*head)[X];
    out[C] = ((*head)[X] - (*tail)[X])* (*tail)[Y] - ((*head)[Y] - (*tail)[Y])* (*tail)[X];
    return GZ_SUCCESS;
}

int findGeneralPlaneEq(const GzCoord* (&tri)[3], GzCoord &out, float &d) {
    // solve for A, B, C, D in Ax+By+Cz+D=0

    GzCoord e1 = {
        (*tri[0])[X] - (*tri[1])[X],
        (*tri[0])[Y] - (*tri[1])[Y],
        (*tri[0])[Z] - (*tri[1])[Z]
    };
    GzCoord e2 = {
        (*tri[0])[X] - (*tri[2])[X],
        (*tri[0])[Y] - (*tri[2])[Y],
        (*tri[0])[Z] - (*tri[2])[Z] };

    // cross product, gives us A,B,C
    out[A] = e1[Y] * e2[Z] - e1[Z] * e2[Y];
    out[B] = e1[Z] * e2[X] - e1[X] * e2[Z];
    out[C] = e1[X] * e2[Y] - e1[Y] * e2[X];

    // solve for D by plugging in a vert
    d = -out[A] * (*tri[0])[X] - out[B] * (*tri[0])[Y] - out[C] * (*tri[0])[Z];

    return GZ_SUCCESS;
}

int sortVerticesCW(const GzCoord* (&tri)[3], GzCoord* (&normals)[3]) {

    // determine top vert
    // sort by Y ascending
    // if Y-tie: break by putting leftmost first
    // which results in a CW arrangement later

    if ((*tri[1])[Y] < (*tri[0])[Y] ||
        ((*tri[1])[Y] == (*tri[0])[Y] && (*tri[1])[X] < (*tri[0])[X])) {
        std::swap(tri[0], tri[1]);
        std::swap(normals[0], normals[1]);
    }
        

        if ((*tri[2])[Y] < (*tri[1])[Y] ||
            ((*tri[2])[Y] == (*tri[1])[Y] && (*tri[2])[X] < (*tri[1])[X])) {
            std::swap(tri[1], tri[2]);
            std::swap(normals[1], normals[2]);
        }
        

        if ((*tri[1])[Y] < (*tri[0])[Y] ||
            ((*tri[1])[Y] == (*tri[0])[Y] && (*tri[1])[X] < (*tri[0])[X])) {
            std::swap(tri[0], tri[1]);
            std::swap(normals[0], normals[1]);
        }
        

    // find TOP-BOT line

    GzCoord abc = { 0, 0, 0 };
    findGeneralLineEq(tri[0], tri[2], abc);

    // determine if mid is right or left

    // find x-coord of point on bot-top line at mid's Y
    float x = (-abc[C] - abc[B] * (*tri[1])[Y]) / abc[A];

    // topmost vertice first
    if (x < (*tri[1])[X]) {
        // then mid is to the right of bot-top line
        // so put mid second for final order top,mid,bot
        std::swap(tri[0], tri[2]);
        std::swap(normals[0], normals[2]);
    }
    else {
        // then mid is to the left of bot-top line
        // put mid last for final order top,bot,mid
        std::swap(tri[0], tri[2]);
        std::swap(normals[0], normals[2]);

        std::swap(tri[1], tri[2]);
        std::swap(normals[1], normals[2]);
    }
    return GZ_SUCCESS;
}

void computeColorAtPt(GzRender *render, const GzCoord &pt, GzCoord &normal, GzCoord *color) {

    GzCoord E;
    E[X] = pt[X] - render->camera.position[X];
    E[Y] = pt[Y] - render->camera.position[Y];
    E[Z] = pt[Z] - render->camera.position[Z];
    normalizeGzCoord(E);

    GzColor sum_spec = { 0, 0, 0 };
    GzColor sum_diff = { 0, 0, 0 };
    for (int i = 0; i < render->numlights; i++) {

        //l goes from vert to light, so flip sign of light direction
        GzCoord L;
        L[X] = -render->lights[i].direction[X];
        L[Y] = -render->lights[i].direction[Y];
        L[Z] = -render->lights[i].direction[Z];

        float LdotN = dotProduct(L, normal);
        float NdotE = dotProduct(E, normal);

        //check signs to indicate direction
        if (LdotN < 0 && NdotE < 0) {
            //both negative -> flip N and recalculate L dot N
            normal[X] = -normal[X];
            normal[Y] = -normal[Y];
            normal[Z] = -normal[Z];
            LdotN = dotProduct(L, normal);
        }
        else if ((LdotN < 0 && NdotE > 0) || (LdotN > 0 && NdotE < 0)) {
            //opposite signs: skip it
            continue;
        }
        // both positive: keep going

        //spec

        // R = 2 * (L dot N)N - L
        GzCoord R;
        float r_scale = 2 * LdotN;
        R[X] = r_scale * normal[X] - L[X];
        R[Y] = r_scale * normal[Y] - L[Y];
        R[Z] = r_scale * normal[Z] - L[Z];

        normalizeGzCoord(R);

        float RdotE = dotProduct(R, E);
        // clamp RdotE to zero
        if (RdotE < 0.0f)
            RdotE = 0.0f;

        float RdotEpowS = std::powf(RdotE, render->spec);
        sum_spec[X] += render->lights[i].color[X] * RdotEpowS;
        sum_spec[Y] += render->lights[i].color[Y] * RdotEpowS;
        sum_spec[Z] += render->lights[i].color[Z] * RdotEpowS;

        // diff

        sum_diff[X] += render->lights[i].color[X] * LdotN;
        sum_diff[Y] += render->lights[i].color[Y] * LdotN;
        sum_diff[Z] += render->lights[i].color[Z] * LdotN;
    }

    sum_spec[X] *= render->Ks[X];
    sum_spec[Y] *= render->Ks[Y];
    sum_spec[Z] *= render->Ks[Z];

    sum_diff[X] *= render->Kd[X];
    sum_diff[Y] *= render->Kd[Y];
    sum_diff[Z] *= render->Kd[Z];

    (*color)[X] = render->Ka[X] * render->ambientlight.color[X] + sum_diff[X] + sum_spec[X];
    (*color)[Y] = render->Ka[Y] * render->ambientlight.color[Y] + sum_diff[Y] + sum_spec[Y];
    (*color)[Z] = render->Ka[Z] * render->ambientlight.color[Z] + sum_diff[Z] + sum_spec[Z];
}

int lee(GzRender *render, const GzCoord* (&tri)[3], GzCoord* (&normals)[3], GzTextureIndex* (&tex)[3]) {

    //perspective correction on vertex texture index
    for (int triIndex = 0; triIndex < 3; triIndex++) {
        float vPrimeZ = (*tri[triIndex])[Z] / (float(INT_MAX) - (*tri[triIndex])[Z]);
        for (int i = 0; i < 2; i++) {  
            (*tex[triIndex])[i] = (*tex[triIndex])[i] / (vPrimeZ + 1.0f);
        }
    }

    GzCoord* tri0_color = (GzCoord*)malloc(sizeof(GzCoord));
    GzCoord* tri1_color = (GzCoord*)malloc(sizeof(GzCoord));
    GzCoord* tri2_color = (GzCoord*)malloc(sizeof(GzCoord));

    GzCoord* tri_colors[3] = { tri0_color, tri1_color, tri2_color };

    if (render->interp_mode == GZ_COLOR) {

        // if GZ_COLOR, do K's later: make them =1
        GzColor one = { 1, 1, 1 };
        memcpy(render->Ka, one, sizeof(GzColor));
        memcpy(render->Kd, one, sizeof(GzColor));
        memcpy(render->Ks, one, sizeof(GzColor));

        for (int triIndex = 0; triIndex < 3; triIndex++) {
            computeColorAtPt(render, *tri[triIndex], *normals[triIndex], tri_colors[triIndex]);
        }
    }

    // find general eq. for tri edges
    GzCoord edges[3] = { 0, 0, 0 };
    findGeneralLineEq(tri[0], tri[1], edges[0]);
    findGeneralLineEq(tri[1], tri[2], edges[1]);
    findGeneralLineEq(tri[2], tri[0], edges[2]);

    // find bounds for acceptable pixels
    int maxX = (int)max((*tri[0])[X], max((*tri[1])[X], (*tri[2])[X]));
    int minX = (int)min((*tri[0])[X], min((*tri[1])[X], (*tri[2])[X]));
    int maxY = (int)max((*tri[0])[Y], max((*tri[1])[Y], (*tri[2])[Y]));
    int minY = (int)min((*tri[0])[Y], min((*tri[1])[Y], (*tri[2])[Y]));

    // find plane of triangle
    GzCoord abc = { 0, 0, 0 };
    float d;
    findGeneralPlaneEq(tri, abc, d);

    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            int in_triangle = 0;
            for (int i = 0; i < 3; i++) {
                float result =
                    (edges[i])[A] * x +
                    (edges[i])[B] * y +
                    (edges[i])[C];
                // Ax + By + C < 0 not in tri
                if (result < 0)
                    break;
                // = 0 means on the edge
                if (result == 0)
                    // throw out if edge has - slope
                    if (-(edges[i])[A] / (edges[i])[B] < 0)
                        break;
                // else point is on "correct" side of edge
                in_triangle++;
            }
            // if point inside for all edges
            if (in_triangle == 3) {

                // draw pixel if z < zbuf
                int zbuf;
                GzGetDisplay(render->display, x, y, 0, 0, 0, 0, &zbuf);

                // interpolate z
                GzDepth z = (GzDepth)std::round((-abc[A] * x - abc[B] * y - d) / abc[C]);

                if (z < zbuf) {

                    if (render->interp_mode == GZ_FLAT) {
                        GzIntensity r = ctoi(render->flatcolor[RED]);
                        GzIntensity g = ctoi(render->flatcolor[GREEN]);
                        GzIntensity b = ctoi(render->flatcolor[BLUE]);

                        GzIntensity a = 1;

                        GzPutDisplay(render->display, x, y, r, g, b, a, z);
                    }
                    else if (render->interp_mode == GZ_COLOR) {

                        GzCoord point = { x, y, z };
                        GzCoord* bary = (GzCoord*)malloc(sizeof(GzCoord));

                        (*bary)[X] = 0; (*bary)[Y] = 0; (*bary)[Z] = 0;

                        //BARYCENTRIC INTERPOLATION OF COLOUR VALUES AT VERTS
                        barycentricCoords(*tri[0], *tri[1], *tri[2], point, bary);

                        GzTextureIndex pt_uv;
                        pt_uv[X] = (*bary)[X] * (*tex[X])[X] + (*bary)[Y] * (*tex[Y])[X] + (*bary)[Z] * (*tex[Z])[X];
                        pt_uv[Y] = (*bary)[X] * (*tex[X])[Y] + (*bary)[Y] * (*tex[Y])[Y] + (*bary)[Z] * (*tex[Z])[Y];

                        //perspective-correction on pixel tex
                        float vPrimeZ = point[Z] / (float(INT_MAX) - point[Z]);
                        pt_uv[X] = pt_uv[X] * (vPrimeZ + 1.0f);
                        pt_uv[Y] = pt_uv[Y] * (vPrimeZ + 1.0f);

                        GzCoord tex_color = { 0, 0, 0 };
                        if (render->tex_fun(pt_uv[X], pt_uv[Y], tex_color) == GZ_FAILURE) {
                            return GZ_FAILURE;
                        }

                        GzIntensity r = ctoi(tex_color[RED] * (*bary)[X] * (*tri0_color)[RED] + tex_color[RED] * (*bary)[1] * (*tri1_color)[RED] + tex_color[RED] * (*bary)[2] * (*tri2_color)[RED]);
                        GzIntensity g = ctoi(tex_color[GREEN] * (*bary)[X] * (*tri0_color)[GREEN] + tex_color[GREEN] * (*bary)[1] * (*tri1_color)[GREEN] + tex_color[GREEN] * (*bary)[2] * (*tri2_color)[GREEN]);
                        GzIntensity b = ctoi(tex_color[BLUE] * (*bary)[X] * (*tri0_color)[BLUE] + tex_color[BLUE] * (*bary)[1] * (*tri1_color)[BLUE] + tex_color[BLUE] * (*bary)[2] * (*tri2_color)[BLUE]);

                        GzIntensity a = 1;

                        GzPutDisplay(render->display, x, y, r, g, b, a, z);
                        free(bary);
                    }
                    else if (render->interp_mode == GZ_NORMALS) {

                        GzCoord point = { x, y, z };
                        GzCoord* bary = (GzCoord*)malloc(sizeof(GzCoord));

                        (*bary)[X] = 0; (*bary)[Y] = 0; (*bary)[Z] = 0;

                        //BARYCENTRIC INTERPOLATION OF NORMAL
                        barycentricCoords(*tri[0], *tri[1], *tri[2], point, bary);

                        GzCoord pt_normal;
                        pt_normal[X] = (*bary)[X] * (*normals[X])[X] + (*bary)[Y] * (*normals[Y])[X] + (*bary)[Z] * (*normals[Z])[X];
                        pt_normal[Y] = (*bary)[X] * (*normals[X])[Y] + (*bary)[Y] * (*normals[Y])[Y] + (*bary)[Z] * (*normals[Z])[Y];
                        pt_normal[Z] = (*bary)[X] * (*normals[X])[Z] + (*bary)[Y] * (*normals[Y])[Z] + (*bary)[Z] * (*normals[Z])[Z];
                        //normalize in case of loss of precision
                        normalizeGzCoord(pt_normal);

                        GzTextureIndex pt_uv;
                        pt_uv[X] = (*bary)[X] * (*tex[X])[X] + (*bary)[Y] * (*tex[Y])[X] + (*bary)[Z] * (*tex[Z])[X];
                        pt_uv[Y] = (*bary)[X] * (*tex[X])[Y] + (*bary)[Y] * (*tex[Y])[Y] + (*bary)[Z] * (*tex[Z])[Y];

                        //perspective-correction on pixel tex
                        float vPrimeZ = point[Z] / (float(INT_MAX) - point[Z]);
                        pt_uv[X] = pt_uv[X] * (vPrimeZ + 1.0f);
                        pt_uv[Y] = pt_uv[Y] * (vPrimeZ + 1.0f);

                        GzCoord tex_color = { 0, 0, 0 };
                        if (render->tex_fun(pt_uv[X], pt_uv[Y], tex_color) == GZ_FAILURE) {
                            return GZ_FAILURE;
                        }
                        memcpy(render->Ka, tex_color, sizeof(GzColor));
                        memcpy(render->Kd, tex_color, sizeof(GzColor));

                        GzCoord* pt_color = (GzCoord*)malloc(sizeof(GzCoord));
                        computeColorAtPt(render, point, pt_normal, pt_color);

                        GzIntensity r = ctoi((*pt_color)[0]);
                        GzIntensity g = ctoi((*pt_color)[1]);
                        GzIntensity b = ctoi((*pt_color)[2]);

                        GzIntensity a = 1;
                        
                        GzPutDisplay(render->display, x, y, r, g, b, a, z);
                        free(bary);
                        free(pt_color);
                    }
                }
            }
        }
    }
    free(tri0_color);
    free(tri1_color);
    free(tri2_color);
    return GZ_SUCCESS;
}
//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool rayTriangleIntersect(
	GzCoord &orig, GzCoord &dir,
	GzCoord &v0, const GzCoord &v1, const GzCoord &v2,
	float &out)
{
	GzCoord e1, e2;  //Edge1, Edge2
	GzCoord P, Q, T;
	float det, inv_det, u, v;
	float t;

	//Find vectors for two edges sharing V1
	e1[0] = v1[0] - v0[0];
	e1[1] = v1[1] - v0[1];
	e1[2] = v1[2] - v0[2];

	e2[0] = v2[0] - v0[0];
	e2[1] = v2[1] - v0[1];
	e2[2] = v2[2] - v0[2];

	//Begin calculating determinant - also used to calculate u parameter
	crossProduct(dir, e2, &P);
	//CROSS(P, D, e2);

	//if determinant is near zero, ray lies in plane of triangle
	det = dotProduct(e1, P);

	//NOT CULLING
	if (det > -0.000001 && det < 0.000001) return false;
	inv_det = 1.f / det;

	//calculate distance from V1 to ray origin
	T[0] = orig[0] - v0[0];
	T[1] = orig[1] - v0[1];
	T[2] = orig[2] - v0[2];


	//Calculate u parameter and test bound
	u = dotProduct(T, P) * inv_det;
	//The intersection lies outside of the triangle
	if (u < 0.f || u > 1.f) return false;

	//Prepare to test v parameter
	//CROSS(Q, T, e1);
	crossProduct(T, e1, &Q);

	//Calculate V parameter and test bound
	v = dotProduct(dir, Q) * inv_det;
	//The intersection lies outside of the triangle
	if (v < 0.f || u + v  > 1.f) return false;

	t = dotProduct(e2, Q) * inv_det;

	if (t > 0.000001) { //ray intersection
		out = t;
		return true;
	}

	// No hit, no win
	return false;
}

int raycast_render(GzRender *render, GzWorldSpaceTriangles *tris)
{
	//Loop through all pixels:
	for (int i = 0; i < render->display->xres; i++)
	{
		for (int j = 0; j < render->display->yres; j++)
		{
			//Create Ray:
			GzRay ray;
			ray.position[0] = render->camera.position[0];
			ray.position[1] = render->camera.position[1];
			ray.position[2] = render->camera.position[2];

			float d = 50;//1 / (tan(degToRad(render->camera.FOV) / 2.0f));
			float W = 2;//2 * d*(tan(degToRad(render->camera.FOV) / 2.0f)); //Horizontal Field of View
			float H = 2;// 2 * d*(tan(degToRad(render->camera.FOV) / 2.0f)); //Assuming same field of view in vertical

			GzCoord cam_Z = {
				render->camera.lookat[X] - render->camera.position[X],
				render->camera.lookat[Y] - render->camera.position[Y],
				render->camera.lookat[Z] - render->camera.position[Z],
			};
			normalizeGzCoord(cam_Z);

			float updotz = dotProduct(render->camera.worldup, cam_Z);
			GzCoord cam_Y = {
				render->camera.worldup[X] - updotz * cam_Z[X],
				render->camera.worldup[Y] - updotz * cam_Z[Y],
				render->camera.worldup[Z] - updotz * cam_Z[Z],
			};
			normalizeGzCoord(cam_Y);

			GzCoord cam_X;
			crossProduct(cam_Y, cam_Z, &cam_X);

			GzCoord Center = {
				render->camera.position[0] + (cam_Z[X] * d),
				render->camera.position[1] + (cam_Z[Y] * d),
				render->camera.position[2] + (cam_Z[Z] * d)
			};

			GzCoord BottomLeft = {
				Center[X] - ((cam_X)[X] * (W / 2)) + ((cam_Y)[X] * (H / 2)),
				Center[Y] - ((cam_X)[Y] * (W / 2)) + ((cam_Y)[Y] * (H / 2)),
				Center[Z] - ((cam_X)[Z] * (W / 2)) + ((cam_Y)[Z] * (H / 2))
			};

			GzCoord pixelPos = {
				BottomLeft[X] + ((cam_X)[X] * (i + .5)*(W / render->display->xres)) - ((cam_Y)[X] * (j + .5)*(H / render->display->yres)),
				BottomLeft[Y] + ((cam_X)[Y] * (i + .5)*(W / render->display->xres)) - ((cam_Y)[Y] * (j + .5)*(H / render->display->yres)),
				BottomLeft[Z] + ((cam_X)[Z] * (i + .5)*(W / render->display->xres)) - ((cam_Y)[Z] * (j + .5)*(H / render->display->yres))

			};

			//Creating eye rays
			//http://web.cse.ohio-state.edu/~hwshen/681/Site/Slides_files/basic_algo.pdf

			pixelPos[0] = pixelPos[0] - ray.position[0];
			pixelPos[1] = pixelPos[1] - ray.position[1];
			pixelPos[2] = pixelPos[2] - ray.position[2];

			normalizeGzCoord(pixelPos);

			ray.direction[0] = pixelPos[0];
			ray.direction[1] = pixelPos[1];
			ray.direction[2] = pixelPos[2];
			normalizeGzCoord(ray.direction);

			//for each triangle
			for (int k = 0; k < tris->tris.size() - 1; k++)
			{
				GzCoord triA;
				triA[0] = tris->tris[k]->vertices[0]->pos[0];//(tri0[0])[0];
				triA[1] = tris->tris[k]->vertices[0]->pos[1];
				triA[2] = tris->tris[k]->vertices[0]->pos[2];

				GzCoord triB;
				triB[0] = tris->tris[k]->vertices[1]->pos[0];
				triB[1] = tris->tris[k]->vertices[1]->pos[1];
				triB[2] = tris->tris[k]->vertices[1]->pos[2];

				GzCoord triC;
				triC[0] = tris->tris[k]->vertices[2]->pos[0];
				triC[1] = tris->tris[k]->vertices[2]->pos[1];
				triC[2] = tris->tris[k]->vertices[2]->pos[2];

				//double t = ray_intersect(ray, triA, triB, triC);
				float t_test = 0;
				bool test = rayTriangleIntersect(ray.position, ray.direction, triA, triB, triC, t_test);

				//if (t != 0)
				if (test)
				{
					//This is where I should check for ray intersection
					//Check if the ray intersects a surface,
					//if so fetch the colour of the surface and assign a new direction to the ray from the point of intersection
					
					//Create a recursive method to reflect ray and fetch colour
					float* newColor =  TracePath(ray, 0);		
					
					GzPutDisplay(render->display, i, j, 4095.0f, 0.0f, 0.0f, 1.0f, 1.0f);
				}
			}
			

		}
	}

	return 0;
}


float* TracePath(GzRay ray, int depth)
{
	if (depth == MAXREFLECTANCE)
	{
		GzCoord black; 
		black[0] = 0.0f; //r g b
		black[1] = 0.0f;
		black[2] = 0.0f;

		return black; 
	}

	//Fetch the colour of the ray position 
	GzColor c; // ray.colour(some way to assign the colour of the object)

	GzRay newray;
	newray.position[0] = ray.position[0]; // newray origin
	newray.position[1] = ray.position[1];
	newray.position[2] = ray.position[2];

	//newray.direction = pick a random direction 


	// Compute the BRDF for this ray (assuming Lambertian reflection)
	
	GzCoord normalWhereObjectWasHit; //create a normal and assign
	float cos_theta = dotProduct(newray.direction, normalWhereObjectWasHit);
//	float* BRDF = 2 * m.reflectance * cos_theta;*/
	float* reflected = TracePath(newray, depth + 1);

	// Apply the Rendering Equation here.
	//return emittance + (BRDF * reflected);*/
}
