// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>
#include <algorithm>

#include <omp.h> 
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <GL/glut.h>

#include "GeoHashGrid.hpp"
#include "GeoCoordinate.hpp"
#include "SphKernel.hpp"
#include "SphCubicSpline.hpp"
#include "SphGaussian.hpp"
#include "SphQuintic.hpp"
#include "SphSpiky.hpp"

using namespace std;

using namespace ENigMA::geometry;
using namespace ENigMA::sph;

#define kPi 3.1415926535f

const int NX = 120;         // x
const int NY = 120;         // y

int NUMPARTICLES = NX * NY;

const float dt = 1E-2f;      // time step
const float h = 10.0f;
const float m0 = 10.0f;
const float rho0 = 10.0f;
const float visc0 = 1.0f;

struct CParticle
{
    CGeoCoordinate<float> pos;
    CGeoVector<float> vel;

    float mass;
    float p;
    float rho;
    float visc;

    float factor;
};

struct CRectangle
{
    CGeoCoordinate<float> p1, p2;
};

CRectangle walls;

std::vector<CParticle> p;

int winwidth = 480;
int winheight = 480;

float scale[6][4];

int frames = 0;
int t0 = 0, te;

CGeoHashGrid<float> aHashGrid;
CSphCubicSpline<float> aCubicSplineKernel(2);
CSphGaussian<float> aGaussianKernel(2);
CSphQuintic<float> aQuinticKernel(2);
CSphSpiky<float> aSpikyKernel(2);

void getColor(float v, float& r, float& g, float& b)
{

    if (v <= scale[0][0])
    {
        r = scale[0][1];
        g = scale[0][2];
        b = scale[0][3];
    }
    else if (v >= scale[4][0])
    {
        r = scale[4][1];
        g = scale[4][2];
        b = scale[4][3];
    }
    else
    {

        int i = 0;

        for (i = 1; i <= 4; ++i)
            if (scale[i - 1][0] <= v && v < scale[i][0])
                break;

        // linear interpolation
        float w = (v - scale[i - 1][0]) / (scale[i][0] - scale[i - 1][0]);

        r = (1 - w)*scale[i - 1][1] + w*scale[i][1];
        g = (1 - w)*scale[i - 1][2] + w*scale[i][2];
        b = (1 - w)*scale[i - 1][3] + w*scale[i][3];

    }
}

void buildHashGrid()
{

    aHashGrid.reset();

    // Build hash grid
    for (int pi = 0; pi < NUMPARTICLES; ++pi)
        aHashGrid.addGeometricObject(pi, p[pi].pos);

    aHashGrid.build();

}

float upperRectangleFunction()
{
    return walls.p1.y();
}

float lowerRectangleFunction()
{
    return walls.p2.y();
}

float upperCircleFunction(CGeoCoordinate<float>& center, float r, float x)
{
    return center.y() - sqrt((r * r) - pow((x - center.x()), 2));
}

float lowerCircleFunction(CGeoCoordinate<float>& center, float r, float x)
{
    return center.y() + sqrt((r * r) - pow((x - center.x()), 2));
}

float intersectionArea(CGeoCoordinate<float>& center, float r)
{

    float area = 0;

    //A variable storing the nearest horizontal edge of the rectangle. 
    float nearestRectangleEdge = 0;

    //Determine what is nearer to the circle center - the rectangle top edge or the rectangle bottom edge
    if (fabs(center.y() - walls.p2.y()) > fabs(center.y() - walls.p1.y()))
        nearestRectangleEdge = walls.p2.y();
    else
        nearestRectangleEdge = walls.p1.y();

    //The bounds of our integration
    float leftBound = 0.0;
    float rightBound = 0.0;

    if (center.y() >= walls.p1.y() && center.y() <= walls.p2.y())
    {
        //Take care if the circle's center lies within the rectangle. 
        leftBound = std::max(-r + center.x(), walls.p1.x());
        rightBound = std::min(r + center.x(), walls.p2.x());
    }
    else if (r >= fabs(nearestRectangleEdge - center.y()))
    {
        //If the circle's center lies outside of the rectangle, we can choose optimal bounds.
        leftBound = std::max(static_cast<float>(-sqrt(r * r - fabs(pow(nearestRectangleEdge - center.y(), 2))) + center.y()), walls.p1.x());
        rightBound = std::min(static_cast<float>(sqrt(r * r - fabs(pow(nearestRectangleEdge - center.y(), 2))) + center.y()), walls.p2.x());
    }

    float upperBound;
    float lowerBound;

    float resolution = 0.01f;

    //Loop trough the intersection area and sum up the area
    for (float i = leftBound + resolution; i <= rightBound; i += resolution)
    {
        upperBound = std::max(upperRectangleFunction(), upperCircleFunction(center, r, i - resolution / 2));
        lowerBound = std::min(lowerRectangleFunction(), lowerCircleFunction(center, r, i - resolution / 2));

        area += (lowerBound - upperBound) * resolution;
    }

    return area;

}

void calculateDensityAndPressure()
{

    float B = 1000.0;

#pragma omp parallel for num_threads(4)
    for (int pi = 0; pi < NUMPARTICLES; ++pi)
    {

        p[pi].rho = 0.0;

        std::vector<Integer> sParticles;
        aHashGrid.find(sParticles, p[pi].pos, h);

        for (Integer n = 0; n < sParticles.size(); ++n)
        {

            Integer pj = sParticles[n];

            CGeoVector<float> r = p[pj].pos - p[pi].pos;

            //p[pi].rho += p[pj].mass * aCubicSplineKernel.W(r, h);
            //p[pi].rho += p[pj].mass * aGaussianKernel.W(r, h);
            //p[pi].rho += p[pj].mass * aQuinticKernel.W(r, h);
            p[pi].rho += p[pj].mass * aSpikyKernel.W(r, h);

        }

        float area = intersectionArea(p[pi].pos, h);

        if (area > 0.0)
        {
            p[pi].factor = 3.14 * h * h / area;
            p[pi].rho *= p[pi].factor;
        }

        if (p[pi].rho - rho0 > 0)
            p[pi].p = B * (float)(pow(p[pi].rho / rho0, 7.0) - 1.0);
        else
            p[pi].p = 0.0;

    }

}

void moveParticles()
{

    float dpx, dpy;
    float dviscx, dviscy;

    float gx = 0.0;
    float gy = -9.8;

#pragma omp parallel for num_threads(4)
    for (int pi = 0; pi < NUMPARTICLES; ++pi)
    {

        dpx = 0.0;
        dpy = 0.0;

        dviscx = 0.0;
        dviscy = 0.0;

        std::vector<Integer> sParticles;
        aHashGrid.find(sParticles, p[pi].pos, h);

        for (Integer n = 0; n < sParticles.size(); ++n)
        {

            Integer pj = sParticles[n];

            if (pi == pj)
                continue;

            CGeoVector<float> r = p[pj].pos - p[pi].pos;

            CGeoVector<float> gradW = aSpikyKernel.gradientW(r, h, 1E-8);

            dpx += p[pj].mass * (p[pi].p / (p[pi].rho * p[pi].rho) + p[pj].p / (p[pj].rho * p[pj].rho)) * gradW.x();
            dpy += p[pj].mass * (p[pi].p / (p[pi].rho * p[pi].rho) + p[pj].p / (p[pj].rho * p[pj].rho)) * gradW.y();

            float lapW = aSpikyKernel.laplacianW(r, h);

            dviscx += p[pj].mass * p[pi].visc / p[pi].rho * (p[pj].vel.x() - p[pi].vel.x()) / p[pj].rho * lapW;
            dviscy += p[pj].mass * p[pi].visc / p[pi].rho * (p[pj].vel.y() - p[pi].vel.y()) / p[pj].rho * lapW;

        }

        p[pi].vel.x() += dt * (gx - dpx + dviscx);
        p[pi].vel.y() += dt * (gy - dpy + dviscy);

    }

}

void advectParticles()
{

    buildHashGrid();

    calculateDensityAndPressure();

    //moveParticles();

}

void drawParticles()
{

    float r, g, b;

    float s;
    float smin = numeric_limits<float>::max();
    float smax = numeric_limits<float>::min();

    for (int pi = 0; pi < NUMPARTICLES; ++pi)
    {
        //s = sqrt(p[pi].vel.x() * p[pi].vel.x() + p[pi].vel.y() * p[pi].vel.y());
        //s = p[pi].p;
        s = p[pi].rho;
        //s = p[pi].factor;

        smin = std::min(smin, s);
        smax = std::max(smax, s);
    }

    std::cout << "smin = " << smin << std::endl;
    std::cout << "smax = " << smax << std::endl;

    glPointSize(6.0f);

    glBegin(GL_POINTS);

    for (int pi = 0; pi < NUMPARTICLES; ++pi)
    {
        //s = (sqrt(p[pi].vel.x() * p[pi].vel.x() + p[pi].vel.y() * p[pi].vel.y()) - smin) / (smax - smin);
        //s = (p[pi].p - smin) / (smax - smin);
        s = (p[pi].rho - smin) / (smax - smin);
        //s = (p[pi].factor - smin) / (smax - smin);

        getColor(s, r, g, b);
        glColor3f(r, g, b);
        glVertex2i((int)(p[pi].pos.x()*winwidth / (NX + 1)), (int)(p[pi].pos.y()*winheight / (NY + 1)));
    }

    glEnd();

}

void display()
{

    glClear(GL_COLOR_BUFFER_BIT);

    glPushMatrix();

    // setup 2d pixel plotting camera
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, (GLfloat)winwidth, 0.0f, (GLfloat)winheight, 0.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0, 0, winwidth, winheight);
    glDepthRange(0.0f, 1.0f);

    drawParticles();

    glPopMatrix();

    glutSwapBuffers();

}

void update()
{

    advectParticles();

    // redraw
    glutPostRedisplay();

    frames++;

    te = glutGet(GLUT_ELAPSED_TIME);

    // every second approximately
    if (te - t0 >= 1000)
    {

        char title[80];
        sprintf(title, "Smooth Particle Hydrodynamics Kernel test    %.1f fps", (1000.0*frames / (te - t0)));
        glutSetWindowTitle(title);

        frames = 0;
        t0 = te;

    }

}

void mouse(int button, int state, int x, int y)
{


}

void reshape(int w, int h)
{

    winwidth = w;
    winheight = h;

}

void init()
{

    // Add particles
    int pi = 0;

    for (int i = 0; i < NX; ++i)
    {

        for (int j = 0; j < NY; ++j)
        {

            CParticle pi;

            pi.pos.x() = 1 + i;
            pi.pos.y() = 1 + j;

            p.push_back(pi);

        }

    }

    for (int pi = 0; pi < NUMPARTICLES; ++pi)
    {

        p[pi].vel.x() = 0.0;
        p[pi].vel.y() = 0.0;

        p[pi].mass = m0;
        p[pi].visc = visc0;
        p[pi].rho = rho0;

    }

    walls.p1.x() = 0.0f;
    walls.p1.y() = 0.0f;

    walls.p2.x() = NX + 1;
    walls.p2.y() = NY + 1;

    scale[0][0] = 0.0 / 4.0;
    scale[0][1] = 0;
    scale[0][2] = 0;
    scale[0][3] = 1;

    scale[1][0] = 1.0 / 4.0;
    scale[1][1] = 0;
    scale[1][2] = 1;
    scale[1][3] = 1;

    scale[2][0] = 2.0 / 4.0;
    scale[2][1] = 0;
    scale[2][2] = 1;
    scale[2][3] = 0;

    scale[3][0] = 3.0 / 4.0;
    scale[3][1] = 1;
    scale[3][2] = 1;
    scale[3][3] = 0;

    scale[4][0] = 4.0 / 4.0;
    scale[4][1] = 1;
    scale[4][2] = 0;
    scale[4][3] = 0;

}

int main(int argc, char* argv[])
{

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(winwidth, winheight);

    glutCreateWindow("Smooth Particle Hydrodynamics Kernel test");

    glutDisplayFunc(display);
    glutIdleFunc(update);
    //glutMouseFunc(mouse);
    glutReshapeFunc(reshape);

    glClearColor(0, 0, 0, 0);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);

    init();

    glutMainLoop();

    return 0;

}
