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
#include "SphCubicSpline.hpp"
#include "SphGaussian.hpp"
#include "SphQuintic.hpp"
#include "SphSpiky.hpp"
#include "SphConvex.hpp"

using namespace std;

using namespace ENigMA::geometry;
using namespace ENigMA::sph;

const int NX = 50;         // x
const int NY = 50;         // y

float dt = 2E-3f;      // time step
float h = 0.05f;
float m0 = 1.0f;
float rho0 = 1000.0f;
float visc0 = 1.0f;

struct Particle
{
    CGeoCoordinate<float> pos;
    CGeoVector<float> vel;

    float mass;
    float p;
    float rho;
    float visc;
};

std::vector<Particle> p;

int winwidth = 480;
int winheight = 480;

float scale[6][4];

int frames = 0;
int t0 = 0, te;

CGeoHashGrid<float> aHashGrid;
CSphCubicSpline<float> aCubicSplineKernel(2);
CSphSpiky<float> aSpikyKernel(2);
CSphGaussian<float> aGaussianKernel(2);
CSphQuintic<float> aQuinticKernel(2);
CSphConvex<float> aConvexKernel(2);

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
    for (int pi = 0; pi < p.size(); ++pi)
        aHashGrid.addGeometricObject(pi, p[pi].pos);

    aHashGrid.build();

}

void calculateDensityAndPressure()
{

    int gamma = 7;
    float cs = 343.0;

#pragma omp parallel for num_threads(4)
    for (int pi = 0; pi < p.size(); ++pi)
    {

        p[pi].rho = 0.0;

        std::vector<Integer> sParticles;
        aHashGrid.find(sParticles, p[pi].pos, h);

        for (Integer n = 0; n < sParticles.size(); ++n)
        {

            Integer pj = sParticles[n];

            CGeoVector<float> r = p[pj].pos - p[pi].pos;

            p[pi].rho += p[pj].mass * aGaussianKernel.W(r, h);

        }

        if (p[pi].rho - rho0 > 0)
            p[pi].p = cs * cs / gamma * (float)(pow(p[pi].rho / rho0, gamma) - 1.0);
        else
            p[pi].p = 0.0;

    }

}

void moveParticles()
{

    float dpx, dpy;
    float dviscx, dviscy;

    float gx = 0.0;
    float gy = 0.0;

#pragma omp parallel for num_threads(4)
    for (int pi = 0; pi < p.size(); ++pi)
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

            CGeoVector<float> gradW = aGaussianKernel.gradientW(r, h, 1E-8);

            dpx += p[pj].mass * (p[pi].p / (p[pi].rho * p[pi].rho) + p[pj].p / (p[pj].rho * p[pj].rho)) * gradW.x();
            dpy += p[pj].mass * (p[pi].p / (p[pi].rho * p[pi].rho) + p[pj].p / (p[pj].rho * p[pj].rho)) * gradW.y();

            float lapW = aGaussianKernel.laplacianW(r, h);

            dviscx += p[pj].mass * p[pi].visc / p[pi].rho * (p[pj].vel.x() - p[pi].vel.x()) / p[pj].rho * lapW;
            dviscy += p[pj].mass * p[pi].visc / p[pi].rho * (p[pj].vel.y() - p[pi].vel.y()) / p[pj].rho * lapW;
        }

        p[pi].vel.x() += dt * (gx - dpx + dviscx);
        p[pi].vel.y() += dt * (gy - dpy + dviscy);

        p[pi].pos.x() += dt * p[pi].vel.x();
        p[pi].pos.y() += dt * p[pi].vel.y();

        // Bound positions
        if (p[pi].pos.x() <= 0.0f)
        {
            p[pi].pos.x() = 0.01;
            p[pi].vel.x() = 0.0;
            p[pi].vel.y() = 0.0;
        }

        if (p[pi].pos.x() >= 1.0f)
        {
            p[pi].pos.x() = 0.99;
            p[pi].vel.x() = 0.0;
            p[pi].vel.y() = 0.0;
        }

        if (p[pi].pos.y() <= 0.0f)
        {
            p[pi].pos.y() = 0.01;
            p[pi].vel.x() = 0.0;
            p[pi].vel.y() = 0.0;
        }

        if (p[pi].pos.y() >= 1.0f)
        {
            p[pi].pos.y() = 0.99;
            p[pi].vel.x() = 0.0;
            p[pi].vel.y() = 0.0;
        }

        if (p[pi].pos.y() >= 0.98f)
        {
            p[pi].vel.x() = 1.0f;
            p[pi].vel.y() = 0.0f;
        }
    }

}

void advectParticles()
{

    buildHashGrid();

    calculateDensityAndPressure();

    moveParticles();

}

void drawParticles()
{

    float r, g, b;

    float s;
    float smin = numeric_limits<float>::max();
    float smax = numeric_limits<float>::min();

    for (int pi = 0; pi < p.size(); ++pi)
    {
        s = p[pi].rho;
        smin = std::min(smin, s);
        smax = std::max(smax, s);
    }

    glPointSize(6.0f);

    glBegin(GL_POINTS);

    for (int pi = 0; pi < p.size(); ++pi)
    {
        s = (p[pi].rho - smin) / (smax - smin);
        getColor(s, r, g, b);
        glColor3f(r, g, b);
        glVertex2f(p[pi].pos.x(), p[pi].pos.y());
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
    glOrtho(0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f);
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
        sprintf(title, "Smooth Particle Hydrodynamics demo    %.1f fps", (1000.0*frames / (te - t0)));
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

    // Add fluid particles
    for (int i = 0; i < NX; ++i)
    {

        for (int j = 0; j < NY; ++j)
        {

            Particle pi;

            pi.pos.x() = +i / (float) (NX - 1);
            pi.pos.y() = +j / (float) (NY - 1);

            pi.vel.x() = 0.0;
            pi.vel.y() = 0.0;

            pi.mass = m0;
            pi.visc = visc0;
            pi.rho = rho0;

            p.push_back(pi);

        }

    }

    // Bottom wall
    for (int i = 0; i < 2 * NX; ++i)
    {
        Particle pi;

        pi.pos.x() = +i / (float)(2 * NX - 1);
        pi.pos.y() = -1 / (float)(2 * NY - 1);

        pi.vel.x() = 0.0;
        pi.vel.y() = 0.0;

        pi.mass = m0;
        pi.visc = visc0;
        pi.rho = rho0;

        p.push_back(pi);
    }

    // Side walls
    for (int j = 0; j < 2 * NY; ++j)
    {
        {
            Particle pi;

            pi.pos.x() = -1 / (float)(2 * NX - 1);
            pi.pos.y() = +j / (float)(2 * NY - 1);

            pi.vel.x() = 0.0;
            pi.vel.y() = 0.0;

            pi.mass = m0;
            pi.visc = visc0;
            pi.rho = rho0;

            p.push_back(pi);
        }

        {
            Particle pi;

            pi.pos.x() = 1.0 + 1 / (float)(2 * NX - 1);
            pi.pos.y() = +j / (float)(2 * NY - 1);

            pi.vel.x() = 0.0;
            pi.vel.y() = 0.0;

            pi.mass = m0;
            pi.visc = visc0;
            pi.rho = rho0;

            p.push_back(pi);
        }
    }

    // Corners
    {
        Particle pi;

        pi.pos.x() = -1 / (float)(2 * NX - 1);
        pi.pos.y() = -1 / (float)(2 * NY - 1);

        pi.vel.x() = 0.0;
        pi.vel.y() = 0.0;

        pi.mass = 10 * m0;
        pi.visc = visc0;
        pi.rho = rho0;

        p.push_back(pi);
    }

    {
        Particle pi;

        pi.pos.x() = 1 + 1 / (float)(2 * NX - 1);
        pi.pos.y() = -1 / (float)(2 * NY - 1);

        pi.vel.x() = 0.0;
        pi.vel.y() = 0.0;

        pi.mass = 10 * m0;
        pi.visc = visc0;
        pi.rho = rho0;

        p.push_back(pi);
    }

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

    glutCreateWindow("Smooth Particle Hydrodynamics demo");

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
