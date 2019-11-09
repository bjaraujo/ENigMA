
#include <iostream>

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <GL/glut.h>

#include <LbmLidDrivenSolver.hpp>

using namespace ENigMA::lbm;

static CLbmLidDrivenSolver<double, 2, 9> aLbmSolver(128, 128, 1);

using namespace std;

#define NUMPARTICLES    32768

const int NX = 128;         // x
const int NY = 128;         // y
const double Re = 1000;     // top cover speed
const double U = 0.1;       // lid velocity
const double dt = 1.0;      // time step
const double s = 50.0;      // only for particle advection

double p[NUMPARTICLES][2];

GLuint texnum;

int winwidth = 480;
int winheight = 480;

double dx, dy, Lx, Ly, P0, tau_f, niu, error;

double scale[6][4];

int frames = 0;
int t0 = 0, te;

bool addVelocity = false;
bool addBoundary = false;

int oldx, oldy;

void getColor(double v, double& r, double& g, double& b)
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

        for (i=1; i <= 4; ++i)
            if (scale[i-1][0] <= v && v < scale[i][0])
                break;

        // linear interpolation
        double w = (v - scale[i-1][0]) / (scale[i][0] - scale[i-1][0]);

        r = (1-w)*scale[i-1][1] + w*scale[i][1];
        g = (1-w)*scale[i-1][2] + w*scale[i][2];
        b = (1-w)*scale[i-1][3] + w*scale[i][3];

    }
}

void advectParticles()
{

    for (int pi = 0; pi < NUMPARTICLES; ++pi)
    {

        int i = (int)p[pi][0] % NX;
        int j = (int)p[pi][1] % NY;

        if (i <= 0 || i >= NX || j <= 0 || j >= NY)
        {
            p[pi][0] = i = 1 + (int)(rand() % NX - 1);
            p[pi][1] = j = 1 + (int)(rand() % NY - 1);
        }

        p[pi][0] += aLbmSolver.getVelocity(i, j, 0) * dt * s; 
        p[pi][1] += aLbmSolver.getVelocity(i, j, 1) * dt * s;

    }

}

void drawParticles()
{

    glBegin(GL_POINTS);

    for (int pi=0; pi < NUMPARTICLES; ++pi)
    {
        glColor3f(1.0, 1.0, 0.0);
        glVertex2i((int)(p[pi][0]*winwidth/(NX+1)), (int)(p[pi][1]*winheight/(NY+1)));
    }
    glEnd();

}

void drawVelocities()
{

    unsigned char bitmap[(NX+1)*(NY+1)*4];    // rgba unsigned bytes

    double m, r, g, b;

    for (int j = 0; j <= NY; ++j)
    {
        for (int i = 0; i <= NX; ++i)
        {

            if (aLbmSolver.getBoundary(i, j) == 1)
            {

                r = g = b = 0;

            }
            else
            {

                double uu = aLbmSolver.getVelocity(i, j, 0);
                double vv = aLbmSolver.getVelocity(i, j, 1);

                m = sqrt(uu * uu + vv * vv);

                getColor(m / U, r, g, b);

            }

            bitmap[(i+j*(NX + 1))*4 + 0] = static_cast<unsigned char>(r*255);
            bitmap[(i+j*(NX + 1))*4 + 1] = static_cast<unsigned char>(g*255);
            bitmap[(i+j*(NX + 1))*4 + 2] = static_cast<unsigned char>(b*255);
            bitmap[(i+j*(NX + 1))*4 + 3] = 255;

        }
    }

    glBindTexture(GL_TEXTURE_2D, texnum);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NX+1, NY+1, 0, GL_RGBA, GL_UNSIGNED_BYTE, bitmap);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  
    glEnable(GL_TEXTURE_2D);

    glBegin(GL_QUADS);

        glColor3f(1.0, 1.0, 1.0);

        glTexCoord2f(0.0, 0.0);
        glVertex2i(0, 0);

        glTexCoord2f(1.0, 0.0);
        glVertex2i(winwidth, 0);

        glTexCoord2f(1.0, 1.0);
        glVertex2i(winwidth, winheight);

        glTexCoord2f(0.0, 1.0);
        glVertex2i(0, winheight);

    glEnd();

    glDisable(GL_TEXTURE_2D);

}

void display()
{

    glClear(GL_COLOR_BUFFER_BIT);

    glPushMatrix();

    // setup 2d pixel plotting camera
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, (GLfloat) winwidth, 0.0f, (GLfloat) winheight, 0.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0, 0, winwidth, winheight);
    glDepthRange(0.0f, 1.0f);

    drawVelocities();

    drawParticles();
    
    glPopMatrix();

    glutSwapBuffers();

}

void update()
{

    double dx = 1.0;
    double dy = 1.0;
    
    double Lx = dx * (double) NX;
    double Ly = dy * (double) NY;
    
    double niu = U * Lx / Re;
    double tau_f = 3.0 * niu + 0.5;

    aLbmSolver.evolve(tau_f);

    advectParticles();

    // redraw
    glutPostRedisplay();

    frames++;

    te = glutGet(GLUT_ELAPSED_TIME);
    
    // every second approximately
    if (te-t0 >= 1000)
    {

        char title[80];
        sprintf(title, "Lattice Boltzmann lid-driven cavity demo    %.1f fps", (1000.0*frames/(te-t0)));
        glutSetWindowTitle(title);

        frames = 0;
        t0 = te;

    }

}

void mouse(int button, int state, int x, int y)
{

    addVelocity = false;
    addBoundary = false;

    x *= (int)((double)(NX+1)/winwidth);
    y *= (int)((double)(NY+1)/winheight);

    if (state == GLUT_DOWN)
    {

        if (button == GLUT_LEFT_BUTTON)
        {
            if (x > 0 && x < NX && y > 0 && y < NY)
            {
                addVelocity = true;
                oldx = x;
                oldy = y;
            }
        }

        if (button == GLUT_RIGHT_BUTTON)
        {
            if (x >= 0 && x <= NX && y >= 0 && y <= NY)
            {
                addBoundary = true;
                oldx = x;
                oldy = y;
            }
        }
    }

}

void motion(int x, int y)
{

    double m, um[2] = {0, 0};

    x *= (int)((double)(NX+1)/winwidth);
    y *= (int)((double)(NY+1)/winheight);

    if (addBoundary && (x >= 0 && x <= NX && y >= 0 && y <= NY))
    {

        int i = (int)x;
        int j = NY - (int)y;

        aLbmSolver.setBoundary(i, j, 1);

    }

    if (addVelocity && (x > 0 && x < NX && y > 0 && y < NY))
    {

        um[0] = (x-oldx);
        um[1] = (oldy-y);

        m = sqrt(um[0]*um[0]+um[1]*um[1]);

        int i = (int)x;
        int j = NY - (int)y;

        um[0] /= (1+10*m);
        um[1] /= (1+10*m);

        aLbmSolver.setVelocity(i, j, 0, um[0]);
        aLbmSolver.setVelocity(i, j, 1, um[1]);

    }

}

void reshape(int w, int h)
{

    winwidth = w;
    winheight = h;

}

void output_velocities()
{

    cout << "------- Nornalized center velocity -------" << endl;

    for (int y = 0; y < NY; y++)
        cout << (double) y / NY << " " << aLbmSolver.getVelocity(NX/2, y, 0) / U << endl;

}

void reset()
{

    for (int i = 0; i <= NX; ++i)
        for (int j = 0; j <= NY; ++j)
            aLbmSolver.setBoundary(i, j, 0);

}

void keyboard(unsigned char key, int x, int y)
{

    if (key == 27)
        exit(0);
    
    if (key == 'o')
        output_velocities();

    if (key == 'r')
        reset();

}

void init()
{

    aLbmSolver.init();

    for (int i = 0; i <= NX; ++i)
    {

        for (int j = 0; j <= NY; ++j)
        {

            aLbmSolver.setVelocity(i, j, 0, 0);
            aLbmSolver.setVelocity(i, j, 1, 0);
            aLbmSolver.setDensity(i, j, 1);

            aLbmSolver.setBoundary(i, j, 0);

            if (i == 0 || i == NX || j == 0)
                aLbmSolver.setBoundary(i, j, 1);

        }

        aLbmSolver.setVelocity(i, NY, 0, U);

    }

    for (int pi=0; pi < NUMPARTICLES; ++pi)
    {
        // dont spawn on the boundary
        p[pi][0] = 1 + (int)(rand() % NX - 1);
        p[pi][1] = 1 + (int)(rand() % NY - 1);
    }

    scale[0][0] = 0.0/4.0;
    scale[0][1] = 0;
    scale[0][2] = 0;
    scale[0][3] = 1;

    scale[1][0] = 1.0/4.0;
    scale[1][1] = 0;
    scale[1][2] = 1;
    scale[1][3] = 1;

    scale[2][0] = 2.0/4.0;
    scale[2][1] = 0;
    scale[2][2] = 1;
    scale[2][3] = 0;

    scale[3][0] = 3.0/4.0;
    scale[3][1] = 1;
    scale[3][2] = 1;
    scale[3][3] = 0;

    scale[4][0] = 4.0/4.0;
    scale[4][1] = 1;
    scale[4][2] = 0;
    scale[4][3] = 0;

}

int main(int argc, char* argv[])
{

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(winwidth, winheight);
    
    glutCreateWindow("Lattice Boltzmann demo");

    glutDisplayFunc(display);
    glutIdleFunc(update);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    glClearColor(0, 0, 0, 0);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    
    glGenTextures(1, &texnum);

    init();

    glutMainLoop();

    return 0;
    
}



