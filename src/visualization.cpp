// base
#include "visualization.h"

// glu
#include <GL/glu.h>

// Qt
#include <QTimer>

// sph
#include "particle.h"
#include "sph.h"


Visualization::Visualization(QWidget *parent) :
   QGLWidget(parent),
   mWidth(0),
   mHeight(0),
   mDrawParticlesChecked(false),
   mDrawVoxelsChecked(false),
   mFrame(0)
{
   QTimer* timer = new QTimer(this);

   connect(
      timer,
      SIGNAL(timeout()),
      this,
      SLOT(updateGL())
   );

   timer->start(16);
   mElapsed.start();
}


Visualization::~Visualization()
{
}


void Visualization::closeEvent(QCloseEvent * /*e*/)
{
   emit closed();
}


void Visualization::initializeGL()
{
   glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

   glEnable(GL_BLEND);

   // line smoothing
   glEnable(GL_LINE_SMOOTH);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

   glPointSize(0.5f);

}


void Visualization::drawBox(float maxX, float maxY, float maxZ)
{
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   /*
             ____
           /|a  d|\
         e/ |    | \h
          | |    | |
          | |b__c| |
          | /    \ |
          |/______\|
          f        g
   */

   vec3 a(0.0f, maxY, maxZ);
   vec3 b(0.0f, 0.0f, maxZ);
   vec3 c(maxX, 0.0f, maxZ);
   vec3 d(maxX, maxY, maxZ);

   vec3 e(0.0f, maxY, 0.0f);
   vec3 f(0.0f, 0.0f, 0.0f);
   vec3 g(maxX, 0.0f, 0.0f);
   vec3 h(maxX, maxY, 0.0f);

   // draw boundaries
   glBegin(GL_LINES);

   glColor4f(1.0f, 1.0f, 1.0f, 0.2f);

   // back
   glVertex3f(a.x, a.y, a.z);
   glVertex3f(b.x, b.y, b.z);

   glVertex3f(b.x, b.y, b.z);
   glVertex3f(c.x, c.y, c.z);

   glVertex3f(c.x, c.y, c.z);
   glVertex3f(d.x, d.y, d.z);

   glVertex3f(d.x, d.y, d.z);
   glVertex3f(a.x, a.y, a.z);

   // connections
   glVertex3f(a.x, a.y, a.z);
   glVertex3f(e.x, e.y, e.z);

   glVertex3f(d.x, d.y, d.z);
   glVertex3f(h.x, h.y, h.z);

   glVertex3f(b.x, b.y, b.z);
   glVertex3f(f.x, f.y, f.z);

   glVertex3f(c.x, c.y, c.z);
   glVertex3f(g.x, g.y, g.z);

   // front
   glVertex3f(e.x, e.y, e.z);
   glVertex3f(f.x, f.y, f.z);

   glVertex3f(f.x, f.y, f.z);
   glVertex3f(g.x, g.y, g.z);

   glVertex3f(g.x, g.y, g.z);
   glVertex3f(h.x, h.y, h.z);

   glVertex3f(h.x, h.y, h.z);
   glVertex3f(e.x, e.y, e.z);

   glEnd();
}


void Visualization::drawParticles()
{
   glBlendFunc(GL_ONE, GL_ONE);

   QTime elapsed;
   elapsed.start();

   int count = mSph->getParticleCount();
   Particle* particles = mSph->getParticles();

   glBegin(GL_POINTS);

   glColor4f(0.0f, 0.2f, 0.8f, 1.0f);

   for (int i = 0; i < count; i++)
   {
      Particle* particle = &particles[i];

      const vec3& pos = particle->mPosition;
      glVertex3f(pos.x, pos.y, pos.z);
   }

   glEnd();

   emit updateElapsed(elapsed.elapsed());
}


void Visualization::drawVoxels()
{
   glBlendFunc(GL_ONE, GL_ONE);

   int x = 0;
   int y = 0;
   int z = 0;
   int index = 0;
   int count = 0;
   float cellSize = mSph->getCellSize();
   mSph->getGridCellCounts(x, y, z);

   QList<Particle*>* grid = mSph->getGrid();

   glBegin(GL_TRIANGLES);

   for (int zi = 0; zi < z; zi++)
   {
      for (int yi = 0; yi < y; yi++)
      {
         for (int xi = 0; xi < x; xi++)
         {
            index =
                 (xi        )
               + (yi * x    )
               + (zi * x * y);

            count = grid[index].count();

            if (count > 0)
            {
               glColor4f(count * 0.02f, 0.0f, 0.0f, 1.0f);

               drawVoxel(
                  xi      * cellSize,
                 (xi + 1) * cellSize,
                  yi      * cellSize,
                 (yi + 1) * cellSize,
                  zi      * cellSize,
                 (zi + 1) * cellSize
               );
            }
         }
      }
   }

   glEnd();
}


void Visualization::drawVoxel(
   float xMin,
   float xMax,
   float yMin,
   float yMax,
   float zMin,
   float zMax
)
{
   /*
             ____
           /|a  d|\
         e/_|____|_\h
          | |    | |
          | |b__c| |
          | /    \ |
          |/______\|
          f        g
   */

   vec3 a(xMin, yMax, zMax);
   vec3 b(xMin, yMin, zMax);
   vec3 c(xMax, yMin, zMax);
   vec3 d(xMax, yMax, zMax);

   vec3 e(xMin, yMax, zMin);
   vec3 f(xMin, yMin, zMin);
   vec3 g(xMax, yMin, zMin);
   vec3 h(xMax, yMax, zMin);

   // back
   glVertex3f(a.x, a.y, a.z);
   glVertex3f(c.x, c.y, c.z);
   glVertex3f(b.x, b.y, b.z);

   glVertex3f(a.x, a.y, a.z);
   glVertex3f(d.x, d.y, d.z);
   glVertex3f(c.x, c.y, c.z);

   // front
   glVertex3f(h.x, h.y, h.z);
   glVertex3f(e.x, e.y, e.z);
   glVertex3f(f.x, f.y, f.z);

   glVertex3f(h.x, h.y, h.z);
   glVertex3f(f.x, f.y, f.z);
   glVertex3f(g.x, g.y, g.z);


   // left
   glVertex3f(e.x, e.y, e.z);
   glVertex3f(a.x, a.y, a.z);
   glVertex3f(b.x, b.y, b.z);

   glVertex3f(e.x, e.y, e.z);
   glVertex3f(b.x, b.y, b.z);
   glVertex3f(f.x, f.y, f.z);

   // right
   glVertex3f(h.x, h.y, h.z);
   glVertex3f(g.x, g.y, g.z);
   glVertex3f(d.x, d.y, d.z);

   glVertex3f(d.x, d.y, d.z);
   glVertex3f(g.x, g.y, g.z);
   glVertex3f(c.x, c.y, c.z);

   // bottom
   glVertex3f(g.x, g.y, g.z);
   glVertex3f(f.x, f.y, f.z);
   glVertex3f(b.x, b.y, b.z);

   glVertex3f(b.x, b.y, b.z);
   glVertex3f(g.x, g.y, g.z);
   glVertex3f(c.x, c.y, c.z);

   // top: aeh, ahd
   glVertex3f(a.x, a.y, a.z);
   glVertex3f(e.x, e.y, e.z);
   glVertex3f(h.x, h.y, h.z);

   glVertex3f(a.x, a.y, a.z);
   glVertex3f(h.x, h.y, h.z);
   glVertex3f(d.x, d.y, d.z);
}


bool Visualization::isDrawVoxelsChecked() const
{
   return mDrawVoxelsChecked;
}


void Visualization::setDrawVoxelsChecked(bool drawVoxelsChecked)
{
   mDrawVoxelsChecked = drawVoxelsChecked;
}


bool Visualization::isDrawParticlesChecked() const
{
   return mDrawParticlesChecked;
}


void Visualization::setDrawParticlesChecked(bool drawParticlesChecked)
{
   mDrawParticlesChecked = drawParticlesChecked;
}


void Visualization::paintGL()
{
   float maxX;
   float maxY;
   float maxZ;
   mSph->getParticleBounds(maxX, maxY, maxZ);
   float invScaleX = 1.0f / maxX;
   float invScaleY = 1.0f / maxY;
   float invScaleZ = 1.0f / maxZ;

   glClear(GL_COLOR_BUFFER_BIT);

   glViewport(0, 0, width(), height());

   // setup projection
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity ();

   glFrustum(
      -0.5f, 1.5f,
      -0.5f, 1.5f,
      0.0f, 10.0f
   );

   gluLookAt(
      0.0f, 0.0f, -0.5f,
      0.5f, 0.5f, 1.0f,
      0.0f, 1.0f, 0.0f
   );

   // setup modelview
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glRotatef(180.0f, 0.0f, 1.0f, 0.0f);
   glTranslatef(-0.5f, -0.5f, 0.0f);
   glScalef(invScaleX, invScaleY, invScaleZ);

   // draw bounding box
   drawBox(maxX, maxY, maxZ);

   // draw particles
   if (isDrawParticlesChecked())
      drawParticles();

   if (isDrawVoxelsChecked())
      drawVoxels();
}


SPH *Visualization::getSph() const
{
   return mSph;
}


void Visualization::setSph(SPH *sph)
{
   mSph = sph;
}


