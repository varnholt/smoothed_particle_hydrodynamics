#ifndef VISUALIZATION_H
#define VISUALIZATION_H

// Qt
#include <QGLWidget>
#include <QTime>

// forward declarations
class SPH;


class Visualization : public QGLWidget
{
   Q_OBJECT

public:

   explicit Visualization(QWidget *parent = 0);
   ~Visualization();

   SPH *getSph() const;
   void setSph(SPH *sph);

   bool isDrawParticlesChecked() const;
   bool isDrawVoxelsChecked() const;


signals:

   void closed();
   void updateElapsed(int);


public slots:

   void setDrawVoxelsChecked(bool drawVoxelsChecked);
   void setDrawParticlesChecked(bool drawParticlesChecked);


protected:

   void closeEvent(QCloseEvent* e);

   //! initialize gl
   virtual void initializeGL();

   //! paint
   virtual void paintGL();

   void drawBox(
      float maxX,
      float maxY,
      float maxZ
   );

   void drawParticles();

   void drawVoxels();

   void drawVoxel(
      float xMin,
      float xMax,
      float yMin,
      float yMax,
      float zMin,
      float zMax
   );

private:

   QTime mElapsed;
   SPH* mSph;

   float mWidth;
   float mHeight;
   float mAspect;

   bool mDrawParticlesChecked;
   bool mDrawVoxelsChecked;

   int mFrame;
};

#endif // VISUALIZATION_H
