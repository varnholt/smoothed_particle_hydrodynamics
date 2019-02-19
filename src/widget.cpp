// base
#include "widget.h"
#include "ui_widget.h"

// Qt
#include <QApplication>

// sph
#include "visualization.h"


Widget::Widget(QWidget *parent) :
   QWidget(parent),
   ui(new Ui::Widget),
   mElapsedVisualization(0),
   mElapsedVoxelize(0),
   mElapsedFindNeighbors(0),
   mElapsedComputeDensity(0),
   mElapsedComputePressure(0),
   mElapsedComputeAcceleration(0),
   mElapsedIntegrate(0)
{
   ui->setupUi(this);

   connect(
      ui->mRun,
      SIGNAL(clicked()),
      this,
      SIGNAL(startClicked())
   );

   connect(
      ui->mQuit,
      SIGNAL(clicked()),
      this,
      SIGNAL(shutDownClicked())
   );

   connect(
      ui->mVisualization,
      SIGNAL(updateElapsed(int)),
      this,
      SLOT(updateElapsedVisualization(int))
   );

   connect(
      ui->mApply,
      SIGNAL(clicked()),
      ui->mAttributes,
      SLOT(writeValuesToSimulation())
   );

   connect(
      ui->mDrawParticles,
      SIGNAL(clicked(bool)),
      ui->mVisualization,
      SLOT(setDrawParticlesChecked(bool))
   );

   connect(
      ui->mDrawVoxels,
      SIGNAL(clicked(bool)),
      ui->mVisualization,
      SLOT(setDrawVoxelsChecked(bool))
   );

   ui->mVisualization->setDrawParticlesChecked(ui->mDrawParticles->isChecked());
   ui->mVisualization->setDrawVoxelsChecked(ui->mDrawVoxels->isChecked());
}


Widget::~Widget()
{
   delete ui;
}



SphConfig *Widget::Config() const
{
   return ui->mAttributes;
}


Visualization *Widget::getVisualization()
{
   return ui->mVisualization;
}


void Widget::updateElapsed()
{
   ui->mText->setPlainText(
      QString(
         "voxelize:\t\t%1\n"
         "find neighbors:\t\t%2\n"
         "compute density:\t%3\n"
         "compute pressure:\t%4\n"
         "compure acceleration:\t%5\n"
         "integrate:\t\t%6\n"
         "draw:\t\t%7\n"
      )
      .arg(mElapsedVoxelize)
      .arg(mElapsedFindNeighbors)
      .arg(mElapsedComputeDensity)
      .arg(mElapsedComputePressure)
      .arg(mElapsedComputeAcceleration)
      .arg(mElapsedIntegrate)
      .arg(mElapsedVisualization)
   );
}


void Widget::updateElapsedSph(
   int timeVoxelize,
   int timeFindNeighbors,
   int timeComputeDensity,
   int timeComputePressure,
   int timeComputeAcceleration,
   int integrate
)
{
   mElapsedVoxelize = timeVoxelize;
   mElapsedFindNeighbors = timeFindNeighbors;
   mElapsedComputeDensity = timeComputeDensity;
   mElapsedComputePressure = timeComputePressure;
   mElapsedComputeAcceleration = timeComputeAcceleration;
   mElapsedIntegrate = integrate;

   updateElapsed();
}


void Widget::updateElapsedVisualization(int elapsed)
{
   mElapsedVisualization = elapsed;
   updateElapsed();
}



