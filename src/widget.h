#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>


class Visualization;
class SphConfig;


namespace Ui {
class Widget;
}

class Widget : public QWidget
{
   Q_OBJECT


public:

   explicit Widget(QWidget *parent = 0);
   ~Widget();

   SphConfig* Config() const;

   Visualization* getVisualization();

   void updateElapsed();


signals:

   void startClicked();
   void shutDownClicked();


public slots:

   void updateElapsedSph(int, int, int, int, int, int);
   void updateElapsedVisualization(int elapsed);


private:

   Ui::Widget *ui;

   int mElapsedVoxelize;
   int mElapsedFindNeighbors;
   int mElapsedComputeDensity;
   int mElapsedComputePressure;
   int mElapsedComputeAcceleration;
   int mElapsedIntegrate;

   int mElapsedVisualization;
};

#endif // WIDGET_H
