// base
#include "sphconfig.h"

// sph
#include "sph.h"
#include "vec3.h"


// i want to avoid having a 'config class' the sph part will use.
// sph.cpp should be mostly standalone (apart from particle.*).
// that's why i went for this quite ugly solution


SphConfig::SphConfig(QWidget* parent)
 : QTreeWidget(parent)
{
//   QFont font;
//   font.setFamily(QString::fromUtf8("Courier"));
//   setFont(font);

   QStringList labels;
   labels << "name";
   labels << "value";
   setHeaderLabels(labels);
   setColumnCount(2);

   QList<QTreeWidgetItem*> items;

   mGravityX = new QTreeWidgetItem((QTreeWidget*)0, QStringList("gravity x"));
   mGravityY = new QTreeWidgetItem((QTreeWidget*)0, QStringList("gravity y"));
   mGravityZ = new QTreeWidgetItem((QTreeWidget*)0, QStringList("gravity z"));

   mStiffness = new QTreeWidgetItem((QTreeWidget*)0, QStringList("stiffness"));
   mViscosity = new QTreeWidgetItem((QTreeWidget*)0, QStringList("viscosity"));
   mDamping = new QTreeWidgetItem((QTreeWidget*)0, QStringList("damping"));
   mTimeStep = new QTreeWidgetItem((QTreeWidget*)0, QStringList("timestep"));
   mCflLimit = new QTreeWidgetItem((QTreeWidget*)0, QStringList("cfl"));

   items << mGravityX;
   items << mGravityY;
   items << mGravityZ;
   items << mStiffness;
   items << mViscosity;
   items << mDamping;
   items << mTimeStep;
   items << mCflLimit;

   for (int i = 0; i < items.size(); i++)
      items[i]->setFlags(mGravityX->flags() | Qt::ItemIsEditable);

   insertTopLevelItems(0, items);
}



void SphConfig::readValuesFromSimulation()
{
   vec3 gravity = mSph->getGravity();
   float stiffness = mSph->getStiffness();
   float viscosity = mSph->getViscosityScalar();
   float damping = mSph->getDamping();
   float timeStep = mSph->getTimeStep();
   float cfl = mSph->getCflLimit();

   mGravityX->setText(1, QString::number(gravity.x, 'f'));
   mGravityY->setText(1, QString::number(gravity.y, 'f'));
   mGravityZ->setText(1, QString::number(gravity.z, 'f'));
   mStiffness->setText(1, QString::number(stiffness, 'f'));
   mViscosity->setText(1, QString::number(viscosity, 'f'));
   mDamping->setText(1, QString::number(damping, 'f'));
   mTimeStep->setText(1, QString::number(timeStep, 'f'));
   mCflLimit->setText(1, QString::number(cfl, 'f'));
}


void SphConfig::writeValuesToSimulation()
{
   vec3 gravity;
   float gravityX = mGravityX->text(1).toFloat();
   float gravityY = mGravityY->text(1).toFloat();
   float gravityZ = mGravityZ->text(1).toFloat();
   gravity.set(gravityX, gravityY, gravityZ);
   float stiffness = mStiffness->text(1).toFloat();
   float viscosity = mViscosity->text(1).toFloat();
   float damping = mDamping->text(1).toFloat();
   float timestep = mTimeStep->text(1).toFloat();
   float cfl = mCflLimit->text(1).toFloat();

   mSph->setGravity(gravity);
   mSph->setStiffness(stiffness);
   mSph->setViscosityScalar(viscosity);
   mSph->setDamping(damping);
   mSph->setTimeStep(timestep);
   mSph->setCflLimit(cfl);
}


void SphConfig::setSph(SPH *sph)
{
   mSph = sph;
}



