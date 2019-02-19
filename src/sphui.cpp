#include "sphui.h"
#include "ui_sphui.h"

SPHUI::SPHUI(QWidget *parent) :
   QFrame(parent),
   ui(new Ui::SPHUI)
{
   ui->setupUi(this);
}

SPHUI::~SPHUI()
{
   delete ui;
}
