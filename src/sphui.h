#ifndef SPHUI_H
#define SPHUI_H

#include <QFrame>

namespace Ui {
class SPHUI;
}

class SPHUI : public QFrame
{
   Q_OBJECT

public:
   explicit SPHUI(QWidget *parent = 0);
   ~SPHUI();

private:
   Ui::SPHUI *ui;
};

#endif // SPHUI_H
