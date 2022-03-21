#include "dialogchk.h"

dialogCHK::dialogCHK(QObject *parent) : QObject(parent)
{
    FRMdialogCHK = new QGroupBox(tr("Dialog"));
    BTNcreate = new QPushButton();
    connect(BTNcreate,SIGNAL(clicked()), this, SLOT(create()));
    BTNdelete = new QPushButton();
    connect(BTNdelete,SIGNAL(clicked()), this, SLOT(delet()));
}

void dialogCHK::create(){

}

void dialogCHK::delet(){

}
