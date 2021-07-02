//  Copyright 2008-2021, Rafael Lopez, Alfredo Aguado, Octavio Roncero
//  This file is part of the package VIEWER
//
//  Author: Rafael Lopez
//  rafael.lopez@uam.es
//
//  Universidad Autonoma de Madrid, May 2018

#ifndef WIDGETSUBCLASSES_H
#define WIDGETSUBCLASSES_H

#include <QApplication>
#include <QObject>
#include <QVector3D>
#include <QWidget>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QLineEdit>
#include <QMatrix4x4>
#include <QPushButton>
#include <QScreen>
#include <QSlider>
#include <QScrollArea>

#include "ColorButton.h"

#include <QDebug>

//  The following subclassing of QlineEdit is necessary to get a correct behavior when return key is pressed,
//  because using the QlineEdit returnPressed signal causes undesired toggle of Show/Hide button
//
class LineEdit : public QLineEdit
{
    Q_OBJECT
public:
    void keyPressEvent(QKeyEvent *event);
signals:
    void textchanged();
};


//  The following subclassing of QDoubleSpinBox is necessary to get a correct behavior when return key is pressed,
//  because using the QDoubleSpinBox returnPressed signal causes undesired toggle of Show/Hide button
//
class DoubleSpinBox : public QDoubleSpinBox
{
    Q_OBJECT
public:
    void keyPressEvent(QKeyEvent *event);
signals:
    void valueChanged(double);
};

//  The following subclassing of editclosePushButton implements a slot for toggling the text (Edit/Close) of a QPushButton
//
class editclosePushButton : public ColorButton
{
    Q_OBJECT
public:
    editclosePushButton();
    ~editclosePushButton();

signals:
    void isEdit(bool);
    void isClose(bool);
public slots:
    void toggletext();
private:
    int ind;
    QString strings[2];
    QString current;
};

//  The following subclassing of editclosePushButton implements a slot for changing the text (Show/Hide) of a QPushButton
//
class showhidePushButton : public ColorButton
{
    Q_OBJECT
public:
    showhidePushButton();
    ~showhidePushButton();

signals:

public slots:
    void inittext(int);
    void toggletext();
private:
    int ind;
    QString strings[2];
    QString current;
};

//  The following subclassing of labelSlider implements a slider with a title displaying numerical value
//
class labelSlider : public QWidget
{
    Q_OBJECT
public:
    labelSlider();
    ~labelSlider();
    void setlabeltext(QString);
    void setRange(int,int);
    void setSingleStep(int);
    void setPageStep(int);
    void setValue(int);

    QSlider *SLDslider;

signals:
    void sliderReleased();
    void valueChanged(int);

public slots:
    void setlabelvalue(int);

private:
    QLabel *LBLtext;
    QLabel *LBLvalue;

};

//  The following subclassing of lineEditSlider implements a slider with a line editor box
class lineEditSlider : public QWidget
{
    Q_OBJECT
public:
    lineEditSlider();
    ~lineEditSlider();
    float getscalevalueFloat(int ix,int i0,int i1,float f0,float f1);
    int  getscalevalueInt(float fx,int i0,int i1,float f0,float f1);
    void setValidator(QValidator *);
    void setlabeltext(QString);
    void setRange(int,int);
    void setSingleStep(int);
    void setPageStep(int);
    void setTXTtext(QString);
    void setValue(float);

    QSlider *SLDslider;

signals:
    void sliderReleased();
    void valueChanged(float);

public slots:
    void SLDslider_changed(int);
    void SLDslider_released();
    void TXTedit_changed();

private:
    float topcolor;
    float vmax;
    LineEdit *TXTedit;
    QLabel *LBLtext;
};

//  The following subclassing of QSpinBox is necessary to get a correct behavior when return key is pressed,
//  because using the QSpinBox returnPressed signal causes undesired toggle of Show/Hide button
//
class SpinBox : public QSpinBox
{
    Q_OBJECT
public:
    void keyPressEvent(QKeyEvent *event);
signals:
    void valueChanged(int);
};

// The following class packs data of a center
//
class centerData
{
public:
    centerData();
    ~centerData();

    int type;       // Center type: 0: atom; 1: cp; 2: mesp local extremum
    int molecule;
    int number;
    int znuc;
    int x;
    int y;
    QString symbol;
    QVector3D xyz;
    QMatrix4x4 *mtrsf;
};

// The following class packs data of a center
//
class togglingGroupBox : public QGroupBox
{
    Q_OBJECT
public:

signals:
    void isvisible(bool);

public slots:
    void toggleVisible();
};


//  The following subclassing of combineSignalsOr implements a combination of two signals with OR operator and emits a signal with the result
class combineSignalsOr : public QWidget
{
    Q_OBJECT
public:
    combineSignalsOr();
    ~combineSignalsOr();

signals:
    void signalsCombined(bool);

public slots:
    void setSignal1(bool);
    void setSignal2(bool);

private:
    bool signal1;
    bool signal2;
};


class myScrollArea : public QScrollArea
{
    Q_OBJECT
public:
    explicit myScrollArea(QWidget *parent = 0);
    ~myScrollArea();
    void updatesize(QSize);
signals:
    void closed();

protected:
    void closeEvent(QCloseEvent *event);
    virtual void reject();

private:
    int height;
    int maxheight;
    int maxwidth;
    int width;
};


#endif // WIDGETSUBCLASSES_H
