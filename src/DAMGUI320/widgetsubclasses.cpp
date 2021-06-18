//  Copyright 2008-2019, Rafael Lopez, Alfredo Aguado, Octavio Roncero
//  This file is part of the package VIEWER
//
//  Author: Rafael Lopez
//  rafael.lopez@uam.es
//
//  Universidad Autonoma de Madrid, May 2018

#include "widgetsubclasses.h"

// Function for LineEdit
//
void LineEdit::keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_Enter || event->key() == Qt::Key_Return)
    {
        emit textchanged();
    }
    else {      // If other key is pressed, QLineEdit processes it as usual
        QLineEdit::keyPressEvent( event );
    }
}

// Function for DoubleSpinBox
//

void DoubleSpinBox::keyPressEvent(QKeyEvent *event)
{


    if(event->key() == Qt::Key_Enter || event->key() == Qt::Key_Return)
    {
        emit valueChanged(this->value());
    }
    else {      // If other key is pressed, QDoubleSpinBox processes it as usual
        QDoubleSpinBox::keyPressEvent( event );
    }
}

// Function for SpinBox
//

void SpinBox::keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_Enter || event->key() == Qt::Key_Return)
    {
        emit valueChanged(this->value());
    }
    else {      // If other key is pressed, QSpinBox processes it as usual
        QSpinBox::keyPressEvent( event );
    }
}

// Class editclosePushButton
//
editclosePushButton::editclosePushButton(){
    ind = 0;
    strings[0] = QString(tr("Edit"));
    strings[1] = QString(tr("Close"));
    current = strings[ind];
    setText(current);
}

editclosePushButton::~editclosePushButton(){

}

void editclosePushButton::toggletext()
{
    ind = (ind+1)%2;
    current = strings[ind];
    setText(current);
    if (ind == 0){
        emit isEdit(true);
        emit isClose(false);
    }
    else{
        emit isEdit(false);
        emit isClose(true);
    }
}

// Class showhidePushButton
//
showhidePushButton::showhidePushButton(){
    ind = 0;
    strings[0] = QString(tr("Show"));
    strings[1] = QString(tr("Hide"));
    current = strings[ind];
    setText(current);
}

showhidePushButton::~showhidePushButton(){

}

void showhidePushButton::toggletext()
{
    ind = (ind+1)%2;
    current = strings[ind];
    setText(current);
}

void showhidePushButton::inittext(int i)
{
    ind = i%2;
    current = strings[ind];
    setText(current);
}

// Class centerData
//

centerData::centerData(){
    molecule = -1;
    number = -2;
    type = -1;
    znuc = -1;
    symbol = "";
    x = 0;
    y = 0;
    xyz = QVector3D(0.,0.,0.);
}

centerData::~centerData()
{

}

// Function for class togglingGroupBox
//

void togglingGroupBox::toggleVisible(){
    setVisible(!isVisible());
    emit isvisible(isVisible());
}


// Class labelSlider
//
labelSlider::labelSlider(){
    LBLtext = new QLabel("");
    LBLvalue = new QLabel("");
    SLDslider = new QSlider(Qt::Horizontal);
    SLDslider->setRange(0,100);
    SLDslider->setSingleStep(1);
    SLDslider->setPageStep(10);
    SLDslider->setTickPosition(QSlider::TicksBelow);
    SLDslider->setValue(50.f);
    connect(SLDslider, SIGNAL(sliderReleased()), this, SIGNAL(sliderReleased()));
    connect(SLDslider, SIGNAL(valueChanged(int)), this, SIGNAL(valueChanged(int)));
    connect(SLDslider, SIGNAL(valueChanged(int)), this, SLOT(setlabelvalue(int)));

    QHBoxLayout *layout1 = new QHBoxLayout();
    layout1->addWidget(LBLtext,0,Qt::AlignLeft);
    layout1->addWidget(LBLvalue,0,Qt::AlignRight);

    QVBoxLayout *layout = new QVBoxLayout();
    layout->addLayout(layout1);
    layout->addWidget(SLDslider);
    layout->addStretch();

    this->setLayout(layout);
}

labelSlider::~labelSlider(){

}

void labelSlider::setlabeltext(QString str){
    LBLtext->setText(str);
}

void labelSlider::setlabelvalue(int i){
    LBLvalue->setNum(i);
}

void labelSlider::setRange(int i,int j){
    SLDslider->setRange(i,j);
}

void labelSlider::setSingleStep(int i){
    SLDslider->setSingleStep(i);
}

void labelSlider::setPageStep(int i){
    SLDslider->setPageStep(i);
}

void labelSlider::setValue(int i){
    SLDslider->setValue(i);
    LBLvalue->setText(QString("%1").arg(i));
}

// Class lineEditSlider
//
lineEditSlider::lineEditSlider(){
    LBLtext = new QLabel("");
    TXTedit = new LineEdit();
    SLDslider = new QSlider(Qt::Horizontal);

    TXTedit->setText("0.");
    vmax = 0.;
    TXTedit->setAlignment(Qt::AlignRight);
    connect(TXTedit, SIGNAL(textchanged()), this, SLOT(TXTedit_changed()));

    SLDslider->setRange(0,100);
    SLDslider->setSingleStep(1);
    SLDslider->setPageStep(10);
    SLDslider->setTickPosition(QSlider::TicksBelow);
    SLDslider->setValue(50.f);
    connect(SLDslider, SIGNAL(sliderReleased()), this, SLOT(SLDslider_released()));
    connect(SLDslider, SIGNAL(valueChanged(int)), this, SLOT(SLDslider_changed(int)));

    QHBoxLayout *layout1 = new QHBoxLayout();
    layout1->addWidget(LBLtext,0,Qt::AlignLeft);
    layout1->addWidget(TXTedit,Qt::AlignRight);

    QVBoxLayout *layout = new QVBoxLayout();
    layout->addLayout(layout1);
    layout->addWidget(SLDslider);
    layout->addStretch();

    this->setLayout(layout);
}

lineEditSlider::~lineEditSlider(){

}
float lineEditSlider::getscalevalueFloat(int ix,int i0,int i1,float f0,float f1)
{
    if (ix > i1) ix = i1;
    if (ix < i0) ix = i0;
    return (f1-f0) * (ix-i0) / (i1-i0) + f0;
}

int lineEditSlider::getscalevalueInt(float fx,int i0,int i1,float f0,float f1)
{
    if (fx > f1) fx = f1;
    if (fx < f0) fx = f0;
    return (int)(((fx-f0)*(i1-i0)/(f1-f0))+i0);
}

void lineEditSlider::setlabeltext(QString str){
    LBLtext->setText(str);
}

void lineEditSlider::setValidator(QValidator *val){
    TXTedit->setValidator(val);
}

void lineEditSlider::setRange(int i,int j){
    SLDslider->setRange(i,j);
}

void lineEditSlider::setSingleStep(int i){
    SLDslider->setSingleStep(i);
}

void lineEditSlider::setPageStep(int i){
    SLDslider->setPageStep(i);
}

void lineEditSlider::setTXTtext(QString str){
    TXTedit->setText(str);
}

void lineEditSlider::setValue(float v){
    vmax = v;
    int ix = getscalevalueInt(TXTedit->text().toFloat(), SLDslider->minimum(), SLDslider->maximum(), 0., vmax);
    SLDslider->setValue(ix);
}

void lineEditSlider::SLDslider_changed(int ix){
    topcolor = getscalevalueFloat(ix,SLDslider->minimum(), SLDslider->maximum(), 0., vmax);
    TXTedit->setText(QString::number(topcolor));
}

void lineEditSlider::SLDslider_released(){
    int ix = SLDslider->sliderPosition();
    topcolor = getscalevalueFloat(ix,SLDslider->minimum(), SLDslider->maximum(), 0., vmax);
    TXTedit->setText(QString::number(topcolor));
    emit valueChanged(topcolor);
    emit sliderReleased();
}

void lineEditSlider::TXTedit_changed(){
    topcolor = TXTedit->text().toFloat();
    int ix = getscalevalueInt(TXTedit->text().toFloat(), SLDslider->minimum(), SLDslider->maximum(), 0., vmax);
//    topcolor = qMin(TXTedit->text().toFloat(),vmax);
    SLDslider->setTracking(false);
    SLDslider->setSliderPosition(ix);
    SLDslider->setTracking(true);
    emit valueChanged(topcolor);
    emit sliderReleased();
}

// Class combineSignalsOr
//

combineSignalsOr::combineSignalsOr(){
    signal1 = false;
    signal2 = false;
}

combineSignalsOr::~combineSignalsOr(){

}

void combineSignalsOr::setSignal1(bool a){
    signal1 = a;
    emit signalsCombined(signal1 || signal2);
}

void combineSignalsOr::setSignal2(bool a){
    signal2 = a;
    emit signalsCombined(signal1 || signal2);
}


/*******************************************************************************************************/
/********************************  Class myScrollArea  implementation  *******************************/
/*******************************************************************************************************/

myScrollArea::myScrollArea(QWidget *parent) : QScrollArea(parent)
{
//    QRect rec = QApplication::desktop()->availableGeometry();
    QScreen *screen = QGuiApplication::primaryScreen();
    QRect  rec = screen->geometry();
    maxheight = rec.height();
    maxwidth = rec.width();
    height = 0;
    width = 0;
}

myScrollArea::~myScrollArea(){

}

void myScrollArea::reject(){
    emit closed();
}

void myScrollArea::closeEvent(QCloseEvent *event){
    event->ignore();
    emit closed();
    event->accept();
}

void myScrollArea::updatesize(QSize a){
    height = qMin(maxheight,a.height()+15);
    width = qMin(maxwidth,a.width()+15);
    resize(QSize(width,height));
}

