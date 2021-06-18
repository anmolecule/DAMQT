#include <QtCore/qmath.h>
#include <QDialog>
#include <QCheckBox>
#include <QColorDialog>
#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QFontDialog>
#include <QGroupBox>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QVBoxLayout>

#include "isosurface.h"
#include "CIsoSurface.h"
#include <math.h>
#include "GlobalInfo.h"

isosurface::isosurface(QWidget *parent) : QWidget(parent)
{
//#if defined(Q_WS_WIN) || defined(Q_OS_WIN)
//    iswindows = true;    // To be used by fortran programs
//    mpi = false;
//#else //Q_WS_X11, Q_WS_MAC
//    iswindows = false;
//#endif
//    qDebug() << "iswindows = " << iswindows;
//    qDebug() << "mpi = " << mpi;
//    qDebug() << "maxnumprocessors = " << maxnumprocessors;
    allindices.clear();
    allvertices.clear();
    connections.clear();
    griddimensions.clear();
    gridindices.clear();
    gridindicesoffset.clear();
    gridvertices.clear();
    BTNexec = nullpointer;
    BTNstop = nullpointer;
    BTNsurfcolor = nullpointer;
    CHKmpi = nullpointer;
    CHKnormalgrad = nullpointer;
    CHKshowgrid = nullpointer;
    CHKtranslucence = nullpointer;
    FRMhighquality = nullpointer;
    FRMisosurface = nullpointer;
    FRMsurfcolor = nullpointer;
    FRMsurftype = nullpointer;
    initialposition = QPoint(200,300);
    LBLalpha = nullpointer;
    LBLcontourvalue = nullpointer;
    LBLfilename = nullpointer;
    LBLmpi = nullpointer;
    LBLopacity = nullpointer;
    LBLscale = nullpointer;
    LBLsensitive = nullpointer;
    LBLstatus = nullpointer;
    myDoubleValidator = nullpointer;
    myProcess = nullpointer;
    RBTscalelin = nullpointer;
    RBTscalelog = nullpointer;
    RBTsolidsurf = nullpointer;
    RBTwiresurf = nullpointer;
    SLDcontourvalue = nullpointer;
    SLDopacity = nullpointer;
    SPBmpi = nullpointer;
    SPBopacity = nullpointer;
    SPBsensitive = nullpointer;
    TXTcontourvalue = nullpointer;
    TXTisosurffilename = nullpointer;

    contourvalue = 0.;
    logdlt = 3;
    mincontourvalue = -1.0;
    maxcontourvalue = 1.0;
    nprocessors = 1;
    opacity = 1.0;

    surfcolor = QColor(245,0,0);
    editoropen = false;
    isdensity = true;
    solidsurf = true;
    settranslucence(false);
    setvisible(true);
    showgridbound = false;
    normalgrad = false;
    compatderiv = false;

    basename = "";
}

isosurface::~isosurface(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
}

//  ------------------------------------------------------------------------------------------------------------------
//
//          Surface editor
//
//  ------------------------------------------------------------------------------------------------------------------

togglingGroupBox * isosurface::editisosurface(){
    if (!FRMisosurface)
        FRMisosurface = new togglingGroupBox();
    //        Contour value
    if (!myDoubleValidator)
        myDoubleValidator = new QDoubleValidator();
    myDoubleValidator->setLocale(QLocale::English);
    myDoubleValidator->setRange(mincontourvalue,maxcontourvalue,6);

    if (!LBLcontourvalue)
        LBLcontourvalue = new QLabel(tr("Value")+":");

    if (!TXTcontourvalue)
        TXTcontourvalue = new LineEdit();
    TXTcontourvalue->setValidator(myDoubleValidator);
    TXTcontourvalue->setAlignment(Qt::AlignRight);
    TXTcontourvalue->setText(QString("%1").arg(contourvalue));
    connections << connect(TXTcontourvalue, SIGNAL(textchanged()), this, SLOT(TXTcontourvalue_changed()));

    if (!LBLscale)
        LBLscale = new QLabel(tr("Scale")+":");

    if (!RBTscalelog)
        RBTscalelog = new QRadioButton(tr("Logarithmic"));
    RBTscalelog->setChecked(true);

    if (!RBTscalelin)
        RBTscalelin = new QRadioButton(tr("Linear"));
    connections << connect(RBTscalelog, SIGNAL(toggled(bool)), this, SLOT(RBTscale_changed()));

    if (!SLDcontourvalue)
        SLDcontourvalue = new QSlider(Qt::Horizontal);
    SLDcontourvalue->setRange(0,1000);
    SLDcontourvalue->setSliderPosition(getscalevalueInt(contourvalue,SLDcontourvalue->minimum(),SLDcontourvalue->maximum(),
                                    mincontourvalue,maxcontourvalue, RBTscalelog->isChecked()));
    connections << connect(SLDcontourvalue, SIGNAL(valueChanged(int)), this, SLOT(SLDcontourvalue_changed(int)));
    connections << connect(SLDcontourvalue, SIGNAL(sliderReleased()), this, SLOT(SLDcontourvalue_released()));

    if (!LBLsensitive)
        LBLsensitive = new QLabel(tr("Sensitiveness")+":");
    LBLsensitive->setVisible(true);

    if (!SPBsensitive)
        SPBsensitive = new QSpinBox();
    SPBsensitive->setMinimum(1);
    SPBsensitive->setMaximum(3);
    SPBsensitive->setMaximumWidth(50);
    SPBsensitive->setValue(3);
    SPBsensitive->setVisible(true);
    logdlt = SPBsensitive->value();
    connections << connect(SPBsensitive,SIGNAL(valueChanged(int)),this,SLOT(SPBsensitive_changed(int)));

//        Surface color and type
    if (!FRMsurfcolor)
        FRMsurfcolor = new QGroupBox(tr("Color and type"));

    if (!BTNsurfcolor)
        BTNsurfcolor = new ColorButton(QIcon(":/images/colores48.png"),tr("Color"));
    BTNsurfcolor->setColor(&surfcolor);
    connections << connect(BTNsurfcolor, SIGNAL(clicked()), this, SLOT(BTNsurfcolor_clicked()));

    if (!FRMsurftype)
        FRMsurftype = new QGroupBox(tr("Surface type"));

    if (!RBTsolidsurf)
        RBTsolidsurf = new QRadioButton(tr("Solid surface"));
    RBTsolidsurf->setChecked(solidsurf);

    if (!RBTwiresurf)
        RBTwiresurf = new QRadioButton(tr("Wire frame"));
    RBTwiresurf->setChecked(!solidsurf);
    connections << connect(RBTsolidsurf, SIGNAL(toggled (bool)), this, SLOT(RBTsurftype_changed()));

//        Opacity
    if (!LBLopacity)
        LBLopacity = new QLabel();
    LBLopacity->setText(tr("Opacity:"));

    if (!LBLalpha)
        LBLalpha = new QLabel();
    LBLalpha->setText(QString("%1").arg((int)100.f * opacity));

    if (!SLDopacity)
        SLDopacity = new QSlider(Qt::Horizontal);
    SLDopacity->setRange(0,100);
    SLDopacity->setSingleStep(1);
    SLDopacity->setPageStep(10);
    SLDopacity->setTickPosition(QSlider::TicksBelow);
    SLDopacity->setValue(100.f * opacity);
    connections << connect(SLDopacity, SIGNAL(valueChanged(int)), this, SLOT(SLDopacity_changed(int)));
    connections << connect(SLDopacity, SIGNAL(valueChanged(int)), LBLalpha, SLOT(setNum(int)));
    connections << connect(SLDopacity, SIGNAL(sliderReleased()), this, SLOT(SLDopacity_released()));

//        Translucency correction
    if (!CHKtranslucence)
        CHKtranslucence = new QCheckBox(tr("Translucence correction"));
    CHKtranslucence->setChecked(false);
    connections << connect(CHKtranslucence, SIGNAL(stateChanged(int)), this, SLOT(CHKtranslucence_changed()));

//        Normals computed from interpolated gradient
    if (!CHKnormalgrad)
        CHKnormalgrad = new QCheckBox(tr("Compute normals from interpolated gradient"));
    CHKnormalgrad->setChecked(normalgrad);
    CHKnormalgrad->setVisible(compatderiv);
    connections << connect(CHKnormalgrad, SIGNAL(stateChanged(int)), this, SLOT(CHKnormalgrad_changed()));

    //        Show grid
    if (!CHKshowgrid)
        CHKshowgrid = new QCheckBox(tr("Show grid boundaries"));
    CHKshowgrid->setChecked(false);
    showgridbound = false;
    connections << connect(CHKshowgrid, SIGNAL(stateChanged(int)), this, SLOT(CHKshowgrid_changed()));


    if (!FRMhighquality)
        FRMhighquality = new QGroupBox(tr("High quality surface"));
    FRMhighquality->setToolTip(tr("Computes the isosurface with analytical gradient and stores it in a file to be loaded with the Add surface option"));
    FRMhighquality->setStyleSheet("QGroupBox::title{padding 2 2000}");

    if (!LBLfilename)
        LBLfilename = new QLabel(tr("File: "),FRMhighquality);

    if (!TXTisosurffilename)
        TXTisosurffilename = new QLineEdit(FRMhighquality);
    TXTisosurffilename->setText(basename);

    if (!CHKmpi)
        CHKmpi = new QCheckBox(tr("MPI"),FRMhighquality);

    if (!LBLmpi)
        LBLmpi = new QLabel(tr("Number of processors"),FRMhighquality);

    if (!SPBmpi)
        SPBmpi = new QSpinBox(FRMhighquality);
    SPBmpi->setRange(1, maxnumprocessors);
    SPBmpi->setValue(nprocessors);
    SPBmpi->setMaximumWidth(60);
    SPBmpi->setToolTip(tr("Number of processors"));
    if (mpi){
        CHKmpi->setChecked(true);
        SPBmpi->setEnabled(true);
    }
    else{
        CHKmpi->setHidden(true);
        CHKmpi->setChecked(false);
        SPBmpi->setHidden(true);
        SPBmpi->setEnabled(false);
    }
    connections << connect(CHKmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKmpi_changed(int)));
    connections << connect(SPBmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBmpi_changed(int)));

    if (!BTNexec)
        BTNexec = new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),FRMhighquality);
    connections << connect(BTNexec, SIGNAL(clicked()), this, SLOT(BTNexec_clicked()));

    if (!BTNstop)
        BTNstop = new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),FRMhighquality);
    connections << connect(BTNstop, SIGNAL(clicked()), this, SLOT(processStop()));

    if (!LBLstatus)
        LBLstatus = new QLabel(tr(""),FRMhighquality);

//        Layouts

    QHBoxLayout *layout1 = new QHBoxLayout();
    layout1->addStretch();
    layout1->addWidget(LBLcontourvalue);
    layout1->addWidget(TXTcontourvalue);
    layout1->addStretch();

    QHBoxLayout *layout2 = new QHBoxLayout();
    layout2->addStretch();
    layout2->addWidget(LBLscale);
    layout2->addWidget(RBTscalelog);
    layout2->addWidget(RBTscalelin);
    layout2->addStretch();

    QHBoxLayout *layout3 = new QHBoxLayout();
    layout3->addStretch();
    layout3->addWidget(LBLsensitive);
    layout3->addWidget(SPBsensitive);
    layout3->addStretch();

    QVBoxLayout *layout4 = new QVBoxLayout();
    layout4->addStretch();
    layout4->addLayout(layout1);
    layout4->addWidget(SLDcontourvalue);
    layout4->addLayout(layout2);
    layout4->addLayout(layout3);
    layout4->addStretch();

    QHBoxLayout *layout5 = new QHBoxLayout();
    layout5->addStretch();
    layout5->addWidget(BTNsurfcolor);
    layout5->addStretch();

    QHBoxLayout *layout6 = new QHBoxLayout(FRMsurftype);
    layout6->addStretch();
    layout6->addWidget(RBTsolidsurf);
    layout6->addWidget(RBTwiresurf);
    layout6->addStretch();

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addWidget(LBLopacity,0,Qt::AlignLeft);
    layout7->addWidget(LBLalpha,0,Qt::AlignRight);

    QVBoxLayout *layout8 = new QVBoxLayout();
    layout8->addLayout(layout7);
    layout8->addWidget(SLDopacity);
    layout8->addStretch();

    QHBoxLayout *layout9 = new QHBoxLayout();
    layout9->addStretch();
    layout9->addWidget(CHKtranslucence);
    layout9->addStretch();

    QVBoxLayout *layout10 = new QVBoxLayout(FRMsurfcolor);
    layout10->addLayout(layout5);
    layout10->addWidget(FRMsurftype);
    layout10->addLayout(layout8);
    layout10->addLayout(layout9);

    QHBoxLayout *Layout11 = new QHBoxLayout();
    Layout11->addWidget(LBLfilename);
    Layout11->addWidget(TXTisosurffilename);

    QHBoxLayout *Layout12 = new QHBoxLayout();
    Layout12->addWidget(CHKmpi);
    Layout12->addWidget(LBLmpi);
    Layout12->addWidget(SPBmpi);

    QHBoxLayout *Layout13 = new QHBoxLayout();
    Layout13->addStretch();
    Layout13->addWidget(BTNstop,0,Qt::AlignRight);
    Layout13->addWidget(BTNexec,0,Qt::AlignRight);

    QHBoxLayout *Layout14 = new QHBoxLayout();
    Layout14->addWidget(LBLstatus,0,Qt::AlignLeft);
    Layout13->addStretch();

    QVBoxLayout *Layout15 = new QVBoxLayout(FRMhighquality);
    Layout15->addLayout(Layout11);
    Layout15->addLayout(Layout12);
    Layout15->addLayout(Layout13);
    Layout15->addLayout(Layout14);

    QVBoxLayout *layout = new QVBoxLayout(FRMisosurface);
    layout->addStretch();
    layout->addLayout(layout4);
    layout->addWidget(FRMsurfcolor);
    layout->addWidget(CHKnormalgrad);
    layout->addWidget(CHKshowgrid);
    layout->addWidget(FRMhighquality);
    layout->addStretch();

    editoropen = true;
    return FRMisosurface;
}

void isosurface::closeeditor(){
    if (!editoropen){
        return;
    }
    delete LBLstatus;
    LBLstatus = nullpointer;
    delete LBLfilename;
    LBLfilename = nullpointer;
    delete TXTisosurffilename;
    TXTisosurffilename = nullpointer;
    delete SPBmpi;
    SPBmpi = nullpointer;
    delete LBLmpi;
    LBLmpi = nullpointer;
    delete CHKmpi;
    CHKmpi = nullpointer;
    delete BTNexec;
    BTNexec = nullpointer;
    delete BTNstop;
    BTNstop = nullpointer;
    delete BTNsurfcolor;
    BTNsurfcolor = nullpointer;
    delete CHKnormalgrad;
    CHKnormalgrad = nullpointer;
    delete CHKshowgrid;
    CHKshowgrid = nullpointer;
    delete CHKtranslucence;
    CHKtranslucence = nullpointer;
    delete LBLalpha;
    LBLalpha = nullpointer;
    delete LBLcontourvalue;
    LBLcontourvalue = nullpointer;
    delete LBLopacity;
    LBLopacity = nullpointer;
    delete LBLscale;
    LBLscale = nullpointer;
    delete LBLsensitive;
    LBLsensitive = nullpointer;
    delete myDoubleValidator;
    myDoubleValidator = nullpointer;
    delete RBTscalelin;
    RBTscalelin = nullpointer;
    delete RBTscalelog;
    RBTscalelog = nullpointer;
    delete RBTsolidsurf;
    RBTsolidsurf = nullpointer;
    delete RBTwiresurf;
    RBTwiresurf = nullpointer;
    delete SLDcontourvalue;
    SLDcontourvalue = nullpointer;
    delete SLDopacity;
    SLDopacity = nullpointer;
    delete SPBsensitive;
    SPBsensitive = nullpointer;
    delete TXTcontourvalue;
    TXTcontourvalue = nullpointer;
    delete FRMhighquality;
    FRMhighquality = nullpointer;
    delete FRMsurftype;
    FRMsurftype = nullpointer;
    delete FRMsurfcolor;
    FRMsurfcolor = nullpointer;
    delete FRMisosurface;
    FRMisosurface = nullpointer;

    editoropen = false;
}

//  ------------------------------------------------------------------------------------------------------------------
//
//          Functions and slots for surface editor
//
//  ------------------------------------------------------------------------------------------------------------------

// Function getscalevalueInt
//    Returns slider position ix corresponding to function value fx.
//      i0, i1: range of slider scale
//      f0, f1: function values at i0 and i1
//      logs: if true: logatihmic scale, if false: linear scale
int isosurface::getscalevalueInt(float fx,int i0,int i1,float f0,float f1, bool logs)
{
    if (fx > f1) fx = f1;
    if (fx < f0) fx = f0;
    int ix;

    if (logs==true){
        float dlt = pow(10.,-logdlt);
        ix = round(((log(fx-f0+dlt)-log(dlt))/(log(f1-f0+dlt)-log(dlt)))*(i1-i0) + i0);
    }else{
        ix = (int)(((fx-f0)*(i1-i0)/(f1-f0))+i0);
    }
    return ix;
}

// Function getscalevalueFloat
//    Returns function value fx corresponding to slider position ix.
//      i0, i1: range of slider scale
//      f0, f1: function values at i0 and i1
//      logs: if true: logatihmic scale, if false: linear scale
float isosurface::getscalevalueFloat(int ix,int i0,int i1,float f0,float f1, bool logs)
{
    if (ix > i1) ix = i1;
    if (ix < i0) ix = i0;
    float fx;

    if (logs==true){
        float dlt = pow(10.,-logdlt);
        fx = exp((log(f1-f0+dlt)-log(dlt))*(ix-i0)/(i1-i0)+log(dlt)) + f0 - dlt;
    }else{
        fx = (f1-f0) * (ix-i0) / (i1-i0) + f0;
    }
    return fx;
}

void isosurface::toggleshowsurf(){
    setvisible(!isvisible());
    emit generatesurface();
    emit updatedisplay();
}

void isosurface::BTNexec_clicked(){
    QString path = QFileInfo(fullname).path();
    QString suffix = QFileInfo(name).suffix();

    if (!QFile::exists(QString(path+"/"+ProjectName+"_2016.damqt"))){
        QMessageBox msgBox;
        msgBox.setText(tr("High quality isosurface"));
        msgBox.setInformativeText(tr("File %1 does not exist. Cannot do operation. ").arg(path+"/"+ProjectName+"_2016.damqt"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    if (!isdensity){
        if (!QFile::exists(QString(path+"/"+ProjectName+"_2016.dmqtv"))){
            QMessageBox msgBox;
            msgBox.setText(tr("High quality isosurface"));
            msgBox.setInformativeText(tr("File %1 does not exist. Cannot do operation. ").arg(path+"/"+ProjectName+"_2016.dmqtv"));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return;
        }
    }
    QString stdinput;
    if (isdensity){
        stdinput = path+"/"+basename+"-DAMISODEN320.inp";
    }
    else{
        stdinput = path+"/"+basename+"-DAMISOPOT320.inp";
    }

    QFile inputfile(stdinput);
    if (!inputfile.isOpen()){
        inputfile.open(QFile::Text | QFile::WriteOnly);
    }
    QTextStream infile(&inputfile);
    infile << "&OPTIONS\n";
    infile << " lmaxrep = 15\n";
    infile << " contourval = " << QString("%1").arg(contourvalue) << " \n";
    if (isdensity){
        infile << " gridname = \"" << basename+"-d.plt" << "\" \n";
        infile << " filename = \"" << TXTisosurffilename->text()+"-d" << "\" \n";
    }
    else{
        infile << " gridname = \"" << basename+"-v.plt" << "\" \n";
        infile << " filename = \"" << TXTisosurffilename->text()+"-v" << "\" \n";
    }
    if (iswindows){
        infile << " iswindows = T\n";
    }
    else{
        infile << " iswindows = F\n";
    }
    infile << " lbinary = T\n";
    infile << "&END\n";
    infile << "\""+path+"/"+ProjectName+"\"";
    inputfile.close();
    QString stdoutput;
    QString strprocess;
    if (isdensity){
        if (mpi && CHKmpi->isChecked()){
            stdoutput = path+"/"+basename+"-DAMISODEN320_mpi.out";
            processname = "DAMISODEN320_mpi.exe";
            QString execName = get_execName(processname, QString("DAM320_mpi"));
            strprocess = QString("%1 %2 %3 %4 %5").arg(mpicommand).arg("-np")
                        .arg(SPBmpi->value()).arg(mpiflags).arg(execName);
//            qDebug() << "strprocess = " << strprocess;
        }
        else{
            stdoutput = path+"/"+basename+"-DAMISODEN320.out";
            processname = "DAMISODEN320.exe";
            QString execName = get_execName(processname, QString("DAM320"));
            strprocess = execName;
//            qDebug() << "strprocess = " << strprocess;
        }
    }
    else{
        if (mpi && CHKmpi->isChecked()){
            stdoutput = path+"/"+basename+"-DAMISOPOT320_mpi.out";
            processname = "DAMISOPOT320_mpi.exe";
            QString execName = get_execName(processname, QString("DAM320_mpi"));
            strprocess = QString("%1 %2 %3 %4 %5").arg(mpicommand).arg("-np")
                        .arg(SPBmpi->value()).arg(mpiflags).arg(execName);
        }
        else{
            stdoutput = path+"/"+basename+"-DAMISOPOT320.out";
            processname = "DAMISOPOT320.exe";
            QString execName = get_execName(processname, QString("DAM320"));
            strprocess = execName;
        }
    }
    if (myProcess){
        delete myProcess;
        myProcess = nullpointer;
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    myProcess->start(strprocess);
    return;
}

void isosurface::BTNsurfcolor_clicked(){
    QColor currcolor = surfcolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        setsurfcolor(newcolor);
        BTNsurfcolor->setColor(&newcolor);
    }
    int r = newcolor.red();
    int g = newcolor.green();
    int b = newcolor.blue();
    QString lblstring;
    if ((r > 210 || b > 0)&& g > 210 ){
        lblstring = QString("QPushButton { color : rgb(0, 0, 0); background-color : rgb(%1,%2,%3); }").arg(r).arg(g).arg(b);
    }
    else{
        lblstring = QString("QPushButton { color : rgb(255, 255, 255); background-color : rgb(%1,%2,%3); }").arg(r).arg(g).arg(b);
    }
    emit generatesurface();
    emit updatelabelcolor(lblstring);
    emit updatedisplay();
}

void isosurface::CHKnormalgrad_changed(){
    normalgrad = CHKnormalgrad->isChecked();
    emit generatesurface();
    emit updatedisplay();
}

void isosurface::CHKshowgrid_changed(){
    if (CHKshowgrid->isChecked())
        showgridbound = true;
    else
        showgridbound = false;
    emit updatedisplay();
}

void isosurface::CHKmpi_changed(int state)
{
    if (CHKmpi->isChecked()){
        LBLmpi->setEnabled(true);
        SPBmpi->setEnabled(true);
    }
    else{
        LBLmpi->setEnabled(false);
        SPBmpi->setEnabled(false);
    }
}

void isosurface::CHKtranslucence_changed(){
    settranslucence(!gettranslucence());
    emit updatedisplay();
}

void isosurface::emitupdateRightMenu(){
    emit updateRightMenu();
}

QString isosurface::get_execName(QString processname, QString subdir){
#if defined(Q_WS_WIN) || defined(Q_OS_WIN)
    QString execName = processname;
    if (!QFileInfo::exists(execName)){
        execName = QCoreApplication::applicationDirPath()+"/"+processname;
    }
#else
    QString execName = QCoreApplication::applicationDirPath()+"/"+processname;
    if (!QFileInfo::exists(execName)){
        execName = processname;
    }
#endif
    if (!QFileInfo::exists(execName))
        execName = QCoreApplication::applicationDirPath()+"/../"+subdir+"/"+processname;
    if (!QFileInfo::exists(execName)){
        QString message1, message2, message3;
        QString direc = QString(QCoreApplication::applicationDirPath());
        direc.truncate(direc.lastIndexOf(QChar('/')));
        message1 = QString(tr("Executable file %1 does not exist\n\n").arg(processname));
        message2 = QString(tr("Check that the program is installed in any of the following directories: \n\n %1 \n %2 \n\n")
                        .arg(QCoreApplication::applicationDirPath()+"/")
                        .arg(direc+"/"+subdir+"/"));
        message3 = QString(tr("or in any other directory available in your $PATH"));
        int messagelen = qMax(qMax(message1.length(),message2.length()),message3.length());
        QMessageBox msg;
        msg.setText(message1+message2+message3);
        msg.setIcon(QMessageBox::Critical);
        QSpacerItem* horizontalSpacer = new QSpacerItem(messagelen * 4, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
        QGridLayout* layout = (QGridLayout*)msg.layout();
        layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
        msg.exec();
        return QString("");
    }
    return execName;
}


//    Starts an external process
void isosurface::processStart()
{
//    statusBar()->showMessage(tr("Computing..."));
    BTNexec->setEnabled(false);
    LBLstatus->setText(tr("Computing..."));
}


//    Stops the external process currently running
void isosurface::processStop()
{
    #if defined(Q_WS_WIN) || defined(Q_OS_WIN) || QT_VERSION < 0x050000
        myProcess->kill();
    #else
        if (myProcess->arguments().count() > 2){
            QString procname = myProcess->arguments().at(2);
            if (procname.contains("_mpi")){      // In case of mpi runs, kill all associated processes
                procname = procname.split('.').first();
                QProcess getprocesses;
                QString pgrep;
                pgrep = QString("pgrep -f %1").arg(procname);
                getprocesses.start(pgrep);
                getprocesses.waitForFinished();
                QByteArray procnumbers = getprocesses.readAllStandardOutput();
                QString procstring(procnumbers);
                QStringList procslist = procstring.split("\n");
                QMessageBox msgBox;
                msgBox.setInformativeText(QString(tr("Do you want to kill all processes named %1?").arg(procname)+"\n"
                       + tr("The following processes will be killed: "))+procstring.replace("\n"," "));
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
                msgBox.setDefaultButton(QMessageBox::Cancel);
                msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
                msgBox.setButtonText(QMessageBox::No, tr("No"));
                msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
                msgBox.setIcon(QMessageBox::Warning);
                int ret = msgBox.exec();
                if (ret == QMessageBox::Yes){
                    for (int i = 1; i < procslist.length()-1 ; i++){  // Kill the processes one by one (the last element in procslist is empty)
                        QProcess killprocess;
                        QString strprocess;
                        strprocess = QString("kill -9 %1").arg(procslist.at(i));
                        killprocess.start(strprocess);
                        killprocess.waitForFinished();
                    }
                }
                myProcess->kill();
            }
        }
        else{
            myProcess->kill();
        }
    #endif
    BTNexec->setEnabled(true);
    LBLstatus->setText(tr("Process stopped"));
}


//    Returns the error code of the external process currently running
void isosurface::processError(QProcess::ProcessError error)
{
    if(error==QProcess::FailedToStart){
        QString message;
        message = QString("Error %1").arg(error)
                + QString(tr("Process failed to start program %1\n").arg(processname))
                + QString(tr("Check that the program is installed in the directory: \n %1 \n")
                              .arg(QCoreApplication::applicationDirPath()))
                + QString(tr("or in any other directory in your $PATH"));
        QMessageBox::critical(this,QString("Error %1").arg(error),message);
    }
    else{
        QMessageBox::critical(this,QString("Error %1").arg(error),tr("Error when running program %1").arg(processname));
    }
    BTNexec->setEnabled(true);
    LBLstatus->setText(tr("Process failed"));
}

//    Slot to be run when a process ends
void isosurface::processOutput(int exitCode, QProcess::ExitStatus exitStatus){
    if(!(exitStatus==QProcess::NormalExit)){
        QMessageBox::information(this,QString(tr("High quality isosurface")),
                    QString(tr("Process %1 crashed, exit code = %2")).arg(processname).arg(exitCode));
    }
    else{
        LBLstatus->setText(tr("Calculation finished"));
    }
    BTNexec->setEnabled(true);
}


void isosurface::RBTscale_changed(){
    int ix = getscalevalueInt(TXTcontourvalue->text().toFloat(), SLDcontourvalue->minimum(), SLDcontourvalue->maximum(),
                                  mincontourvalue, maxcontourvalue, RBTscalelog->isChecked());
    SLDcontourvalue->setSliderPosition(ix);
    LBLsensitive->setVisible(RBTscalelog->isChecked());
    SPBsensitive->setVisible(RBTscalelog->isChecked());
}

void isosurface::RBTsurftype_changed(){
    if (RBTsolidsurf->isChecked())
        solidsurf = true;
    else
        solidsurf = false;
    emit updatedisplay();
}

void isosurface::SLDcontourvalue_changed(int ix){
    float contour = getscalevalueFloat(ix,SLDcontourvalue->minimum(), SLDcontourvalue->maximum(),
                                  mincontourvalue, maxcontourvalue, RBTscalelog->isChecked());
    TXTcontourvalue->setText(QString::number(contour));
}

void isosurface::SLDcontourvalue_released(){
    int ix = SLDcontourvalue->sliderPosition();
    contourvalue = getscalevalueFloat(ix,SLDcontourvalue->minimum(), SLDcontourvalue->maximum(),
                                  mincontourvalue, maxcontourvalue, RBTscalelog->isChecked());
    TXTcontourvalue->setText(QString::number(contourvalue));
    emit generatesurface();
    emit updatedisplay();
}

void isosurface::SLDopacity_changed(int ix){
    setopacity(float(ix)/100.);
}

void isosurface::SLDopacity_released(){
    opacity = 0.01f * SLDopacity->value();
    emit generatesurface();
    emit updatedisplay();
}

void isosurface::SPBmpi_changed(int nprocessors)
{
    SPBmpi->setValue(nprocessors);
//    qDebug() << "SPBmpi->value() = " << SPBmpi->value();
}

void isosurface::SPBsensitive_changed(int){
    logdlt = SPBsensitive->value();
    TXTcontourvalue_changed();
}

void isosurface::TXTcontourvalue_changed(){
    contourvalue = TXTcontourvalue->text().toFloat();
    int ix = getscalevalueInt(TXTcontourvalue->text().toFloat(), SLDcontourvalue->minimum(), SLDcontourvalue->maximum(),
                              mincontourvalue, maxcontourvalue, RBTscalelog->isChecked());   
    SLDcontourvalue->setTracking(false);
    SLDcontourvalue->setSliderPosition(ix);
    SLDcontourvalue->setTracking(true);
    emit generatesurface();
    emit updatedisplay();
}

//  ------------------------------------------------------------------------------------------------------------------
//
//      General purpose functions
//
//  ------------------------------------------------------------------------------------------------------------------

void isosurface::generategridbounds(float *a){
    float xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = a[0];
    xmax = a[1];
    ymin = a[2];
    ymax = a[3];
    zmin = a[4];
    zmax = a[5];
    gridvertices.clear();
    QVector4D color = QVector4D(1.,0.5,0.,1.);
    VertexNormalData v;
    v.normal.setX(0.);
    v.normal.setY(0.);
    v.normal.setZ(1.);
    v.color.setX(color.x());
    v.color.setY(color.y());
    v.color.setZ(color.z());
    v.color.setW(color.w());
    v.position.setX(xmin);
    v.position.setY(ymin);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymin);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymax);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmin);
    v.position.setY(ymax);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmin);
    v.position.setY(ymin);
    v.position.setZ(zmax);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymin);
    v.position.setZ(zmax);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymax);
    v.position.setZ(zmax);
    gridvertices.append(v);
    v.position.setX(xmin);
    v.position.setY(ymax);
    v.position.setZ(zmax);
    gridvertices.append(v);
    gridindicesoffset.clear();
    gridindices.clear();
    gridindicesoffset << 0 << 2 << 4 << 6 << 8 << 10 << 12 << 14 << 16 << 18 << 20 << 22 << 24;
    gridindices << 0 << 1 << 1 << 2 << 2 << 3 << 3 << 0 << 4 << 5 << 5 << 6 << 6 << 7 << 7 << 4
                << 0 << 4 << 1 << 5 << 2 << 6 << 3 << 7;
}

QString isosurface::get_basename(){
//    QString base = QFileInfo(name).completeBaseName();
    QString base = name;
    int i = base.indexOf("-d_surf");
    if (i < 0){
        isdensity = false;
        i = base.indexOf("-v_surf");
        if (i < 0){
            QMessageBox msgBox;
            msgBox.setText(tr("High quality isosurface"));
            msgBox.setInformativeText(tr("Extension of file %1 must be -d_surf* or -v_surf*. Cannot do operation.").arg(base));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return name;
        }
    }
    base.truncate(i);
    return base;
}

QPoint isosurface::getinitialposition(){
    return initialposition;
}

bool isosurface::getnormalgrad(){
    return normalgrad;
}

bool isosurface::getshowgridbound(){
    return showgridbound;
}

bool isosurface::getsolidsurf(){
    return solidsurf;
}

bool isosurface::gettranslucence(){
    return translucence;
}

bool isosurface::isvisible(){
    return visible;
}

float isosurface::getcontourvalue(){
    return contourvalue;
}

float isosurface::getmaxcontourvalue(){
    return maxcontourvalue;
}

float isosurface::getmincontourvalue(){
    return mincontourvalue;
}

float isosurface::getopacity(){
    return opacity;
}

QColor isosurface::getsurfcolor(){
    return surfcolor;
}

void isosurface::setcontourvalue(float a){
    contourvalue = a;
}

QString isosurface::getfullname(){
    return fullname;
}

QString isosurface::getname(){
    return name;
}


QVector <GLuint> isosurface::getallindices(){
    return allindices;
}

QVector <VertexNormalData> isosurface::getallvertices(){
    return allvertices;
}

void isosurface::setcompatderiv(bool a){
    compatderiv = a;
}

void isosurface::setfullname(QString a){
    fullname = a;
}

void isosurface::setinitialposition(QPoint a){
    initialposition = a;
}

void isosurface::setmaxcontourvalue(float a){
    maxcontourvalue = a;
}

void isosurface::setmincontourvalue(float a){
    mincontourvalue = a;
}

void isosurface::setname(QString a){
    name = a;
    basename = get_basename();
}

void isosurface::setopacity(float a){
    opacity = a;
}

void isosurface::set_ProjectFolder(QString name){
    ProjectFolder = name;
}

void isosurface::set_ProjectName(QString name){
    ProjectName = name;
}

void isosurface::setsolidsurf(bool a){
    solidsurf = a;
    RBTsolidsurf->setChecked(solidsurf);
}

void isosurface::setsurfcolor(QColor a){
    surfcolor = a;
}

void isosurface::settranslucence(bool a){
    translucence = a;
}

void isosurface::setvisible(bool a){
    visible = a;
}


/*******************************************************************************************************/
/********************************  Class editIsoSurfaceDialog  implementation  *******************************/
/*******************************************************************************************************/

editIsoSurfaceDialog::editIsoSurfaceDialog(QWidget *parent) : QDialog(parent)
{

}

editIsoSurfaceDialog::~editIsoSurfaceDialog(){

}

void editIsoSurfaceDialog::reject(){
    emit closed();
}
void editIsoSurfaceDialog::closeEvent(QCloseEvent *event){
    event->ignore();
    this->setVisible(false);
    emit closed();
    event->accept();
}
