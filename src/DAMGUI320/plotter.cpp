//  Copyright 2008-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
//  Guillermo Ramirez, David Zorrilla, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
// 
//  This file is part of DAMQT.
// 
//  DAMQT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  DAMQT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with DAMQT.  If not, see <http://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------
//
//
//  File:   plotter.cpp
//
//      Last version: September 2018
//

#include <cmath>
#include <math.h>       /* round */

#include "plotter.h"

#include <QtDebug>

#include <QtGlobal>

#define ANGSTROMTOBOHR 1.88971616463

Plotter::Plotter(QWidget *parent)
    : QWidget(parent)
{
    arrowsseparation = 150;
    arrowssize = 15;
    arrowsskew = 6;
    arrowswidth = 6;
    aspectratio = true;
    backgroundcolor = Qt::white;
    basinscolor = Qt::red;
    basinspenwidth = 3;
    bonds = new QVector< QVector2D>();
    bondswidth = 2;
    bondsthr = 1.2;
    automaticticks = true;
    bondscolor = QColor(0,0,255);
    bondswidth = 2;
    BottomMargin = 80;
    BTNcolor = nullpointer;
    BTNelinescolor = nullpointer;
    canvassize = size();
    CHKhide = nullpointer;
    CHKplotarrows = nullpointer;
    centerslabelshft = 0;
    centerradius = 3;
    centerscolor = QColor(0,0,0);
    centerslabel = new QVector<QString>();
    centerslabelcolor = QColor(0,0,255);
    centerlabelspos = new QVector<QVector2D>();
    centerlabelrect = new QVector<QRect>();
    cntpenwidth = 2;
    contourlevels = new QVector<double>();
    contourtolerance = 0.9;
    contourplot = false;
    cpscolor << QColor(255,0,0) << QColor(0,255,0) << QColor(195,195,195) << QColor(255,128,0);
    cpsradius = 3;
    curvecolors = new QVector<QColor>();
    curvecontnum = new QVector<int>();
    curvehide = new QVector<bool>();
    curvelabelrect = new QVector<QRect>();
    curvelabels = new QVector<QString>();
    fradplot = false;
    curvespenwidth = 2;
    dialog = nullpointer;
    dialogelines = nullpointer;
    drawarrows = false;
    elinespenwidth = 2;
    efieldcolors = new QVector<QColor>();
    efieldpenwidth = new QVector<int>();
    existcntatlbls = false;
    fieldplot = false;
    fontcontour = QFont("Liberation Mono", 7);
    fontcenterslabel = QFont("Liberation Mono", 20);
    fontcurvelabels = QFont("Liberation Mono", 20);
    fontsghistlabels = QFont("Liberation Mono", 20);
    FRMelinescolor = nullpointer;
    LBLefpenwidth = nullpointer;
    RBTallelines = nullpointer;
    RBTatomelines = nullpointer;
    RBTsingleeline = nullpointer;
    imagesize = QSize(1000,1000);
    imagequality = 20;
    indcenter = -1;
    indcurve = -1;
    indeline = -1;
    invertarrows = false;
    LeftMargin = 80;
    movecurvelabel = false;
    movesghistlabel = false;
    movetitle = false;
    moveXlabel = false;
    moveYlabel = false;
    multicolor = false;
    negativecontourcolor = colorForIds[2];
    negativecontourstyle = 1;
    painterinuse = false;
    Pal = new QPalette(palette());
    planeABC = QVector3D(1.,0.,0.);
    plotarrows = new QVector<bool>();
    plotbasins = false;
    plotbonds = true;
    plotcenters = true;
    plotcps = false;
    positivecontourcolor = colorForIds[0];
    positivecontourstyle = 0;
    rcen = new QVector<QVector3D>();
    RightMargin = 80;
    sghistcolors = new QVector<QColor>();
    sghistpenstyle = new QVector<int>();
//    sghistcontnum = new QVector<int>();
    sghisthaslabel = new QVector<bool>();
    sghisthide = new QVector<bool>();
    sghistlabelrect = new QVector<QRect>();
    sghistlabels = new QVector<QString>();
    sghistplot = false;
    sghistpenwidth = 2;
    showcenterslabel = true;
    showcontourlabels = true;
    showcoords = false;
    showcurvelabels = true;
    showsghistlabels = true;
    showXscalebottom = true;
    showXscaletop = false;
    showYscaleleft = true;
    showYscaleright = false;
    showXticksbottom = true;
    showXtickstop = false;
    showYticksleft = true;
    showYticksright = false;
    showXlabel = true;
    showtitle = false;
    showtitlebox = false;
    showXlabel = false;
    showXlabelbox = false;
    showYlabel = false;
    showYlabelbox = false;
    SPBelinespenwidth = nullpointer;
    surf2d = nullpointer;
    Xcifras = 1;
    Ycifras = 1;
    texto = nullpointer;
    TopMargin = 80;
    transpbackground = false;
    transpbckgrcapture = true;
    updateYpos = false;
    umax = 0.;
    umin = 0.;
    vmax = 0.;
    vmin = 0.;
    uvznuc = new QVector<QVector3D>();
    whoisabove = 2;
    Ylabelvert = true;
    zerocontourcolor =  colorForIds[1];
    zerocontourstyle = 2;
    zerotolerancepow = -5;
    for (int i = 0 ; i < 4 ; i++){
        showcps[i] = true;
    }

    setAutoFillBackground(true);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setFocusPolicy(Qt::StrongFocus);
    rubberBandIsShown = false;

    zoomInButton = new QToolButton(this);
    zoomInButton->setIcon(QIcon(":/images/zoom_in.png"));
    zoomInButton->adjustSize();
    connect(zoomInButton, SIGNAL(clicked()), this, SLOT(zoomIn()));

    zoomOutButton = new QToolButton(this);
    zoomOutButton->setIcon(QIcon(":/images/zoom_out.png"));
    zoomOutButton->adjustSize();
    connect(zoomOutButton, SIGNAL(clicked()), this, SLOT(zoomOut()));

    setPlotSettings(PlotSettings());
}

void Plotter::setPlotSettings(const PlotSettings &settings)
{
    zoomStack.clear();
    zoomStack.append(settings);
    curZoom = 0;
    zoomInButton->hide();
    zoomOutButton->hide();
    update();
}

void Plotter::zoomOut()
{
    if (curZoom > 0) {
        --curZoom;
        zoomOutButton->setEnabled(curZoom > 0);
        zoomInButton->setEnabled(true);
        zoomInButton->show();
        if(fradplot){
            for (int i = 0 ; i < curvelabelrect->count() ; i++ ){ 
                QFontMetrics fm( fontcurvelabels );
                QSize fmsize = fm.size( Qt::TextSingleLine, curvelabels->at(i) );
                (curvelabelrect->operator [](i)).moveCenter(QPoint(4*(canvassize.width())/5,TopMargin+(i+1)*fmsize.height()));
            }
        }
        else if(sghistplot){
            int id = 1;
            for (int i = 0 ; i < sghistlabelrect->count() ; i++ ){
                if (!sghisthaslabel->at(i))
                    continue;
                QFontMetrics fm( fontsghistlabels );
                QSize fmsize = fm.size( Qt::TextSingleLine, sghistlabels->at(i) );
                (sghistlabelrect->operator [](i)).moveCenter(QPoint(4*(canvassize.width())/5,
                                    TopMargin+(id++)*fmsize.height()));
            }
        }
        xmin = zoomStack[curZoom].minX;
        xmax = zoomStack[curZoom].maxX;
        ymin = zoomStack[curZoom].minY;
        ymax = zoomStack[curZoom].maxY;
        LeftMargin = zoomStack[curZoom].LeftMargin;
        RightMargin = zoomStack[curZoom].RightMargin;
        TopMargin = zoomStack[curZoom].TopMargin;
        BottomMargin = zoomStack[curZoom].BottomMargin;
        updateYpos = true;
        showcoords = false;
        zoom_changed();
        emit zoomChanged();
        refreshPixmap();
        updateallcenterslabelsrect();
        updatecontourlabels();

        update();

        emit moveToTop();
    }
}

void Plotter::zoomIn()
{
    if (curZoom < zoomStack.count() - 1) {
        ++curZoom;
        zoomInButton->setEnabled(curZoom < zoomStack.count() - 1);
        zoomOutButton->setEnabled(true);
        zoomOutButton->show();

        if(fradplot){
            for (int i = 0 ; i < curvelabelrect->count() ; i++ ){ 
                QFontMetrics fm( fontcurvelabels );
                QSize fmsize = fm.size( Qt::TextSingleLine, curvelabels->at(i) );
                (curvelabelrect->operator [](i)).moveCenter(QPoint(4*(canvassize.width())/5,TopMargin+(i+1)*fmsize.height()));
            }
        }
        else if(sghistplot){
            int id = 1;
            for (int i = 0 ; i < sghistlabelrect->count() ; i++ ){
                if (!sghisthaslabel->at(i))
                    continue;
                QFontMetrics fm( fontsghistlabels );
                QSize fmsize = fm.size( Qt::TextSingleLine, sghistlabels->at(i) );
                (sghistlabelrect->operator [](i)).moveCenter(QPoint(4*(canvassize.width())/5,TopMargin+(id++)*fmsize.height()));
            }
        }
        xmin = zoomStack[curZoom].minX;
        xmax = zoomStack[curZoom].maxX;
        ymin = zoomStack[curZoom].minY;
        ymax = zoomStack[curZoom].maxY;
        LeftMargin = zoomStack[curZoom].LeftMargin;
        RightMargin = zoomStack[curZoom].RightMargin;
        TopMargin = zoomStack[curZoom].TopMargin;
        BottomMargin = zoomStack[curZoom].BottomMargin;
        updateYpos = true;
        showcoords = false;
        zoom_changed();
        emit zoomChanged();
        refreshPixmap();
        updateallcenterslabelsrect();
        updatecontourlabels();

        update();

        emit moveToTop();
    }
}

void Plotter::zoom_changed(){
    QFontMetrics fm = QFontMetrics( getfontYlabel() );
    QSize fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    setLeftMargin(qMax(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    setRightMargin(qMax(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    fmsize = fm.size( Qt::TextSingleLine, Xlabel );
    setTopMargin(std::max<int>(getTopMargin(),fm.height()+20));
    setBottomMargin(qMax(50,2*fm.size( Qt::TextSingleLine, "0" ).height()+20));
}

void Plotter::updatecontourlabels()
{
    if (contourplot){
        double minx = zoomStack[curZoom].minX;
        double maxx = zoomStack[curZoom].maxX;
        double miny = zoomStack[curZoom].minY;
        double maxy = zoomStack[curZoom].maxY;
        int lm = LeftMargin;
        int tm = TopMargin;
        double dltx = canvassize.width()-lm-RightMargin;
        double dlty = canvassize.height()-tm-BottomMargin;
        for (int i = 0 ; i < contourlabels.length() ; i++){
            ContourLabel clabel = contourlabels[i];
            int x = (clabel.getu()-minx)*dltx/(maxx-minx) + lm;
            int y = (clabel.getv()-maxy)*dlty/(miny-maxy) + tm;
            clabel.setx(x);
            clabel.sety(y);
            QFontMetrics fm = QFontMetrics( fontcontour );
            QSize fmsize = fm.size( Qt::TextSingleLine, clabel.getlabel() );
            clabel.setrect(QRect(clabel.getx() - fmsize.width()/2, clabel.gety() - fmsize.height()/2 ,
                                 fmsize.width(), fmsize.height())) ;
            contourlabels.replace(i,clabel);
        }

    }
}


void Plotter::setContourData(int id, const QVector<QPointF> &data)
{
    contourMap[id] = data;
}

void Plotter::setCurveData(int id, const QVector<QPointF> &data)
{
    curveMap[id] = data;
}

void Plotter::clearCurve(int id, QMap<int, QVector<QPointF> > map)
{
    map.remove(id);
    refreshPixmap();
    update();
}

void Plotter::setEfieldData(int id, const QVector<QPointF> &data)
{
    efieldMap[id] = data;
}


void Plotter::paintEvent(QPaintEvent * /* event */)
{
    QStylePainter painter(this);
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(LeftMargin+1, TopMargin,
        canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);

    if (transpbackground == true)
        painter.setBackgroundMode(Qt::TransparentMode);
    painter.drawPixmap(0, 0, pixmap);
    if (rubberBandIsShown) {
        painter.setPen(QPen(gridcolor));
        painter.drawRect(rubberBandRect.normalized().adjusted(0, 0, -1, -1));
    }

    if (hasFocus()) {
        QStyleOptionFocusRect option;
        option.initFrom(this);
        option.backgroundColor = Qt::white;
        painter.drawPrimitive(QStyle::PE_FrameFocusRect, option);
    }
	
    if (showtitle){
        painter.setFont(fonttitle);
        painter.setPen(fonttitlecolor);
        painter.drawText(titlerect, Qt::AlignCenter, title);
        if (showtitlebox)
            painter.drawRect(titlerect.x()-3,titlerect.y()-3,titlerect.width()+6,titlerect.height()+6);
    }
	
    if (showXlabel){
        painter.setFont(fontXlabel);
        painter.setPen(fontXlabelcolor);
        painter.drawText(Xlabelrect, Qt::AlignCenter, Xlabel);
        if (showXlabelbox)
            painter.drawRect(Xlabelrect.x()-3,Xlabelrect.y()-3,Xlabelrect.width()+6,Xlabelrect.height()+6);
    }
	
    if (showYlabel){
        painter.setFont(fontYlabel);
        painter.setPen(fontYlabelcolor);
        if (Ylabelvert){
                painter.save();
                painter.translate(Ylabelrect.center().x(), Ylabelrect.center().y());
                painter.rotate(270);
                painter.drawText(0, 0, Ylabel);
                painter.restore();
                if (showYlabelbox)
                    painter.drawRect(Ylabelrect.x()-Ylabelrect.width()/4-3, Ylabelrect.y()-Ylabelrect.height()/2-3,
                            Ylabelrect.width()+6, Ylabelrect.height()+6);
        }
        else{
            painter.drawText(Ylabelrect, Qt::AlignCenter, Ylabel);
            if (showYlabelbox)
                painter.drawRect(Ylabelrect.x()-3,Ylabelrect.y()-3,Ylabelrect.width()+6,Ylabelrect.height()+6);
        }
    }
    if (showcontourlabels){
        painter.save();
        painter.setFont(fontcontour);
        painter.setPen(contourslabelcolor);
        for (int i = 0 ; i < contourlabels.count(); ++i){
            if (!rect.contains(contourlabels[i].getrect(),true))
                continue;

            if (!transpbackground == true) painter.setBackgroundMode(Qt::OpaqueMode);
            painter.drawText(contourlabels[i].getrect(), contourlabels[i].getlabel());   
        }
        painter.restore();
    }
    if (showcoords){
        painter.save();
        painter.setBackground(QBrush(Qt::yellow));
        if (!transpbackground == true) painter.setBackgroundMode(Qt::OpaqueMode);
        painter.setPen(Qt::black);
        painter.setFont(fontcontour);
        QFontMetrics fm = QFontMetrics( fontcontour );
        QSize fmsize = fm.size( Qt::TextSingleLine, "x" );
        painter.setBackgroundMode(Qt::TransparentMode);
        painter.drawText(dotx-fmsize.width()/2,doty+fmsize.height()/4, "x");
        if (!transpbackground == true) painter.setBackgroundMode(Qt::OpaqueMode);
        fmsize = fm.size( Qt::TextSingleLine, coordslabel );
        painter.drawText(coordsrect.x()-fmsize.width()/2,coordsrect.y()+2*fmsize.height(), coordslabel);
        painter.restore();
    }
    if (fradplot && showcurvelabels){
        painter.save();
        for (int id = 0 ; id < curveMap.count() ; id++){
            QPen pen;
            if (curvehide->at(id)){
//                pen.setColor(Pal->background().color());
                pen.setColor(Pal->window().color());
            }
            else{
                pen.setColor(curvecolors->at(id));
            }
            painter.setPen(pen);
            painter.setFont(fontcurvelabels);
            painter.drawText(curvelabelrect->at(id), Qt::AlignCenter, curvelabels->at(id));
        }
        painter.restore();
    }
    if (sghistplot && showsghistlabels){
        painter.save();
        for (int id = 0 ; id < curveMap.count() ; id++){
            if (!sghisthaslabel->at(id))
                continue;
            QPen pen;
            if (sghisthide->at(id)){
//                pen.setColor(Pal->background().color());
                pen.setColor(Pal->window().color());
            }
            else{
                pen.setColor(sghistcolors->at(id));
            }
            painter.setPen(pen);
            painter.setFont(fontsghistlabels);
            painter.drawText(sghistlabelrect->at(id), Qt::AlignCenter, sghistlabels->at(id));
        }
        painter.restore();
    }
    if (fieldplot  || (contourplot && existcntatlbls)){
        QPen pen;
        if (plotbonds){
            painter.save();
            pen.setColor(bondscolor);
            pen.setWidth(bondswidth);
            painter.setPen(pen);
            painter.setClipRect(rect);
            QPoint r1;
            QPoint r2;
            for (int i = 0 ; i < bonds->count() ; i +=2){
                r1.setX((double)rect.left() + (((double)bonds->at(i).x() - settings.minX) * ((double)rect.width() ) / settings.spanX()) + 1);
                r1.setY((double)rect.bottom() - (((double)bonds->at(i).y() - settings.minY) * ((double)rect.height() ) / settings.spanY()) + 1);
                r2.setX((double)rect.left() + (((double)bonds->at(i+1).x() - settings.minX) * ((double)rect.width() ) / settings.spanX()) + 1);
                r2.setY((double)rect.bottom() - (((double)bonds->at(i+1).y() - settings.minY) * ((double)rect.height() ) / settings.spanY()) + 1);
                painter.drawLine(r1,r2);
            }
            painter.restore();
        }
        if (plotcenters){
            painter.save();
            pen.setColor(centerslabelcolor);
            painter.setPen(pen);
            painter.setClipRect(rect);
            double x;
            double y;
            painter.setBrush(QBrush(centerscolor));
            for (int i = 0 ; i < uvznuc->count() ; i++){
                x = rect.left() + ((uvznuc->at(i).x() - settings.minX) * (rect.width()) / settings.spanX()) + 1;
                y = rect.bottom() - ((uvznuc->at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
                painter.drawEllipse(QPointF(x,y), centerradius, centerradius);
            }
            painter.restore();
        }
        if (plotcps){
            painter.save();
            painter.setClipRect(rect);
            double x;
            double y;
            painter.setBrush(QBrush(centerscolor));
            for (int i = 0 ; i < cpsuvval.count() ; i++){
                if (!showcps[cpstype[i]])
                    continue;
                painter.setBrush(cpscolor[i]);
                x = rect.left() + ((cpsuvval.at(i).x() - settings.minX) * (rect.width()) / settings.spanX()) + 1;
                y = rect.bottom() - ((cpsuvval.at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
                painter.drawEllipse(QPointF(x,y), cpsradius, cpsradius);
            }
            painter.restore();
        }
        if (showcenterslabel){
            painter.save();
            pen.setColor(centerslabelcolor);
            painter.setFont(fontcenterslabel);
            painter.setPen(pen);
            painter.setClipRect(rect);
            for (int i = 0 ; i < centerlabelrect->count() ; i++){
                painter.drawText(centerlabelrect->at(i), Qt::AlignCenter, centerslabel->at(i));
            }
            painter.restore();
        }
    }
}

void Plotter::capture_image()
{
    QPixmap pixmapcpt;
    QPoint titlecenter = titlerect.center() ;
    QPoint Xlabelcenter = Xlabelrect.center();
    QPoint Ylabelcenter = Ylabelrect.center();

    if (userdefresol){
        canvassize = imagesize;
    }
    else{
        canvassize = size();
    }
    PlotSettings originalsettings = zoomStack[curZoom];
    double scalex = (double)(canvassize.width()) / (double)(width());
    double scaley = (double)(canvassize.height()) / (double)(height());

    zoomStack[curZoom].TopMargin = (double)zoomStack[curZoom].TopMargin * scaley;
    zoomStack[curZoom].BottomMargin = (double)zoomStack[curZoom].BottomMargin * scaley;
    zoomStack[curZoom].LeftMargin = (double)zoomStack[curZoom].LeftMargin * scalex;
    zoomStack[curZoom].RightMargin = (double)zoomStack[curZoom].RightMargin * scalex;

    QRect rect(zoomStack[curZoom].LeftMargin+1, zoomStack[curZoom].TopMargin,
        canvassize.width() - zoomStack[curZoom].LeftMargin - zoomStack[curZoom].RightMargin,
        canvassize.height() - zoomStack[curZoom].TopMargin - zoomStack[curZoom].BottomMargin);

    pixmapcpt = QPixmap(canvassize);
    pixmapcpt.fill(backgroundcolor);
    QPainter painter(&pixmapcpt);
    painter.setClipRegion(QRegion(0,0,canvassize.width(),canvassize.height()));

    if (transpbckgrcapture == true)
        painter.setBackgroundMode(Qt::TransparentMode);
    else
        painter.setBackgroundMode(Qt::OpaqueMode);
    painter.setRenderHint(QPainter::Antialiasing);
//    painter.begin(&image);
    drawGrid(&painter, false, scalex);
    if (fradplot) drawCurves(&painter, scalex, curveMap, 1);
    if (sghistplot) drawCurves(&painter, scalex, curveMap, 4);
    if (whoisabove == 2){
        if (fieldplot) drawCurves(&painter, scalex, efieldMap, 3);
        if (contourplot) drawCurves(&painter, scalex, contourMap, 2);
    }
    else{
        if (contourplot) drawCurves(&painter, scalex, contourMap, 2);
        if (fieldplot) drawCurves(&painter, scalex, efieldMap, 3);
    }
    if (plotbasins){
        drawBasins(&painter, scalex);
    }

    painter.setClipRegion(QRegion(0,0,canvassize.width(),canvassize.height()));
    QPen pen;
    QRect rectaux;
    QPoint originalcenter;
    if (showtitle){
        painter.save();
        pen.setColor(fonttitlecolor);
        painter.setPen(pen);
        QFont font = fonttitle;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        rectaux = titlerect;
        originalcenter = rectaux.center();
        rectaux.setSize(QSize(rectaux.width()*scalex,rectaux.height()*scaley));
        rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
        painter.setFont(font);
        painter.drawText(rectaux, Qt::AlignCenter, title);
        if (showtitlebox) 
            painter.drawRect(rectaux);
        painter.restore();
    }
	
    if (showXlabel){
        painter.save();
        pen.setColor(fontXlabelcolor);
        painter.setPen(pen);
        QFont font = fontXlabel;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        rectaux = Xlabelrect;
        originalcenter = rectaux.center();
        rectaux.setSize(QSize(rectaux.width()*scalex,rectaux.height()*scaley));
        rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
        painter.drawText(rectaux, Qt::AlignCenter, Xlabel);
        if (showXlabelbox) 
            painter.drawRect(rectaux);
        painter.restore();
    }
	
    if (showYlabel){
        painter.save();
        pen.setColor(fontYlabelcolor);
        painter.setPen(pen);
        QFont font = fontYlabel;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        rectaux = Ylabelrect;
        originalcenter = rectaux.center();
        rectaux.setSize(QSize(rectaux.width()*scalex,rectaux.height()*scaley));
        rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
        if (Ylabelvert){
                painter.save();
                painter.translate(rectaux.center().x(), rectaux.center().y());
                painter.rotate(270);
                painter.drawText(0, 0, Ylabel);
                painter.restore();
                if (showYlabelbox) 
                    painter.drawRect(rectaux.x()-rectaux.width()/4, rectaux.y()-rectaux.height()/2,
                            rectaux.width(), rectaux.height());
        }
        else{
            painter.drawText(rectaux, Qt::AlignCenter, Ylabel);
            if (showYlabelbox) 
                painter.drawRect(rectaux);
        }
        painter.restore();
    }
    if (showcontourlabels){
        painter.save();
        pen.setColor(contourslabelcolor);
        painter.setPen(pen);
        QFont font = fontcontour;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        if (transpbackground) {
            painter.setBackgroundMode(Qt::TransparentMode);
        }
        for (int i = 0 ; i < contourlabels.count(); ++i){   
            rectaux = contourlabels[i].getrect();
            originalcenter = rectaux.center();
            rectaux.setSize(QSize(rectaux.width()*1.05*scalex,rectaux.height()*scaley));
            rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
            if (!rect.contains(rectaux,true))
                continue;
            if (!transpbackground )
                painter.fillRect(rectaux, backgroundcolor);
            painter.drawText(rectaux, contourlabels[i].getlabel());
        }
        painter.restore();
    }
    if (showcoords){
        painter.save();
        painter.setBackground(QBrush(Qt::yellow));
        if (!transpbackground == true){
            painter.setBackgroundMode(Qt::OpaqueMode);
        }
        else{
            painter.setBackgroundMode(Qt::TransparentMode);
        }
        painter.setPen(Qt::black);
        QFont font = fontcontour;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        QFontMetrics fm = QFontMetrics( font );
        QSize fmsize = fm.size( Qt::TextSingleLine, "x" );
        painter.drawText(scalex * dotx-fmsize.width()/2, scaley * doty+fmsize.height()/4, "x");
        fmsize = fm.size( Qt::TextSingleLine, coordslabel );
        painter.drawText(coordsrect.x()*scalex-fmsize.width()/2, coordsrect.y()*scaley+2*fmsize.height(), coordslabel);
        painter.restore();
    }
    if (fradplot && showcurvelabels){
        painter.save();
        QFont font = fontcurvelabels;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        for (int i = 0 ; i < curvelabelrect->count() ; i++){
            rectaux = curvelabelrect->at(i);
            originalcenter = rectaux.center();
            rectaux.setSize(QSize(rectaux.width()*scalex,rectaux.height()*scaley));
            rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
            if (curvehide->at(i)){
//                painter.fillRect(rectaux,Pal->background().color());
                painter.fillRect(rectaux,Pal->window().color());
                continue;
            }
            pen.setColor(curvecolors->at(i));
            painter.setPen(pen);
            painter.drawText(rectaux, Qt::AlignCenter, curvelabels->at(i));
        }
        painter.restore();
    }
    else if (sghistplot && showsghistlabels){
        painter.save();
        QFont font = fontsghistlabels;
        font.setPointSizeF(font.pointSizeF()*scalex);
        painter.setFont(font);
        for (int i = 0 ; i < sghistlabelrect->count() ; i++){
            if (!sghisthaslabel->at(i))
                continue;
            rectaux = sghistlabelrect->at(i);
            originalcenter = rectaux.center();
            rectaux.setSize(QSize(rectaux.width()*scalex,rectaux.height()*scaley));
            rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
            if (sghisthide->at(i)){
//                painter.fillRect(rectaux,Pal->background().color());
                painter.fillRect(rectaux,Pal->window().color());
                continue;
            }
            pen.setColor(sghistcolors->at(i));
            painter.setPen(pen);
            painter.drawText(rectaux, Qt::AlignCenter, sghistlabels->at(i));
        }
        painter.restore();
    }
    else if (fieldplot  || (contourplot && existcntatlbls)){

        if (plotbonds){
            painter.save();
            pen.setColor(bondscolor);
            pen.setWidth(bondswidth*scalex);
            painter.setPen(pen);
            PlotSettings settings = zoomStack[curZoom];
            QRect rect(LeftMargin, TopMargin,
                       canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
            QPoint r1;
            QPoint r2;
            for (int i = 0 ; i < bonds->count() ; i +=2){
                r1.setX(rect.left() + ((bonds->at(i).x() - settings.minX) * (rect.width() ) / settings.spanX()) + 1);
                r1.setY(rect.bottom() - ((bonds->at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1);
                r2.setX(rect.left() + ((bonds->at(i+1).x() - settings.minX) * (rect.width() ) / settings.spanX()) + 1);
                r2.setY(rect.bottom() - ((bonds->at(i+1).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1);
                painter.drawLine(r1,r2);
            }
            painter.restore();
        }
        if (plotcenters){
            painter.save();
            pen.setColor(centerslabelcolor);
            painter.setPen(pen);
            PlotSettings settings = zoomStack[curZoom];
            QRect rect(LeftMargin, TopMargin,
                       canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
            double x;
            double y;
            painter.setBrush(QBrush(centerscolor));
            for (int i = 0 ; i < uvznuc->count() ; i++){
                x = rect.left() + ((uvznuc->at(i).x() - settings.minX) * (rect.width() ) / settings.spanX()) + 1;
                y = rect.bottom() - ((uvznuc->at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
                painter.drawEllipse(QPointF(x,y), centerradius*scalex, centerradius*scalex);
            }
            painter.restore();
        }
        if (plotcps){
            painter.save();
            PlotSettings settings = zoomStack[curZoom];
            QRect rect(LeftMargin, TopMargin,
                       canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
            double x;
            double y;
            painter.setBrush(QBrush(centerscolor));
            for (int i = 0 ; i < cpsuvval.count() ; i++){
                painter.setBrush(cpscolor[i]);
                x = rect.left() + ((cpsuvval.at(i).x() - settings.minX) * (rect.width()) / settings.spanX()) + 1;
                y = rect.bottom() - ((cpsuvval.at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
                painter.drawEllipse(QPointF(x,y), cpsradius*scalex, cpsradius*scalex);
            }
            painter.restore();
        }
        if (showcenterslabel){
            painter.save();
            pen.setColor(centerslabelcolor);
            painter.setPen(pen);
            QFont font = fontcenterslabel;
            font.setPointSizeF(font.pointSizeF()*scalex);
            QFontMetrics fm( font );
            painter.setFont(font);
            for (int i = 0 ; i < centerlabelrect->count() ; i++){
                rectaux = centerlabelrect->at(i);
                originalcenter = rectaux.center();
                rectaux.setSize(QSize(rectaux.width()*scalex,rectaux.height()*scaley));
                rectaux.moveCenter(QPoint(originalcenter.x()*scalex, originalcenter.y()*scaley));
                painter.drawText(rectaux, Qt::AlignCenter, centerslabel->at(i));
            }
            painter.restore();
        }
    }
    painter.end();
    QImage image = pixmapcpt.toImage();
    if( !image.save( plotfile, tipo, imagequality) ){
         QMessageBox::warning(this, tr("Saving image file failed"),tr("It could not save the file")+plotfile);
    }
    canvassize = size();
    zoomStack[curZoom] = originalsettings;
    titlerect.moveCenter(titlecenter);
    Xlabelrect.moveCenter(Xlabelcenter);
    Ylabelrect.moveCenter(Ylabelcenter);
    refreshPixmap();
    update();
}

void Plotter::resizeEvent(QResizeEvent * /* event */)
{
    updateplot();
}

void Plotter::updateplot()
{
    canvassize = size();
    titlerect.moveCenter(QPoint(canvassize.width()/2,titlerect.height()));
    emit titleposChanged();
    QFontMetrics fm( getfontYlabel() );
    QSize fmsize;
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    if (Ylabelvert){
//        setLeftMargin(std::max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
        setLeftMargin(fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20);
    }
    else{
//        setLeftMargin(std::max(50,fmsize.width()+20));
        setLeftMargin(fmsize.width()+20);
    }
    setRightMargin(getLeftMargin());
    fm = QFontMetrics( getfonttitle() );
    fmsize = fm.size( Qt::TextSingleLine, title );
    setTopMargin(std::max(50,2*fmsize.height()+20));
    setBottomMargin(std::max(50,2*fm.size( Qt::TextSingleLine, "0" ).height()+20));
    fm = QFontMetrics( fontXlabel );
    fmsize = fm.size( Qt::TextSingleLine, Xlabel );
    Xlabelrect.moveCenter(QPoint(canvassize.width()/2,canvassize.height()-fmsize.height()));
    emit XlabelposChanged();
    fm = QFontMetrics( fontYlabel );
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    Ylabelrect.moveCenter(QPoint(fmsize.height(),canvassize.height()/2));
    emit YlabelposChanged();
//    Ylabelvert = true;
    if (aspectratio){
        double w = canvassize.rwidth();
        double h = canvassize.rheight();
        double deltax = (w - LeftMargin - RightMargin) / (zoomStack[curZoom].maxX - zoomStack[curZoom].minX);
        double deltay = (h - TopMargin - BottomMargin) / (zoomStack[curZoom].maxY - zoomStack[curZoom].minY);
        if (deltax > deltay){ 
            LeftMargin = (int) (0.5*(w - deltay * (zoomStack[curZoom].maxX - zoomStack[curZoom].minX)));  
            RightMargin = LeftMargin;
        }
        else if (deltax < deltay){
            TopMargin = (int) (0.5*(h - deltax * (zoomStack[curZoom].maxY - zoomStack[curZoom].minY)));
            BottomMargin = TopMargin;
        }
        Ylabelrect.moveCenter(QPoint(LeftMargin - 2*Ylabelrect.width() - 20, canvassize.height()/2));
    }
    for (int i = 0 ; i < zoomStack.count() ; ++i){
        zoomStack[i].BottomMargin = std::min(zoomStack[i].BottomMargin,BottomMargin);
        zoomStack[i].TopMargin = std::min(zoomStack[i].TopMargin, TopMargin);
        zoomStack[i].LeftMargin = std::min(zoomStack[i].LeftMargin, LeftMargin);
        zoomStack[i].RightMargin = std::min(zoomStack[i].RightMargin, RightMargin);
    }
    emit MarginsChanged(BottomMargin, TopMargin, LeftMargin, RightMargin);
//    contourlabels.clear();
    updatecontourlabels();
    showcoords = false;
    int x = canvassize.width() - (zoomInButton->width() + zoomOutButton->width() + 10);
    zoomInButton->move(x, 5);
    zoomOutButton->move(x + zoomInButton->width() + 5, 5);
    if (fradplot && showcurvelabels && curvelabelrect != nullpointer){
        for (int i = 0 ; i < curvelabelrect->count() ; ++i){
            (curvelabelrect->operator [](i)).moveCenter(QPoint(4*(canvassize.width())/5,TopMargin+(i+1)*fmsize.height()));
        }
    }
    if (sghistplot && showsghistlabels && sghistlabelrect != nullpointer){
        int id = 1;
        for (int i = 0 ; i < sghistlabelrect->count() ; ++i){
            if (!sghisthaslabel->at(i))
                continue;
            (sghistlabelrect->operator [](i)).moveCenter(QPoint(4*(canvassize.width())/5,TopMargin+(id++)*fmsize.height()));
        }
    }
    if(fieldplot){
        updateallcenterslabelsrect();
    }
    refreshPixmap();
    update();
}

void Plotter::mousePressEvent(QMouseEvent *event)
{
    QRect rect(LeftMargin, TopMargin,
        canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
	const Qt::KeyboardModifiers modifiers = event->modifiers();
    if (event->button() == Qt::LeftButton && modifiers.testFlag(Qt::ShiftModifier)) {
        if (rect.contains(event->pos())) {
            rubberBandIsShown = true;
            rubberBandRect.setTopLeft(event->pos());
            rubberBandRect.setBottomRight(event->pos());
            updateRubberBandRegion();
            setCursor(Qt::CrossCursor);
        }
    }
    if (event->button() == Qt::RightButton && (contourplot || fieldplot)){
        showcoords = true;
        dotx = event->x();
        doty = event->y();
        double u, v;
        u = zoomStack[curZoom].minX + (dotx - LeftMargin) * (zoomStack[curZoom].maxX - zoomStack[curZoom].minX)
                / (canvassize.width() - LeftMargin - RightMargin );
        v = zoomStack[curZoom].maxY + (doty - TopMargin) * (zoomStack[curZoom].minY - zoomStack[curZoom].maxY)
                / (canvassize.height() - TopMargin - BottomMargin );
        double iu = (u - umin) * (nu-1) / (umax - umin);
        double iv = (v - vmin) * (nv-1) / (vmax - vmin);
        if (!(u < umin || u > umax || v < vmin || v > vmax)){
            if (contourplot)
                coordslabel = QString(" (%1, %2, %3) ").arg(u).arg(v).arg(value(iu,iv));
            else if (fieldplot)
                coordslabel = QString(" (%1, %2) ").arg(u).arg(v);
            QFontMetrics fm = QFontMetrics( fontcontour );
            QSize fmsize = fm.size( Qt::TextSingleLine, coordslabel );
            coordsrect = QRect(dotx, doty, fmsize.width(), fmsize.height() );
            repaint();
        }
    }
    else{
        showcoords = false;

        refreshPixmap();
        update();
    }
    if (showtitle && titlerect.contains(event->pos()))
        movetitle = true;
    else
        movetitle = false;
    if (showXlabel && Xlabelrect.contains(event->pos()))
        moveXlabel = true;
    else
        moveXlabel = false;
    if (showYlabel){
        if ((Ylabelvert && Ylabelrect.contains(event->pos()+QPoint(0,Ylabelrect.height()/2))) 
                || (!Ylabelvert && Ylabelrect.contains(event->pos())))
            moveYlabel = true;
    else
            moveYlabel = false;
    }
    movecurvelabel = false;
    if (fradplot && showcurvelabels){
        for (int i = 0 ; i < curveMap.count() ; ++i){
            if (curvelabelrect->at(i).contains(event->pos())){
                movecurvelabel = true;
                indcurve = i;
                break;
            }
        }
    }
    movesghistlabel = false;
    if (sghistplot && showsghistlabels){
        for (int i = 0 ; i < curveMap.count() ; ++i){
            if (sghistlabelrect->at(i).contains(event->pos())){
                movesghistlabel = true;
                indcurve = i;
                break;
            }
        }
    }
    movecenterlabel = false;
    if ((fieldplot  || (contourplot && existcntatlbls)) && showcenterslabel){
        for (int i = 0 ; i < centerlabelrect->count() ; i++){
            if (centerlabelrect->at(i).contains(event->pos())){
                movecenterlabel = true;
                indcenter = i;
                break;
            }
        }
    }
}

void Plotter::mouseMoveEvent(QMouseEvent *event)
{
    if (movetitle ){
        titlerect.moveCenter(event->pos());
        emit titleposChanged();
        update();
    }
    else if (moveXlabel ){
        Xlabelrect.moveCenter(event->pos());
        emit XlabelposChanged();
        update();
    }
    else if (moveYlabel ){
        if (!Ylabelvert){
            Ylabelrect.moveCenter(event->pos());
        }
        else{
            Ylabelrect.moveCenter(event->pos()+QPoint(Ylabelrect.width()/4,Ylabelrect.height()/2));
        }
        emit YlabelposChanged();
        update();
    }
    else if (movecurvelabel){
        (curvelabelrect->operator [](indcurve)).moveCenter(event->pos());
        refreshPixmap();
        update();
    }
    else if (movesghistlabel){
        (sghistlabelrect->operator [](indcurve)).moveCenter(event->pos());
        refreshPixmap();
        update();
    }
    else if (rubberBandIsShown) {
        updateRubberBandRegion();
        rubberBandRect.setBottomRight(event->pos());
        updateRubberBandRegion();
    }
    else if (movecenterlabel){
        (centerlabelrect->operator [](indcenter)).moveCenter(event->pos());
        double x = ((double)(centerlabelrect->at(indcenter).x() - zoomStack[curZoom].LeftMargin - 1) *  zoomStack[curZoom].spanX()) /
                ((double)(canvassize.width() - zoomStack[curZoom].LeftMargin - zoomStack[curZoom].RightMargin - 1)) + zoomStack[curZoom].minX;
        double y = ((double)(canvassize.height() - zoomStack[curZoom].BottomMargin + 1 - centerlabelrect->at(indcenter).y()) *  zoomStack[curZoom].spanY()) /
                ((double)(canvassize.height() - zoomStack[curZoom].TopMargin - zoomStack[curZoom].BottomMargin - 1)) + zoomStack[curZoom].minY;
        centerlabelspos->operator [](indcenter) = QVector2D(x,y);
        updatecenterlabelrect(indcenter);
        update();
    }
}

void Plotter::mouseReleaseEvent(QMouseEvent *event)
{
    if ((event->button() == Qt::LeftButton) && rubberBandIsShown) {
        rubberBandIsShown = false;
        updateRubberBandRegion();
        unsetCursor();

        QRect rect = rubberBandRect.normalized();
        if (rect.width() < 4 || rect.height() < 4)
            return;
        rect.translate(-LeftMargin, -TopMargin);

        PlotSettings prevSettings = zoomStack[curZoom];
        PlotSettings settings;
        double dx = prevSettings.spanX() / (double)(canvassize.width() - LeftMargin - RightMargin);
        double dy = prevSettings.spanY() / (double)(canvassize.height() - TopMargin - BottomMargin);
        settings.minX = prevSettings.minX + dx * rect.left();
        settings.maxX = prevSettings.minX + dx * rect.right();
        settings.minY = prevSettings.maxY - dy * rect.bottom();
        settings.maxY = prevSettings.maxY - dy * rect.top();
//        settings.adjust();
        xmin = settings.minX;
        xmax = settings.maxX;
        ymin = settings.minY;
        ymax = settings.maxY;
        emit zoomChanged();

        zoomStack.resize(curZoom + 1);
        zoomStack.append(settings);
        zoomIn();
		
        movetitle = false;
        moveXlabel = false;
        moveYlabel = false;
        movecurvelabel = false;
        movesghistlabel = false;
    }
    movecenterlabel = false;
    refreshPixmap();
    repaint();
}

/* Defines mouse double click event */
void Plotter::mouseDoubleClickEvent(QMouseEvent *event)
{
    int x = event->x();
    int y = event->y();
    if (contourplot && showcontourlabels){ 
        if(x > LeftMargin && x < canvassize.width()-RightMargin && y > TopMargin && y < canvassize.height()-BottomMargin){
            double u, v;
            u = zoomStack[curZoom].minX + (x - LeftMargin) * (zoomStack[curZoom].maxX - zoomStack[curZoom].minX)
                    / (canvassize.width() - LeftMargin - RightMargin );
            v = zoomStack[curZoom].maxY + (y - TopMargin) * (zoomStack[curZoom].minY - zoomStack[curZoom].maxY)
                    / (canvassize.height() - TopMargin - BottomMargin );
            double iu = (u - umincnt) * (nu-1) / (umaxcnt - umincnt);
            double iv = (v - vmincnt) * (nv-1) / (vmaxcnt - vmincnt);
            double val = value(iu,iv);
            bool addnewcontour = true;
            for (int j = 0 ; j < contourlabels.count(); ++j){
                if (contourlabels[j].getrect().contains(QPoint(x,y))){
                    contourlabels.remove(j);
                    addnewcontour = false;
                    break;
                }
            }
            if (addnewcontour){        
                bool lzero = false;
                int ilevelzero = -1;
                for (int i = 0 ; i < nlevels ; ++i){
                    if (contourlevels->at(i) == zero && !lzero){
                        ilevelzero = i;
                        lzero = true;
                        continue;
                    } 
                    if (contourlevels->at(i)*val <= zero) continue;
                    double max = std::max<double>(std::abs(contourlevels->at(i)),std::abs(val));
                    double min = std::min<double>(std::abs(contourlevels->at(i)),std::abs(val));
                    if ((max > 10 && std::log10(min)/std::log10(max) > contourtolerance)
                            || (min/max) > contourtolerance){
                        lzero = false;
                        const Qt::KeyboardModifiers modifiers = event->modifiers();
                        if(multicolor && modifiers.testFlag(Qt::ShiftModifier)){
                            contourlabels_dialog(i);
                        }
                        else{
                            ContourLabel clabel;
                            QFontMetrics fm = QFontMetrics( fontcontour );
                            clabel.setlabel(QString("%1").arg(contourlevels->at(i)));
                            clabel.setu(u);
                            clabel.setv(v);
                            clabel.setx(x);
                            clabel.sety(y);
                            QSize fmsize = fm.size( Qt::TextSingleLine, clabel.getlabel() );
                            clabel.setrect(QRect(x - fmsize.width()/2, y - fmsize.height()/2 , fmsize.width(), fmsize.height())) ;
                            contourlabels.append(clabel);
                        }
                        break;
                    }
                }
                if (lzero && ilevelzero >= 0){
                    if (std::abs(val) < std::pow(10.,zerotolerancepow)){
                        const Qt::KeyboardModifiers modifiers = event->modifiers();
                        if(modifiers.testFlag(Qt::ShiftModifier)){
                            contourlabels_dialog(ilevelzero);
                        }
                        else{
                            ContourLabel clabel;
                            QFontMetrics fm = QFontMetrics( fontcontour );
                            clabel.setlabel(QString("0"));
                            clabel.setu(u);
                            clabel.setv(v);
                            clabel.setx(x);
                            clabel.sety(y);
                            QSize fmsize = fm.size( Qt::TextSingleLine, clabel.getlabel() );
                            clabel.setrect(QRect(x - fmsize.width()/2, y - fmsize.height()/2 , fmsize.width(), fmsize.height())) ;
                            contourlabels.append(clabel);
                        }
                    }
                }
            }
            update();
        }
    }
    if (fradplot && showcurvelabels && curvelabelrect != nullpointer){
        for (int i = 0 ; i < curvelabelrect->count() ; ++i){
            if (curvelabelrect->at(i).contains(QPoint(x,y))){
                curvelabels_dialog(i);
                break;
            }
        }
    }
    else if (sghistplot && showsghistlabels && sghistlabelrect != nullpointer){
        for (int i = 0 ; i < sghistlabelrect->count() ; ++i){
            if (!sghisthaslabel->at(i))
                continue;
            if (sghistlabelrect->at(i).contains(QPoint(x,y))){
                sghistlabels_dialog(i);
                break;
            }
        }
    }
    else if (fieldplot){
        if(x > LeftMargin && x < canvassize.width()-RightMargin && y > TopMargin && y < canvassize.height()-BottomMargin){
            int iline = elines(x, y);
            if (iline >= 0){
                elines_dialog(iline);
            }
        }
    }
    if (showtitle && titlerect.contains(QPoint(x,y))){
        emit title_dialog();
    } 
    if (showXlabel && Xlabelrect.contains(QPoint(x,y))){
        emit Xlabel_dialog();
    }
    if (showYlabel && ((Ylabelvert && Ylabelrect.contains(QPoint(x,y)
            +QPoint(Ylabelrect.width()/4,Ylabelrect.height()/2))) 
                || (!Ylabelvert && Ylabelrect.contains(QPoint(x,y))))){
        emit Ylabel_dialog();
    }
    
}

/* Sets a region for zooming */
void Plotter::setzoomRegion()
{
    bool update = false;    
    PlotSettings prevSettings = zoomStack[curZoom];
    PlotSettings settings;
    if (std::abs(xmin-prevSettings.minX) > 1.e-5 && xmin < xmax){
        settings.minX = xmin;
        update = true;
    }
    else
        settings.minX = prevSettings.minX;
    if (std::abs(xmax-prevSettings.maxX) > 1.e-5 && xmin < xmax){
        settings.maxX = xmax;
        update = true;
    }
    else
        settings.maxX = prevSettings.maxX;
    if (std::abs(ymin-prevSettings.minY) > 1.e-3*std::abs(ymin) && ymin < ymax){
        settings.minY = ymin;
        update = true;
    }
    else
        settings.minY = prevSettings.minY;
    if (std::abs(ymax-prevSettings.maxY) > 1.e-3*std::abs(ymax) && ymin < ymax){
        settings.maxY = ymax;
        update = true;
    }
    else
        settings.maxY = prevSettings.maxY;
    settings.BottomMargin = prevSettings.BottomMargin;
    settings.TopMargin = prevSettings.TopMargin;
    settings.LeftMargin = prevSettings.LeftMargin;
    settings.RightMargin = prevSettings.RightMargin;
    if (!update)
        return;
//    settings.adjust();
    zoomStack.resize(curZoom + 1);
    zoomStack.append(settings);
    zoomIn();
}

double Plotter::value(double x, double y)
{
    int nx, ny, nx1, ny1;
    double f11, f12, f21, f22, x1, y1, x2, y2;
    nx = ((int) x);
    ny = ((int) y);
    if (nx >= nu || nx < 0 || ny >= nv || ny < 0)
        return -99999.;
    x1 = (double) (nx);
    y1 = (double) (ny);
    f11 = surf2d->at(ny*nu+nx);
    if (nx < nu-1){
        nx1 = 1;
        x2 = x1 + 1.;
    }
    else{
        nx1 = -1;
        x2 = x1 - 1.;
    }
    if (ny < nv-1){
        ny1 = 1;
        y2 = y1 + 1.;
    }
    else{
        ny1 = -1;
        y2 = y1 - 1.;
    }
    f12 = surf2d->at((ny+ny1)*nu+nx);
    f21 = surf2d->at(ny*nu+nx+nx1);
    f22 = surf2d->at((ny+ny1)*nu+nx+nx1);
    return (f11 * (x2-x) * (y2-y) + f21 * (x-x1) * (y2-y)
            + f12 * (x2-x) * (y-y1)  + f22 * (x-x1) * (y-y1))
            / ((x2-x1) * (y2-y1));
}

void Plotter::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Plus:
        zoomIn();
        break;
    case Qt::Key_Minus:
        zoomOut();
        break;
    case Qt::Key_Left:
        zoomStack[curZoom].scroll(+1, 0);
        refreshPixmap();
        updateallcenterslabelsrect();
        update();
        break;
    case Qt::Key_Right:
        zoomStack[curZoom].scroll(-1, 0);
        refreshPixmap();
        updateallcenterslabelsrect();
        update();
        break;
    case Qt::Key_Down:
        zoomStack[curZoom].scroll(0, +1);
        refreshPixmap();
        updateallcenterslabelsrect();
        update();
        break;
    case Qt::Key_Up:
        zoomStack[curZoom].scroll(0, -1);
        refreshPixmap();
        updateallcenterslabelsrect();
        update();
        break;
    default:
        QWidget::keyPressEvent(event);
    }
}

void Plotter::wheelEvent(QWheelEvent *event)
{

}

void Plotter::updateRubberBandRegion()
{
    QRect rect = rubberBandRect.normalized();
    update(rect.left(), rect.top(), rect.width(), 1);
    update(rect.left(), rect.top(), 1, rect.height());
    update(rect.left(), rect.bottom(), rect.width(), 1);
    update(rect.right(), rect.top(), 1, rect.height());
}

void Plotter::refreshPixmap()
{
    if (painterinuse)  // To prevent attempts to paint by LeftMargins_Changed called by drawGrid (painter is in use) 
        return;
    pixmap = QPixmap(canvassize);
    pixmap.fill(backgroundcolor);
    QPainter painter(&pixmap);
//    painter.initFrom(this);
    painter.begin(this);
    painterinuse = true;
    drawGrid(&painter, true, 1.);
    if (fradplot) drawCurves(&painter, 1., curveMap, 1);
    if (sghistplot) drawCurves(&painter, 1., curveMap, 4);
    if (whoisabove == 2){
        if (fieldplot) drawCurves(&painter, 1., efieldMap, 3);
        if (contourplot) drawCurves(&painter, 1., contourMap, 2);
    }
    else{
        if (contourplot) drawCurves(&painter, 1., contourMap, 2);
        if (fieldplot) drawCurves(&painter, 1., efieldMap, 3);
    }
    if (plotbasins){
        drawBasins(&painter, 1.);
    }
    updateallcenterslabelsrect();
    painter.end();
    painterinuse = false;
}

void Plotter::drawGrid(QPainter *painter, bool display, double scalefont)
{
    PlotSettings settings = zoomStack[curZoom];
    QFont font;
    zoomStack[curZoom].ncifrasy = std::max<int>(0,std::floor(-std::log10(std::abs(settings.spanY() / settings.numYTicks))));
    if (((double)zoomStack[curZoom].ncifrasy + std::floor(std::log10(std::abs(settings.spanY() / settings.numYTicks)))) < 0)
        zoomStack[curZoom].ncifrasy++;
    if (!automaticticks)
        zoomStack[curZoom].ncifrasy = std::max<int>(zoomStack[curZoom].ncifrasy,Ycifras);
    zoomStack[curZoom].ncifrasx = std::max<int>(0,std::floor(-std::log10(std::abs(settings.spanX() / settings.numXTicks))));
    if (((double)zoomStack[curZoom].ncifrasx + std::floor(std::log10(std::abs(settings.spanX() / settings.numXTicks)))) < 0)
        zoomStack[curZoom].ncifrasx++;
    char f, fx, fy;
    fx = 'f';
    if (std::max(std::abs(settings.minX),std::abs(settings.minX+settings.spanX())) < 0.1){
        fx = 'e';
        zoomStack[curZoom].ncifrasx = Xcifras;
    }
    else if(zoomStack[curZoom].ncifrasx < Xcifras){
        zoomStack[curZoom].ncifrasx = Xcifras;
    }
    BottomMargin = zoomStack[curZoom].BottomMargin;
    TopMargin = zoomStack[curZoom].TopMargin;
    fy = 'f';
    if (std::max(std::abs(settings.minY),std::abs(settings.minY+settings.spanY())) < 0.1 ||
           std::max(std::abs(settings.minY),std::abs(settings.minY+settings.spanY())) > 1000. ){
        fy = 'e';
        zoomStack[curZoom].ncifrasy = Ycifras;
    }
    else if(zoomStack[curZoom].ncifrasy < Ycifras){
        zoomStack[curZoom].ncifrasy = Ycifras;
    }
    LeftMargin = zoomStack[curZoom].LeftMargin;
#if QT_VERSION < 0x050B00
    LeftMargin = std::max<int>(LeftMargin,int(5 + 2*QFontMetrics(fontYlabel).xHeight() + QFontMetrics(fontYlabel).width(QString("0"))
            *(QString::number(settings.minY+settings.spanY(),fy,zoomStack[curZoom].ncifrasy).length()+1)));
#else
    LeftMargin = std::max<int>(LeftMargin,int(5 + 2*QFontMetrics(fontYlabel).xHeight() + QFontMetrics(fontYlabel).horizontalAdvance(QString("0"))
            *(QString::number(settings.minY+settings.spanY(),fy,zoomStack[curZoom].ncifrasy).length()+1)));
#endif
    LeftMargin = std::max<int>(LeftMargin,int(5 + 2.5*QFontMetrics(fontYlabel).height()));
    RightMargin = zoomStack[curZoom].RightMargin;
    if (aspectratio){
        double w = canvassize.rwidth();
        double h = canvassize.rheight();
        double deltax = (w - LeftMargin - RightMargin) / (zoomStack[curZoom].maxX - zoomStack[curZoom].minX);
        double deltay = (h - TopMargin - BottomMargin) / (zoomStack[curZoom].maxY - zoomStack[curZoom].minY);
        if (deltax > deltay){ 
            LeftMargin = (int) (0.5*(w - deltay * (zoomStack[curZoom].maxX - zoomStack[curZoom].minX)));  
            RightMargin = LeftMargin;
        }
        else if (deltax < deltay){
            TopMargin = (int) (0.5*(h - deltax * (zoomStack[curZoom].maxY - zoomStack[curZoom].minY)));
            BottomMargin = TopMargin;
        }
        if (updateYpos){
            Ylabelrect.moveCenter(QPoint(LeftMargin - 2*Ylabelrect.width() - 20,
                                         (canvassize.height()+Ylabelrect.height())/2));
            updateYpos = false;
        }
    }
    if (display) emit MarginsChanged(BottomMargin, TopMargin, LeftMargin, RightMargin);
    zoomStack[curZoom].BottomMargin = BottomMargin;
    zoomStack[curZoom].TopMargin = TopMargin;
    zoomStack[curZoom].LeftMargin = LeftMargin;
    zoomStack[curZoom].RightMargin = RightMargin;

    settings = zoomStack[curZoom];

    QRect rect(LeftMargin+1, TopMargin,
               canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);

    if (!rect.isValid())
        return;
    
    QPen gridpen = QPen(gridcolor);
    if (solidgrid)
        gridpen.setStyle(Qt::SolidLine);
    else
        gridpen.setStyle(Qt::DashLine);
    zoomStack[curZoom].ncifrasx = std::max<int>(zoomStack[curZoom].ncifrasx,Xcifras);
    zoomStack[curZoom].ncifrasy = std::max<int>(zoomStack[curZoom].ncifrasy,Ycifras);
    for (int i = 0; i <= settings.numXTicks; ++i) {
        int x = qMin(rect.right(),rect.left() + (i * (rect.width() ) / settings.numXTicks));
        double label = settings.minX + (i * settings.spanX()
                                          / settings.numXTicks);
        painter->setPen(gridpen);
        if (showgrid && i < settings.numXTicks)
                painter->drawLine(x, rect.top(), x, rect.bottom());
        if (showXticksbottom)
                painter->drawLine(x, rect.bottom(), x, rect.bottom() + 5*scalefont);
        if (showXtickstop)
                painter->drawLine(x, rect.top()-5*scalefont, x, rect.top());
        if (label == 0)
            f = 'f';
        else
            f = fx;
        font = fontXlabel;
        font.setPointSizeF(font.pointSizeF()*scalefont);
        painter->setFont(font);
        painter->setPen(fontXlabelcolor);
	
        if (showXscalebottom)
                painter->drawText(x - QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasx).length()/2,
                rect.bottom() + 5*scalefont,
                                QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasx).length(),
                QFontMetrics(font).height(), Qt::AlignHCenter | Qt::AlignTop, QString::number(label,f,zoomStack[curZoom].ncifrasx));
        if (showXscaletop)
                painter->drawText(x - QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasx).length()/2,
                rect.top() - 5*scalefont - QFontMetrics(font).height(),
                                QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasx).length(),
                QFontMetrics(font).height(), Qt::AlignHCenter | Qt::AlignTop, QString::number(label,f,zoomStack[curZoom].ncifrasx));
    }

    for (int j = 0; j <= settings.numYTicks; ++j) {
        int y = rect.bottom() - (j * (rect.height() - 1)
                                   / settings.numYTicks);
        double label = settings.minY + (j * settings.spanY()
                                          / settings.numYTicks);
        painter->setPen(gridpen);
        if (showgrid)
                painter->drawLine(rect.left(), y, rect.right(), y);
        if (showYticksleft)
                painter->drawLine(rect.left() - 5*scalefont, y, rect.left(), y);
        if (showYticksright)
                painter->drawLine(rect.right(), y, rect.right()+5*scalefont, y);
		
        font = fontYlabel;
        font.setPointSizeF(font.pointSizeF()*scalefont);
        painter->setFont(font);
        painter->setPen(fontYlabelcolor);
        if ((int)label == 0)
            f = 'f';
        else
            f = fy;
        if (showYscaleleft)	
                painter->drawText(rect.left() - QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasy).length() - 7*scalefont,
                y - QFontMetrics(font).height()/2,
                QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasy).length(), QFontMetrics(font).height(),
				Qt::AlignRight | Qt::AlignVCenter, QString::number(label,f,zoomStack[curZoom].ncifrasy));
        if (showYscaleright)	
                painter->drawText(rect.right() + 8*scalefont,
                y - QFontMetrics(font).height()/2,
                QFontMetrics(font).maxWidth()*QString::number(label,f,zoomStack[curZoom].ncifrasy).length(), QFontMetrics(font).height(),
                Qt::AlignLeft | Qt::AlignVCenter, QString::number(label,f,zoomStack[curZoom].ncifrasy));
    }
    gridpen.setStyle(Qt::SolidLine);
    painter->setPen(gridpen);
    painter->drawRect(rect.adjusted(0, 0, -1, -1));
    xmax = zoomStack[curZoom].maxX;
    xmin = zoomStack[curZoom].minX;
    ymax = zoomStack[curZoom].maxY;
    ymin = zoomStack[curZoom].minY;
}

void Plotter::drawCurves(QPainter *painter, double scalewidth, QMap<int, QVector<QPointF> > map, int caseplot)
{
//    caseplot == 1: 1D curves (radial factors and derivatives)
//    caseplot == 2: 2D contour plots
//    caseplot == 3: 2D field lines
//    caseplot == 4: 1D MESP sigma hole histogram
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(LeftMargin+1, TopMargin,
               canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
    if (!rect.isValid())
        return;

    painter->setClipRect(rect.adjusted(+1, +1, -1, -1));
    QMapIterator<int, QVector<QPointF> > i(map);
    QVector<QPointF> *arrowsline = new QVector<QPointF>();
    QPen pen;
    while (i.hasNext()) {
        i.next();
        int id = i.key();
        bool showarrows;
        if (plotarrows->size() > id)
            showarrows = plotarrows->at(id);
        else
            showarrows = false;
        QVector<QPointF> data = i.value();
        QPolygonF polyline(data.count());
        arrowsline->clear();
        for (int j = 0; j < data.count(); ++j) {
            double dx = data[j].x() - settings.minX;
            double dy = data[j].y() - settings.minY;
            double x = rect.left() + (dx * (rect.width() ) / settings.spanX());
            double y = rect.bottom() - (dy * (rect.height() ) / settings.spanY())+1;
            polyline[j] = QPointF(x, y);
            if ( (caseplot == 3) && drawarrows && showarrows && (j+1)%(arrowsseparation+arrowsskew) == 0 ){
                QVector2D r1(data[j]);
                QVector2D r2(data[j-arrowsskew+1]);
                if (!invertarrows){
                    arrowsline->append(QPointF(r1.toPointF()));
                    arrowsline->append(QPointF(r2.toPointF()));
                }
                else{
                    arrowsline->append(QPointF(r2.toPointF()));
                    arrowsline->append(QPointF(r1.toPointF()));
                }
            }
        }
        if (caseplot == 1){     // 1D curves (radial factors and derivatives)
            pen = QPen(curvecolors->at(id));
            pen.setWidth(curvespenwidth * scalewidth);
//          In the following if block the label is only written in the 2D viewer (not in the capture)
            if ((std::abs(scalewidth-1.) < 0.009) && showcurvelabels && curvelabelrect){
                painter->save();
                painter->setFont(fontcurvelabels);
                painter->setPen(pen);
                painter->drawText(curvelabelrect->at(id), Qt::AlignCenter, curvelabels->at(id));
                painter->restore();
            }
        }
        if (caseplot == 2){     // 2D contour plots
            if (multicolor){
                pen = QPen(curvecolors->at(curvecontnum->at(id)));
                if (contourplot && contourlevels){
                    if (contourlevels->at(curvecontnum->at(id)) > 0){
                        pen.setStyle(Qt::PenStyle (penstyle[positivecontourstyle]) );
                    }
                    else if (contourlevels->at(curvecontnum->at(id)) == zero){
                        pen.setStyle(Qt::PenStyle (penstyle[zerocontourstyle]) );
                    }
                    else{
                        pen.setStyle(Qt::PenStyle (penstyle[negativecontourstyle]) );
                    }
                }
            }
            else if (curvecontnum->count() > 0){
                if (contourlevels->at(curvecontnum->at(id)) > zero){
                    pen = QPen(positivecontourcolor);
                    pen.setStyle(Qt::PenStyle (penstyle[positivecontourstyle]) );
                }
                else if (contourlevels->at(curvecontnum->at(id)) == zero){
                    pen = QPen(zerocontourcolor);
                    pen.setStyle(Qt::PenStyle (penstyle[zerocontourstyle]) );
                }
                else{
                    pen = QPen(negativecontourcolor);
                    pen.setStyle(Qt::PenStyle (penstyle[negativecontourstyle]) );
                }
            }
            pen.setWidth(cntpenwidth * scalewidth);
        }
        else if (caseplot == 3){    // 2D field lines
            pen = QPen(efieldcolors->at(id));
            if (drawarrows){
                double x1;
                double y1;
                double x2;
                double y2;
                double x3;
                double y3;
                double xaux;
                double yaux;
                for (int k = 0 ; k < arrowsline->count()-1 ; k += 2 ){
                    QVector2D r1(arrowsline->at(k));
                    QVector2D raux(arrowsline->at(k+1));
                    x1 = rect.left() + ((r1.x() - settings.minX) * (rect.width() ) / settings.spanX());
                    y1 = rect.bottom() - ((r1.y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
                    xaux = rect.left() + ((raux.x() - settings.minX) * (rect.width() ) / settings.spanX());
                    yaux = rect.bottom() - ((raux.y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
                    double norm = (QVector2D(x1,y1)-QVector2D(xaux,yaux)).length();
                    if (norm != 0.){
                        double sint = (y1-yaux) / norm;
                        double cost = (x1-xaux) / norm;
                        xaux = x1 - arrowssize * cost * scalewidth;
                        yaux = y1 - arrowssize * sint * scalewidth;
                        x2 = xaux - arrowswidth * sint * scalewidth;
                        y2 = yaux + arrowswidth * cost * scalewidth;
                        x3 = xaux + arrowswidth * sint * scalewidth;
                        y3 = yaux - arrowswidth * cost * scalewidth;
                        QPainterPath triangle;
                        triangle.moveTo(x1,y1);
                        triangle.lineTo(x2,y2);
                        triangle.lineTo(x3,y3);
                        triangle.lineTo(x1,y1);
                        painter->setPen (Qt :: NoPen);
                        painter->fillPath (triangle, QBrush (efieldcolors->at(id)));
                    }
                }
            }
            pen.setWidth(efieldpenwidth->at(id) * scalewidth);
        }
        else if (caseplot == 4){     // 1D curves (sigma hole histogram)
            pen = QPen(sghistcolors->at(id));
            pen.setStyle(Qt::PenStyle (penstyle[sghistpenstyle->at(id)]) );
            pen.setWidth(sghistpenwidth * scalewidth);
//          In the following if block the label is only written in the 2D viewer (not in the capture)
            if ((std::abs(scalewidth-1.) < 0.009) && showsghistlabels && sghistlabelrect && sghisthaslabel->at(id)){
                painter->save();
                painter->setFont(fontsghistlabels);
                painter->setPen(pen);
                painter->drawText(sghistlabelrect->at(id), Qt::AlignCenter, sghistlabels->at(id));
                painter->restore();
            }
        }
        painter->setPen(pen);
        if ((fradplot && curvehide->at(id)) || (sghistplot && sghisthide->at(id))) continue;
        painter->drawPolyline(polyline);
    }
}

void Plotter::drawBasins(QPainter *painter, double scalewidth){
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(LeftMargin+1, TopMargin,
               canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
    if (!rect.isValid())
        return;

    painter->setClipRect(rect.adjusted(+1, +1, -1, -1));
    QPen pen = QPen(basinscolor);
    pen.setWidth(basinspenwidth * scalewidth);
    painter->setPen(pen);
    for (int i = 0 ; i < basinsuv.count() ; i++){
        QPolygonF polyline(basinsuv.at(i).count());
        for (int j = 0; j < basinsuv.at(i).count(); ++j) {
            double dx = basinsuv.at(i).at(j).x() - settings.minX;
            double dy = basinsuv.at(i).at(j).y() - settings.minY;
            double x = rect.left() + (dx * (rect.width() ) / settings.spanX());
            double y = rect.bottom() - (dy * (rect.height() ) / settings.spanY())+1;
            polyline[j] = QPointF(x, y);
        }
        painter->drawPolyline(polyline);
    }
}

int Plotter::elines(int x0, int y0)
{
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(LeftMargin, TopMargin,
               canvassize.width() - LeftMargin - RightMargin, canvassize.height() - TopMargin - BottomMargin);
    if (!rect.isValid())
        return -1;
    QMapIterator<int, QVector<QPointF> > i(efieldMap);
    while (i.hasNext()) {
        i.next();
        QVector<QPointF> data = i.value();
        for (int j = 0; j < data.count(); ++j) {
            double dx = data[j].x() - settings.minX;
            double dy = data[j].y() - settings.minY;
            int x = qRound(rect.left() + (dx * (rect.width() ) / settings.spanX()));
            int y = rect.bottom() - (dy * (rect.height() ) / settings.spanY());
            if (qAbs(x0-x) < 3 && qAbs(y0-y) < 3){
                return i.key();
            }
        }
    }
    return -1;
}

void Plotter::Xticks_changed(int numx){
    if (numx < 2) return;
    zoomStack[curZoom].numXTicks = numx-1;
    zoomStack[curZoom].ncifrasx = Xcifras;
    zoomStack[curZoom].adjustAxisforced(zoomStack[curZoom].minX,zoomStack[curZoom].maxX,zoomStack[curZoom].numXTicks);
    refreshPixmap();
    update();
}

void Plotter::Yticks_changed(int numy){
    if (numy < 2) return;
    zoomStack[curZoom].numYTicks = numy-1;
    zoomStack[curZoom].ncifrasy = Ycifras;
    zoomStack[curZoom].adjustAxisforced(zoomStack[curZoom].minY,zoomStack[curZoom].maxY,zoomStack[curZoom].numYTicks);
    refreshPixmap();
    update();
}

void Plotter::restore_automatic(){
    for (int i = 0 ; i < zoomStack.count() ; ++i){
        zoomStack[i].adjustforced();
    }
}

void Plotter::aceptar_label()
{
    if (indcurve < 0) return;
    curvelabels->replace(indcurve,texto->text());
    QFontMetrics fm( fontcurvelabels );
    (curvelabelrect->replace(indcurve,QRect(curvelabelrect->at(indcurve).topLeft(),
            fm.size( Qt::TextSingleLine, curvelabels->at(indcurve)))));
    curvecolors->replace(indcurve,color);
    cerrar_label();
    refreshPixmap();
    update();
}

void Plotter::aceptar_sghistlabel()
{
    if (indcurve < 0) return;
    sghistlabels->replace(indcurve,texto->text());
    QFontMetrics fm( fontsghistlabels );
    (sghistlabelrect->replace(indcurve,QRect(sghistlabelrect->at(indcurve).topLeft(),
            fm.size( Qt::TextSingleLine, sghistlabels->at(indcurve)))));
    sghistcolors->replace(indcurve,color);
    for (int i = indcurve+1 ; i < sghistcolors->count() ; i++){
        if (sghisthaslabel->at(i))
            break;
        sghistcolors->replace(i,color);
    }
    cerrar_label();
    refreshPixmap();
    update();
}

void Plotter::cerrar_label()
{
    if (texto != nullpointer){
        delete(texto);
        texto = nullpointer;
    }
    if (BTNcolor != nullpointer){
        delete(BTNcolor);
        BTNcolor = nullpointer;
    }
    if (CHKhide != nullpointer){
        delete(CHKhide);
        CHKhide = nullpointer;
    }
    if (dialog != nullpointer){
        delete(dialog);
        dialog = nullpointer;
    }
}

void Plotter::setcenterlabelrect(){
    double x;
    double y;
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(settings.LeftMargin, settings.TopMargin,
               canvassize.width() - settings.LeftMargin - settings.RightMargin, canvassize.height() - settings.TopMargin - settings.BottomMargin);
    if (!rect.isValid())
        return;
    centerlabelrect->clear();
    centerlabelspos->clear();
    for (int i = 0 ; i < centerslabel->count() ; i++){
        centerlabelspos->append(QVector2D(uvznuc->at(i).x(),uvznuc->at(i).y()));
        x = rect.left() + ((uvznuc->at(i).x() - settings.minX) * (rect.width() ) / settings.spanX()) + 1;
        y = rect.bottom() - ((uvznuc->at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
        QFontMetrics fm( getfontcenterslabel() );
        QSize fmsize = fm.size( Qt::TextSingleLine, centerslabel->at(i) );
        centerlabelrect->append(QRect(0,0, fmsize.width(), fmsize.height() ));
        centerlabelrect->last().moveCenter(QPoint(x,y+0.5*fmsize.height()-centerslabelshft));
    }
}

void Plotter::updateallcenterslabelsrect(){
    double x;
    double y;
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(settings.LeftMargin, settings.TopMargin,
               canvassize.width() - settings.LeftMargin - settings.RightMargin, canvassize.height() - settings.TopMargin - settings.BottomMargin);
    if (!rect.isValid())
        return;
    centerlabelrect->clear();
    for (int i = 0 ; i < centerslabel->count() ; i++){
        x = rect.left() + ((centerlabelspos->at(i).x() - settings.minX) * (rect.width() ) / settings.spanX()) + 1;
        y = rect.bottom() - ((centerlabelspos->at(i).y() - settings.minY) * (rect.height() ) / settings.spanY()) + 1;
        QFontMetrics fm( getfontcenterslabel() );
        QSize fmsize = fm.size( Qt::TextSingleLine, centerslabel->at(i) );
        centerlabelrect->append(QRect(0,0, fmsize.width(), fmsize.height() ));
        centerlabelrect->last().moveCenter(QPoint(x,y+0.5*fmsize.height()-centerslabelshft));
    }
}

void Plotter::updatecenterlabelrect(int index){
    PlotSettings settings = zoomStack[curZoom];
    QRect rect(settings.LeftMargin, settings.TopMargin,
               canvassize.width() - settings.LeftMargin - settings.RightMargin, canvassize.height() - settings.TopMargin - settings.BottomMargin);
    if (!rect.isValid())
        return;
    double x = (double)rect.left() + ((centerlabelspos->at(index).x() - settings.minX) * ((double)rect.width() ) / settings.spanX()) + 1;
    double y = (double)rect.bottom() - ((centerlabelspos->at(index).y() - settings.minY) * ((double)rect.height() ) / settings.spanY()) + 1;
    QFontMetrics fm( getfontcenterslabel() );
    QSize fmsize = fm.size( Qt::TextSingleLine, centerslabel->at(index) );
    centerlabelrect->operator [](index).moveCenter(QPoint(x,y+0.5*fmsize.height()-centerslabelshft));
}


void Plotter::BTNcolor_clicked(){
    QColor col = QColorDialog::getColor(color, this);
    if(col.isValid()) {
        color = col;
        BTNcolor->setColor(&color);
    }
}

void Plotter::BTNelinescolor_clicked(){
    QColor col = QColorDialog::getColor(elinescolor, this);
    if(col.isValid()) {
        elinescolor = col;
        BTNelinescolor->setColor(&elinescolor);
    }
}

void Plotter::CHKhide_changed(){
    if (indcurve < 0) return;
    curvehide->replace(indcurve,!curvehide->at(indcurve));
}

void Plotter::CHKsghisthide_changed(){
    if (indcurve < 0) return;
    sghisthide->replace(indcurve,!sghisthide->at(indcurve));
    int i = indcurve+1;
    while ((i < sghisthide->length()) && (i < sghistlabelrect->length()) && (sghistlabelrect->at(i).size() == QSize(0,0))){
        sghisthide->replace(i++,sghisthide->at(indcurve));
    }
}

void Plotter::curvelabels_dialog(int i){
    indcurve = i;
    dialog = new QDialog();
    dialog->setWindowTitle(tr("Label"));
    dialog->setFixedHeight(60);
    dialog->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Fixed);
    dialog->setAttribute(Qt::WA_DeleteOnClose);

    texto = new QLineEdit();
    texto->setText(curvelabels->at(i));
    BTNcolor = new ColorButton();
    BTNcolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNcolor->setText(tr("Color"));
    color = curvecolors->operator [](i);
    BTNcolor->setColor(&color);
    BTNcolor->setEnabled(true);
    connect(BTNcolor, SIGNAL(clicked()), this, SLOT(BTNcolor_clicked()));
    CHKhide = new QCheckBox(tr("Hide"));   
    if (curvehide->at(indcurve)){
        CHKhide->setChecked(true);
    }
    else{
        CHKhide->setChecked(false);
    }
    connect(CHKhide, SIGNAL(stateChanged(int)), this, SLOT(CHKhide_changed()));
    QPushButton *BTNtitleOK = new QPushButton();
    BTNtitleOK->setText(tr("Accept"));
    connect(BTNtitleOK, SIGNAL(clicked()), this, SLOT(aceptar_label()));
    QPushButton *BTNtitleCancel = new QPushButton();
    BTNtitleCancel->setText(tr("Cancel"));
    connect(BTNtitleCancel, SIGNAL(clicked()), this, SLOT(cerrar_label()));
    QPushButton *BTNdeletecurve = new QPushButton();
    BTNdeletecurve->setText(tr("Delete"));
    connect(BTNdeletecurve, SIGNAL(clicked()), this, SLOT(emit_deletefrad()));
    connect(BTNdeletecurve, SIGNAL(clicked()), this, SLOT(cerrar_label()));
    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(texto);
    Layout1->addWidget(BTNcolor);
    Layout1->addWidget(CHKhide);
    Layout1->addWidget(BTNtitleOK);
    Layout1->addWidget(BTNtitleCancel);
    Layout1->addWidget(BTNdeletecurve);
    QVBoxLayout *Layout2 = new QVBoxLayout(dialog);
    Layout2->addLayout(Layout1);
    dialog->exec();
}

void Plotter::sghistlabels_dialog(int i){
    indcurve = i;
    dialog = new QDialog();
    dialog->setWindowTitle(tr("Label"));
    dialog->setFixedHeight(60);
    dialog->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Fixed);
    dialog->setAttribute(Qt::WA_DeleteOnClose);

    texto = new QLineEdit();
    texto->setText(sghistlabels->at(i));
    BTNcolor = new ColorButton();
    BTNcolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNcolor->setText(tr("Color"));
    color = sghistcolors->operator [](i);
    BTNcolor->setColor(&color);
    BTNcolor->setEnabled(true);
    connect(BTNcolor, SIGNAL(clicked()), this, SLOT(BTNcolor_clicked()));
    CHKhide = new QCheckBox(tr("Hide"));
    if (sghisthide->at(indcurve)){
        CHKhide->setChecked(true);
    }
    else{
        CHKhide->setChecked(false);
    }
    connect(CHKhide, SIGNAL(stateChanged(int)), this, SLOT(CHKsghisthide_changed()));
    QPushButton *BTNtitleOK = new QPushButton();
    BTNtitleOK->setText(tr("Accept"));
    connect(BTNtitleOK, SIGNAL(clicked()), this, SLOT(aceptar_sghistlabel()));
    QPushButton *BTNtitleCancel = new QPushButton();
    BTNtitleCancel->setText(tr("Cancel"));
    connect(BTNtitleCancel, SIGNAL(clicked()), this, SLOT(cerrar_label()));
    QPushButton *BTNdeletecurve = new QPushButton();
    BTNdeletecurve->setText(tr("Delete"));
    connect(BTNdeletecurve, SIGNAL(clicked()), this, SLOT(emit_deletesghist()));
    connect(BTNdeletecurve, SIGNAL(clicked()), this, SLOT(cerrar_label()));
    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(texto);
    Layout1->addWidget(BTNcolor);
    Layout1->addWidget(CHKhide);
    Layout1->addWidget(BTNtitleOK);
    Layout1->addWidget(BTNtitleCancel);
    Layout1->addWidget(BTNdeletecurve);
    QVBoxLayout *Layout2 = new QVBoxLayout(dialog);
    Layout2->addLayout(Layout1);
    dialog->exec();
}


void Plotter::aceptarcontour_label()
{
    if (indcurve >= 0)
        curvecolors->replace(indcurve,color);
    cerrar_label();
    refreshPixmap();
    update();
}

void Plotter::contourlabels_dialog(int i){
    indcurve = i;
    dialog = new QDialog();
    dialog->setWindowTitle(tr("Contour"));
    dialog->setFixedHeight(60);
    dialog->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Fixed);
    dialog->setAttribute(Qt::WA_DeleteOnClose);

    BTNcolor = new ColorButton();
    BTNcolor->setIcon(QIcon(":/images/colores48.png"));
    BTNcolor->setText(tr("Color"));
    color = curvecolors->operator [](i);
    BTNcolor->setColor(&color);
    BTNcolor->setEnabled(true);
    connect(BTNcolor, SIGNAL(clicked()), this, SLOT(BTNcolor_clicked()));
    QPushButton *BTNtitleOK = new QPushButton();
    BTNtitleOK->setText(tr("Accept"));
    connect(BTNtitleOK, SIGNAL(clicked()), this, SLOT(aceptarcontour_label()));
    QPushButton *BTNtitleCancel = new QPushButton();
    BTNtitleCancel->setText(tr("Cancel"));
    connect(BTNtitleCancel, SIGNAL(clicked()), this, SLOT(cerrar_label()));

    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(BTNcolor);
    Layout1->addWidget(BTNtitleOK);
    Layout1->addWidget(BTNtitleCancel);
    QVBoxLayout *Layout2 = new QVBoxLayout(dialog);
    Layout2->addLayout(Layout1);

    dialog->exec();
}

void Plotter::aceptar_elinescolor(){
    if (indeline < 0) return;
    QMapIterator<int, QVector<QPointF> > i(efieldMap);
    double dx;
    double dy;
    bool found = false;
    elinespenwidth = SPBelinespenwidth->value();
    if (RBTallelines->isChecked()){
        found = true;
        while (i.hasNext()) {
            i.next();
            efieldcolors->replace(i.key(), elinescolor);
            efieldpenwidth->replace(i.key(), elinespenwidth);
            plotarrows->replace(i.key(), CHKplotarrows->isChecked());
        }
    }
    else{
        while (i.hasNext()) {
            i.next();
            if (i.key() < indeline)
                continue;
            QVector<QPointF> data = i.value();
            dx = data[0].x();
            dy = data[0].y();
            found = true;
            efieldcolors->replace(i.key(), elinescolor);
            efieldpenwidth->replace(i.key(), elinespenwidth);
            plotarrows->replace(i.key(), CHKplotarrows->isChecked());
            break;

        }
        if (found && RBTatomelines->isChecked()){
            i.toFront();
            while (i.hasNext()) {
                i.next();
                QVector<QPointF> data = i.value();
                if (qAbs(data[0].x() - dx) < 1.e-3 && qAbs(data[0].y() - dy) < 1.e-3){
                    efieldcolors->replace(i.key(), elinescolor);
                    efieldpenwidth->replace(i.key(), elinespenwidth);
                    plotarrows->replace(i.key(), CHKplotarrows->isChecked());
                }
            }
        }
    }
    if (found){
        refreshPixmap();
        update();
    }  
    cerrar_elinescolor();
}

void Plotter::cerrar_elinescolor()
{
    if (texto != nullpointer){
        delete(texto);
        texto = nullpointer;
    }
    if (BTNelinescolor != nullpointer){
        delete(BTNelinescolor);
        BTNelinescolor = nullpointer;
    }
    if (RBTallelines != nullpointer){
        delete(RBTallelines);
        RBTallelines = nullpointer;
    }
    if (RBTatomelines != nullpointer){
        delete(RBTatomelines);
        RBTatomelines = nullpointer;
    }
    if (RBTsingleeline != nullpointer){
        delete(RBTsingleeline);
        RBTsingleeline = nullpointer;
    }
    if (LBLefpenwidth != nullpointer){
        delete(LBLefpenwidth);
        LBLefpenwidth = nullpointer;
    }
    if (SPBelinespenwidth != nullpointer){
        delete(SPBelinespenwidth);
        SPBelinespenwidth = nullpointer;
    }
    if (CHKplotarrows != nullpointer){
        delete(CHKplotarrows);
        CHKplotarrows = nullpointer;
    }
    if (FRMelinescolor != nullpointer){
        delete(FRMelinescolor);
        FRMelinescolor = nullpointer;
    }    
    if (dialog != nullpointer){
        delete(dialog);
        dialog = nullpointer;
    }
}

void Plotter::elines_dialog(int i){
    indeline = i;
    dialog = new QDialog();
    dialog->setWindowTitle(tr("Lines"));
    dialog->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
    dialog->setAttribute(Qt::WA_DeleteOnClose);

    FRMelinescolor = new QGroupBox();

    BTNelinescolor = new ColorButton();
    BTNelinescolor->setIcon(QIcon(":/images/colores48.png"));
    BTNelinescolor->setText(tr("Color"));
    elinescolor = efieldcolors->operator [](i);
    BTNelinescolor->setColor(&elinescolor);
    BTNelinescolor->setEnabled(true);
    connect(BTNelinescolor, SIGNAL(clicked()), this, SLOT(BTNelinescolor_clicked()));
    RBTsingleeline = new QRadioButton(tr("Single line"), FRMelinescolor);
    RBTsingleeline->setChecked(false);
    RBTatomelines = new QRadioButton(tr("Collective lines"), FRMelinescolor);
    RBTatomelines->setChecked(true);
    RBTallelines = new QRadioButton(tr("All lines"), FRMelinescolor);
    RBTallelines->setChecked(false);

    LBLefpenwidth = new QLabel(tr("Width"));
    SPBelinespenwidth = new QSpinBox();
    SPBelinespenwidth->setRange(0,8);
    SPBelinespenwidth->setValue(elinespenwidth);

    CHKplotarrows = new QCheckBox(tr("Show arrows"));
    CHKplotarrows->setChecked(plotarrows->at(i));

    QPushButton *BTNtitleOK = new QPushButton();
    BTNtitleOK->setText(tr("Accept"));
    connect(BTNtitleOK, SIGNAL(clicked()), this, SLOT(aceptar_elinescolor()));
    QPushButton *BTNtitleCancel = new QPushButton();
    BTNtitleCancel->setText(tr("Cancel"));
    connect(BTNtitleCancel, SIGNAL(clicked()), this, SLOT(cerrar_elinescolor()));

    QVBoxLayout *Layout1 = new QVBoxLayout(FRMelinescolor);
    Layout1->addWidget(RBTsingleeline);
    Layout1->addWidget(RBTatomelines);
    Layout1->addWidget(RBTallelines);

    QHBoxLayout *Layout2 = new QHBoxLayout();
    Layout2->addWidget(LBLefpenwidth);
    Layout2->addWidget(SPBelinespenwidth);

    QVBoxLayout *Layout3 = new QVBoxLayout();
    Layout3->addWidget(BTNelinescolor);
    Layout3->addLayout(Layout2);
    Layout3->addWidget(CHKplotarrows);
    Layout3->addWidget(BTNtitleOK);
    Layout3->addWidget(BTNtitleCancel);

    QHBoxLayout *Layout = new QHBoxLayout(dialog);
    Layout->addStretch();
    Layout->addWidget(FRMelinescolor);
    Layout->addStretch();
    Layout->addLayout(Layout3);
    Layout->addStretch();

    dialog->exec();
}

void Plotter::emit_deletefrad(){
    emit deletefrad(indcurve);
}

void Plotter::emit_deletesghist(){
    emit deletesghist(indcurve);
}

bool Plotter::getaspectratio(){
    return aspectratio;
}
bool Plotter::getautomaticticks(){
    return automaticticks;
}
bool Plotter::getcontourplot(){
    return contourplot;
}
bool Plotter::getfradplot(){
    return fradplot;
}
bool Plotter::getfieldplot(){
    return fieldplot;
}
bool Plotter::getmulticolor(){
    return aspectratio;
}
bool Plotter::getsghistplot(){
    return sghistplot;
}
bool Plotter::getshowcoords(){
    return showcoords;
}
bool Plotter::getshowgrid(){
    return showgrid;
}
bool Plotter::getshowcontourlabels(){
    return showcontourlabels;
}
bool Plotter::getshowcurvelabels(){
    return showcurvelabels;
}
bool Plotter::getshowsghistlabels(){
    return showsghistlabels;
}
bool Plotter::getshowXscalebottom(){
    return showXscalebottom;
}
bool Plotter::getshowXscaletop(){
    return showXscaletop;
}
bool Plotter::getshowYscaleleft(){
    return showYscaleleft;
}
bool Plotter::getshowYscaleright(){
    return showYscaleright;
}
bool Plotter::getshowXticksbottom(){
    return showXticksbottom;
}
bool Plotter::getshowXtickstop(){
    return showXtickstop;
}
bool Plotter::getshowYticksleft(){
    return showYticksleft;
}
bool Plotter::getshowYticksright(){
    return showYticksright;
}
bool Plotter::getshowtitle(){
    return showtitle;
}
bool Plotter::getshowtitlebox(){
    return showtitlebox;
}
bool Plotter::getshowXlabel(){
    return showXlabel;
}
bool Plotter::getshowXlabelbox(){
    return showXlabelbox;
}
bool Plotter::getshowYlabel(){
    return showYlabel;
}
bool Plotter::getshowYlabelbox(){
    return showYlabelbox;
}
bool Plotter::getsolidgrid(){
    return solidgrid;
}
bool Plotter::gettranspbackground(){
    return transpbackground;
}
bool Plotter::gettranspbckgrcapture(){
    return transpbckgrcapture;
}
bool Plotter::getYlabelvert(){
    return Ylabelvert;
}
double Plotter::getcontourtolerance(){
    return contourtolerance;
};
double Plotter::getumin(){
    return umin;
}
double Plotter::getumax(){
    return umax;
}
double Plotter::getvmin(){
    return vmin;
}
double Plotter::getvmax(){
    return vmax;
}
double Plotter::getxmin(){
    return xmin;
}
double Plotter::getxmax(){
    return xmax;
}
double Plotter::getymin(){
    return ymin;
}
double Plotter::getymax(){
    return ymax;
}
int Plotter::getimagequality(){
    return imagequality;
}
int Plotter::getnlevels(){
    return nlevels;
}
int Plotter::getnu(){
    return nu;
}
int Plotter::getnv(){
    return nv;
}
int Plotter::getXcifras(){
    return Xcifras;
}
int Plotter::getYcifras(){
    return Ycifras;
}
int Plotter::getnegativecontourstyle(){
    return negativecontourstyle;
}
int Plotter::getpositivecontourstyle(){
    return positivecontourstyle;
}
int Plotter::getzerocontourstyle(){
    return zerocontourstyle;
}
int Plotter::getbasinspenwidth(){
    return basinspenwidth;
}
int Plotter::getcntpenwidth(){
    return cntpenwidth;
}
int Plotter::getcurvespenwidth(){
    return curvespenwidth;
}
int Plotter::getTopMargin(){
    return TopMargin;
}
int Plotter::getBottomMargin(){
    return BottomMargin;
}
int Plotter::getLeftMargin(){
    return LeftMargin;
}
int Plotter::getRightMargin(){
    return RightMargin;
}
int Plotter::getzerotolerancepow(){
    return zerotolerancepow;
}
QColor Plotter::getbackgroundcolor(){
    return backgroundcolor;
}
QColor Plotter::getbasinscolor(){
    return basinscolor;
}
QColor Plotter::getcenterscolor(){
    return centerscolor;
}
QColor Plotter::getcenterslabelcolor(){
    return centerslabelcolor;
}
QColor Plotter::getcontourslabelcolor(){
    return contourslabelcolor;
}
QColor Plotter::getfonttitlecolor(){
    return fonttitlecolor;
}
QColor Plotter::getfontXlabelcolor(){
    return fontXlabelcolor;
}
QColor Plotter::getfontYlabelcolor(){
    return fontYlabelcolor;
}
QColor Plotter::getgridcolor(){
    return gridcolor;
}
QColor Plotter::getnegativecontourcolor(){
    return negativecontourcolor;
}
QColor Plotter::getpositivecontourcolor(){
    return positivecontourcolor;
}
QColor Plotter::getzerocontourcolor(){
    return zerocontourcolor;
}
QFont  Plotter::getfontcenterslabel(){
    return fontcenterslabel;
}
QFont  Plotter::getfontcurvelabels(){
    return fontcurvelabels;
}
QFont  Plotter::getfontcontour(){
    return fontcontour;
}
QFont  Plotter::getfontsghistlabels(){
    return fontsghistlabels;
}
QFont  Plotter::getfonttitle(){
    return fonttitle;
}
QFont  Plotter::getfontXlabel(){
    return fontXlabel;
}
QFont  Plotter::getfontYlabel(){
    return fontYlabel;
}

void Plotter::setimagesize(QSize a){
    imagesize.setHeight(a.height());
    imagesize.setWidth(a.width());
}
void Plotter::setaspectratio(bool a){
    aspectratio = a;
}
void Plotter::setautomaticticks(bool a){
    automaticticks = a;
}
void Plotter::setcontourplot(bool a){
    contourplot = a;
}
void Plotter::setfradplot(bool a){
    fradplot = a;
}
void Plotter::setdrawarrows(bool a){
    drawarrows = a;
}
void Plotter::setexistcntatlbls(bool a){
    existcntatlbls = a;
}

void Plotter::setfieldplot(bool a){
    fieldplot = a;
}
void Plotter::setimagequality(int a){
    imagequality = a;
}
void Plotter::setinvertarrows(bool a){
    invertarrows = a;
}
void Plotter::setmulticolor(bool a){
    multicolor = a;   
}
void Plotter::setplotbasins(bool a){
    plotbasins = a;
}
void Plotter::setplotbonds(bool a){
    plotbonds = a;
}
void Plotter::setplotcenters(bool a){
    plotcenters = a;
}
void Plotter::setplotcps(bool a){
    plotcps = a;
}
void Plotter::setsghistplot(bool a){
    sghistplot = a;
}
void Plotter::setshowcenterslabel(bool a){
    showcenterslabel = a;
}
void Plotter::setshowcoords(bool a){
    showcoords = a;
}
void Plotter::setshowgrid(bool a){
    showgrid = a;
}
void Plotter::setshowcontourlabels(bool a){
    showcontourlabels = a;
}
void Plotter::setshowcurvelabels(bool a){
    showcurvelabels = a;
}
void Plotter::setshowsghistlabels(bool a){
    showsghistlabels = a;
}
void Plotter::setshowXscalebottom(bool a){
    showXscalebottom = a;
}
void Plotter::setshowXscaletop(bool a){
    showXscaletop = a;
}
void Plotter::setshowYscaleleft(bool a){
    showYscaleleft = a;
}
void Plotter::setshowYscaleright(bool a){
    showYscaleright = a;
}
void Plotter::setshowXticksbottom(bool a){
    showXticksbottom = a;
}
void Plotter::setshowXtickstop(bool a){
    showXtickstop = a;
}
void Plotter::setshowYticksleft(bool a){
    showYticksleft = a;
}
void Plotter::setshowYticksright(bool a){
    showYticksright = a;
}
void Plotter::setshowtitle(bool a){
    showtitle = a;
}
void Plotter::setshowtitlebox(bool a){
    showtitlebox = a;
}
void Plotter::setshowXlabel(bool a){
    showXlabel = a;
}
void Plotter::setshowXlabelbox(bool a){
    showXlabelbox = a;
}
void Plotter::setshowYlabel(bool a){
    showYlabel = a;
}
void Plotter::setshowYlabelbox(bool a){
    showYlabelbox = a;
}
void Plotter::setsolidgrid(bool a){
    solidgrid = a;
}
void Plotter::settranspbackground(bool a){
    transpbackground = a;
}
void Plotter::settranspbckgrcapture(bool a){
    transpbckgrcapture = a;
}
void Plotter::setuserdefresol(bool a){
    userdefresol = a;
}
void Plotter::setYlabelvert(bool a){
    Ylabelvert = a;
}
void Plotter::setbondsthr(double a){
    bondsthr = a;
};
void Plotter::setcontourtolerance(double a){
    contourtolerance = a;
};
void Plotter::setumin(double a){
    if (a < umax)
        umin = a;
    else if (a > umax){
        umin = umax;
        umax = a;
    }
    else{
        umax = a + 1.;
    }
}
void Plotter::setumax(double a){
    if (a > umin)
        umax = a;
    else if(a < umin){
        umax = umin;
        umin = a;
    }
    else{
        umax = a + 1.;
    }
}
void Plotter::setvmin(double a){
    if (a < vmax)
        vmin = a;
    else if (a > vmax){
        vmin = vmax;
        vmax = a;
    }
    else{
        vmax = a + 1.;
    }
}
void Plotter::setvmax(double a){
    if (a > vmin)
        vmax = a;
    else if(a < vmin){
        vmax = vmin;
        vmin = a;
    }
    else{
        vmax = a + 1.;
    }
}
void Plotter::setumincnt(double a){
    if (a < umaxcnt)
        umincnt = a;
    else if (a > umaxcnt){
        umincnt = umax;
        umaxcnt = a;
    }
    else{
        umaxcnt = a + 1.;
    }
}
void Plotter::setumaxcnt(double a){
    if (a > umincnt)
        umaxcnt = a;
    else if(a < umincnt){
        umaxcnt = umincnt;
        umincnt = a;
    }
    else{
        umaxcnt = a + 1.;
    }
}
void Plotter::setvmincnt(double a){
    if (a < vmaxcnt)
        vmincnt = a;
    else if (a > vmaxcnt){
        vmincnt = vmaxcnt;
        vmaxcnt = a;
    }
    else{
        vmaxcnt = a + 1.;
    }
}
void Plotter::setvmaxcnt(double a){
    if (a > vmincnt)
        vmaxcnt = a;
    else if(a < vmincnt){
        vmaxcnt = vmincnt;
        vmincnt = a;
    }
    else{
        vmaxcnt = a + 1.;
    }
}
void Plotter::setuminelines(double a){
    if (a < umaxelines)
        uminelines = a;
    else if (a > umaxelines){
        uminelines = umax;
        umaxelines = a;
    }
    else{
        umaxelines = a + 1.;
    }
}
void Plotter::setumaxelines(double a){
    if (a > uminelines)
        umaxelines = a;
    else if(a < uminelines){
        umaxelines = uminelines;
        uminelines = a;
    }
    else{
        umaxelines = a + 1.;
    }
}
void Plotter::setvminelines(double a){
    if (a < vmaxelines)
        vminelines = a;
    else if (a > vmaxelines){
        vminelines = vmaxelines;
        vmaxelines = a;
    }
    else{
        vmaxelines = a + 1.;
    }
}
void Plotter::setvmaxelines(double a){
    if (a > vminelines)
        vmaxelines = a;
    else if(a < vminelines){
        vmaxelines = vminelines;
        vminelines = a;
    }
    else{
        vmaxelines = a + 1.;
    }
}
void Plotter::setxmin(double a){
    if (a < xmax)
        xmin = a;
    else if (a > xmax){
        xmin = xmax;
        xmax = a;
    }
    else{
        xmax = a + 1.;
    }
}
void Plotter::setxmax(double a){
    if (a >xmin)
        xmax = a;
    else if(a < xmin){
        xmax = xmin;
        xmin = a;
    }
    else{
        xmax = a + 1.;
    }
}
void Plotter::setymin(double a){
    if (a < ymax)
        ymin = a;
    else if (a > ymax){
        ymin = ymax;
        ymax = a;
    }
    else{
        ymax = a + 1.;
    }
}
void Plotter::setymax(double a){
    if (a > ymin)
        ymax = a;
    else if(a < ymin){
        ymax = ymin;
        ymin = a;
    }
    else{
        ymax = a + 1.;
    }
}

void Plotter::setarrowsseparation(int a){
    arrowsseparation = a;
}
void Plotter::setarrowssize(int a){
    arrowssize = a;
}
void Plotter::setarrowsskew(int a){
    arrowsskew = a;
}
void Plotter::setarrowswidth(int a){
    arrowswidth = a;
}
void Plotter::setshowcps(int i, bool a){
    showcps[i] = a;
}
void Plotter::setbackgroundcolor(QColor a){
    backgroundcolor = a;
}
void Plotter::setbasinscolor(QColor a){
    basinscolor = a;
}
void Plotter::setbonds(){
    Elements *Elem = new Elements();
    bonds->clear();
    for (int i = 1 ; i < rcen->count() ; i++){
        QVector3D ri = rcen->at(i);
        for (int j = 0 ; j < i ; j++){
            QVector3D rj = rcen->at(j);
            double dist = (ri - rj).length();
            double threshold = bondsthr * (Elem->getCovalentRadius((int)(uvznuc->at(i).z()+0.5))
                        + Elem->getCovalentRadius((int)(uvznuc->at(j).z()+0.5)))  * ANGSTROMTOBOHR;
            if (dist < threshold){
                bonds->append(QVector2D(uvznuc->at(i).x(),uvznuc->at(i).y()));
                bonds->append(QVector2D(uvznuc->at(j).x(),uvznuc->at(j).y()));
            }
        }
    }
}
void Plotter::setbondscolor(QColor a){
    bondscolor = a;
}
void Plotter::setbondswidth(int a){
    bondswidth = a;
}
void Plotter::setcenterscolor(QColor a){
    centerscolor = a;
}
void Plotter::setcenterslabel(){
    Elements *Elem = new Elements();
    centerslabel->clear();
    for (int i = 0 ; i < uvznuc->count() ; i++){
        centerslabel->append(Elem->getSymbol((int)(uvznuc->at(i).z()+0.5)));
    }
}
void Plotter::setcenterslabelcolor(QColor a){
    centerslabelcolor = a;
}
void Plotter::setcenterslabelshft(int a){
    centerslabelshft = a;
}
void Plotter::setcenterradius(int a){
    centerradius = a;
}
void Plotter::setcpsradius(int a){
    cpsradius = a;
}
void Plotter::setcpstype(QVector<int> a){
    cpstype.clear();
    cpstype = a;
}
void Plotter::setcontourslabelcolor(QColor a){
    contourslabelcolor = a;
}


void Plotter::setnlevels(int a){
    nlevels = a;
}
void Plotter::setnu(int a){
    nu = a;
}
void Plotter::setnv(int a){
    nv = a;
}
void Plotter::setplaneABC(QVector3D a){
    planeABC = a;
}
void Plotter::setXcifras(int a){
    Xcifras = a;
}
void Plotter::setYcifras(int a){
    Ycifras = a;
}
void Plotter::setnegativecontourstyle(int a){
    negativecontourstyle = a;
}
void Plotter::setpositivecontourstyle(int a){
    positivecontourstyle = a;
}
void Plotter::setzerocontourstyle(int a){
    zerocontourstyle = a;
}
void Plotter::setbasinspenwidth(int a){
    basinspenwidth = a;
}
void Plotter::setcntpenwidth(int a){
    cntpenwidth = a;
}
void Plotter::setcurvespenwidth(int a){
    curvespenwidth = a;
}
void Plotter::setsghistpenwidth(int a){
    sghistpenwidth = a;
}
void Plotter::setTopMargin(int a){
    TopMargin = a;
    zoomStack[curZoom].TopMargin = a;
}
void Plotter::setBottomMargin(int a){
    BottomMargin = a;
    zoomStack[curZoom].BottomMargin = a;
}
void Plotter::setLeftMargin(int a){
    LeftMargin = a;
    zoomStack[curZoom].LeftMargin = a;
}
void Plotter::setRightMargin(int a){
    RightMargin = a;
    zoomStack[curZoom].RightMargin = a;
}
void Plotter::setwhoisabove(int a){
    if (a == 2 || a == 3)
        whoisabove = a;
    else
        whoisabove = 2;
}
void Plotter::setzerotolerancepow(int a){
    zerotolerancepow = a;
}
void Plotter::setzero(double a){
    zero = a;
}

void Plotter::setcpscolor(QVector<QColor> a){
    cpscolor = a;
}
void Plotter::setbasinsuv( QVector<QVector <QVector2D> > a){
    basinsuv.clear();
    basinsuv = a;
}
void Plotter::setcpsuvval(QVector<QVector3D> a){
    cpsuvval.clear();
    cpsuvval = a;
}
void Plotter::setfonttitlecolor(QColor a){
    fonttitlecolor = a;
}
void Plotter::setfontXlabelcolor(QColor a){
    fontXlabelcolor = a;
}
void Plotter::setfontYlabelcolor(QColor a){
    fontYlabelcolor = a;
}
void Plotter::setgridcolor(QColor a){
    gridcolor = a;
}
void Plotter::setnegativecontourcolor(QColor a){
    negativecontourcolor = a;
}
void Plotter::setpositivecontourcolor(QColor a){
    positivecontourcolor = a;
}
void Plotter::setzerocontourcolor(QColor a){
    zerocontourcolor = a;
}
void Plotter::setfontcenterslabel(QFont a){
    fontcenterslabel = a;
}
void Plotter::setfontcurvelabels(QFont a){
    fontcurvelabels = a;
}
void Plotter::setfontcontour(QFont a){
    fontcontour = a;
    for (int i = 0 ; i < contourlabels.length() ; i++){
        ContourLabel clabel = contourlabels[i];
        QFontMetrics fm = QFontMetrics( fontcontour );
        QSize fmsize = fm.size( Qt::TextSingleLine, clabel.getlabel() );
        clabel.setrect(QRect(clabel.getx() - fmsize.width()/2, clabel.gety() - fmsize.height()/2 ,
                             fmsize.width(), fmsize.height())) ;
        contourlabels.replace(i,clabel);
    }
}
void Plotter::setfontsghistlabels(QFont a){
    fontsghistlabels = a;
}
void Plotter::setfonttitle(QFont a){
    fonttitle = a;
}
void Plotter::setfontXlabel(QFont a){
    fontXlabel = a;
}
void Plotter::setfontYlabel(QFont a){
    fontYlabel = a;
}


//***************************************************************************
//********************  PLOTSETTINGS FUNCTIONS    ***************************
//***************************************************************************

PlotSettings::PlotSettings()
{
    minX = 0.0;
    maxX = 10.0;
    numXTicks = 5;

    minY = 0.0;
    maxY = 10.0;
    numYTicks = 5;
    
    LeftMargin = 80;
    RightMargin = 80;
    TopMargin = 80;
    BottomMargin = 80;
}

void PlotSettings::scroll(int dx, int dy)
{
    double stepX = spanX() / numXTicks;
    minX += dx * stepX;
    maxX += dx * stepX;

    double stepY = spanY() / numYTicks;
    minY += dy * stepY;
    maxY += dy * stepY;
}

void PlotSettings::set_numXticks(int numx)
{
    numXTicks = numx-1;
}

void PlotSettings::set_numYticks(int numy)
{
    numYTicks = numy-1;
}

void PlotSettings::adjust()
{
    adjustAxis(minX, maxX, numXTicks);
    adjustAxis(minY, maxY, numYTicks);
}

void PlotSettings::adjustAxis(double &min, double &max, int &numTicks)
{
    const int MinTicks = 4;
    double grossStep = (max - min) / MinTicks;
    double step = std::pow(10.0, std::floor(std::log10(grossStep)));

    if (5 * step < grossStep) {
        step *= 5;
    } else if (2 * step < grossStep) {
        step *= 2;
    }

    numTicks = int(std::ceil(max / step) - std::floor(min / step));
    if (numTicks < MinTicks)
        numTicks = MinTicks;
    min = std::floor(min / step) * step;
    max = std::ceil(max / step) * step;
}

void PlotSettings::adjustforced()
{
    adjustAxisforced(minX, maxX, numXTicks);
    adjustAxisforced(minY, maxY, numYTicks);
}

void PlotSettings::adjustAxisforced(double &min, double &max, int &numTicks)
{
    double grossStep = (max - min) / numTicks;
    double step = std::pow(10.0, std::floor(std::log10(grossStep)));
    min = std::floor(min / step) * step;
    max = std::ceil(max / step) * step;
}


//
//    MainWindow2DViewer SLOTS
//

ContourLabel::ContourLabel()
{
    u = 0.;
    v = 0.;
    x = 0;
    y = 0;
    label=QString("");
    rect = QRect(0,0,0,0);
}

ContourLabel::~ContourLabel(){
}

double ContourLabel::getu(){
    return u;
}

double ContourLabel::getv(){
    return v;
}

int ContourLabel::getx(){
    return x;
}

int ContourLabel::gety(){
    return y;
}

QRect ContourLabel::getrect(){
    return rect;
}

QString ContourLabel::getlabel(){
    return label;
}

void ContourLabel::setlabel(QString str){
    label = str;
}

void ContourLabel::setrect(QRect r){
    rect = r;
}

void ContourLabel::setu(double i){
    u = i;
}

void ContourLabel::setv(double i){
    v = i;
}

void ContourLabel::setx(int i){
    x = i;
}

void ContourLabel::sety(int i){
    y = i;
}


