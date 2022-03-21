//  Copyright 2008-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//  File:   main.cpp
//
//      Last version: September 2018
//
#include <QtGlobal>

#include <QApplication>
#include <QStyleFactory>
#include <QSplashScreen>
#include <QTranslator>
#include <QLibraryInfo>

#include <QtDebug>

#include "mainwindow.h"

#include <stdio.h>
#include <stdlib.h>

#if QT_VERSION >= 0x050000
    void myMessageOutput(QtMsgType type, const QMessageLogContext& context, const QString &msg)
    {
        //in this function, you can write the message to any stream!
       QByteArray localMsg = msg.toLocal8Bit();
       fprintf(stderr, "MESSAGE (%s:%u %s): %s\n", context.file, context.line, context.function, localMsg.constData());
       fflush(stderr);

    }
#else
    void myMessageOutput(QtMsgType type, const char *msg)
    {
        switch (type) {
        case QtDebugMsg:
            fprintf(stderr, "Debug: %s\n", msg);
            break;
        case QtWarningMsg:
            fprintf(stderr, "Warning: %s\n", msg);
            break;
        case QtCriticalMsg:
            fprintf(stderr, "Critical: %s\n", msg);
            break;
        case QtFatalMsg:
            fprintf(stderr, "Fatal: %s\n", msg);
            abort();
        }
    }
#endif

int main(int argc, char *argv[])
{
    Q_INIT_RESOURCE(damqt);
#if QT_VERSION >= 0x050000
        qInstallMessageHandler(myMessageOutput); //install : set the callback
#else
        qInstallMsgHandler(myMessageOutput);
#endif        
    QApplication app(argc, argv);
#if QT_VERSION >= 0x050000
    app.setStyle(QStyleFactory::create("plastique"));
#else 
    qDebug() << "QT_VERSION < 0x050000" ;
    QApplication::setStyle(QStyleFactory::create("plastique")); //windowsxp,macintosh
#endif
    QTranslator qtTranslator;
    qtTranslator.load("qt_" + QLocale::system().name(),
    QLibraryInfo::location(QLibraryInfo::TranslationsPath));
    app.installTranslator(&qtTranslator);
    app.processEvents();
    MainWindow mainWin;
    mainWin.show();
    mainWin.finishsplash();
    return app.exec();
}
