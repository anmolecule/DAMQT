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
//  File:   Sheet.h
//
//      Last version: March 2016
//
#ifndef SHEET_H
#define SHEET_H

#include <qwidget.h>
#include <QItemDelegate>
#include <QModelIndex>
#include <QObject>
#include <QSize>
#include <QLineEdit>
#include <QSpinBox>
#include <QRect>

class QTableWidget;
class QToolButton;
class QGroupBox;
class QDoubleValidator;
class QTableWidgetItem;
class LineEditDelegate;

class Sheet : public QWidget
{
    Q_OBJECT
public:
    Sheet(int rows,int cols,int max, bool type, QWidget *parent = 0, int width = 56);
	void setVisible(bool visible);
	void setHeader(QStringList qsl);
    void clear();
	QString getcellvalue(int row,int col);
	void setcellvalue(QString value,int row,int col);
    void setcellfloatvalue(QString value,int row,int col);
	void resizeRows(int rows);
	QTableWidget *tabla;
	void max_changed(int max);
	void min_changed(int min);
	static const int max_sel = 200; // Maximum number of atoms for atomic tabulations

public slots:
    void BTNadd_clicked();
	void BTNdel_clicked();

signals:
    void modified();

private slots:
    void somethingChanged();

private:
	LineEditDelegate *del;
	QToolButton *BTNadd;
	QToolButton *BTNdel;
	
};

class LineEditDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    LineEditDelegate(QObject *parent = 0);

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const;

    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const;

    void updateEditorGeometry(QWidget *editor,
        const QStyleOptionViewItem &option, const QModelIndex &index) const;

	bool type; // if false: integer type entries (spinbox), if true: double type entries
	int max; // highest allowed value in spinbox (integer entries)
	int min; // lowest allowed value in spinbox (integer entries)
};

#endif

