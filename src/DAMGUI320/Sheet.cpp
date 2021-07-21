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
//	Class for displaying a Sheet admitting data
//  
//  File:   Sheet.cpp
//
//      Last version: September 2018
//
#include <QApplication>
#include <QHeaderView>
#include <QLayout>
#include <QScrollBar>
#include <QToolButton>
#include <QTableWidget>
#include "Sheet.h"

#include <QtDebug>


Sheet::Sheet(int rows,int cols,int max,bool type,QWidget *parent): QWidget(parent)
{
//	lang = new IDIOMA();
	QString path=QApplication::applicationDirPath();
	tabla=new QTableWidget(rows,cols,parent);
	tabla->setSelectionBehavior(QAbstractItemView::SelectRows);
	tabla->setSelectionMode(QAbstractItemView::SingleSelection); //NoSelection
	tabla->setAlternatingRowColors(true);
	for(int i=0; i < cols ; i++){
		tabla->setColumnWidth(i,56);
	}
	del=new LineEditDelegate;
	del->max=max;
	del->type=type;
	tabla->setItemDelegate(del);

	BTNadd=new QToolButton(parent);
	BTNadd->setIcon(QIcon(":/images/add.png"));
	BTNadd->setToolTip(tr("Adds a new entry to table"));
	connect(BTNadd, SIGNAL(clicked()), this, SLOT(BTNadd_clicked()));
	BTNdel=new QToolButton(parent);
	BTNdel->setIcon(QIcon(":/images/delete.png"));
	BTNdel->setToolTip(tr("Removed selected entries from table"));
	connect(BTNdel, SIGNAL(clicked()), this, SLOT(BTNdel_clicked()));

	QHBoxLayout *Layout2=new QHBoxLayout();
	Layout2->addWidget(BTNadd);
	Layout2->addWidget(BTNdel);
	Layout2->setAlignment(Qt::AlignRight);
	QHBoxLayout *Layout3=new QHBoxLayout();
	Layout3->addWidget(tabla);
	Layout3->setAlignment(Qt::AlignCenter);
	QVBoxLayout *Layout1=new QVBoxLayout(parent);
	Layout1->addLayout(Layout3);
	Layout1->addLayout(Layout2);
	Layout1->addStretch();
	resizeRows(tabla->rowCount());
}

void Sheet::max_changed(int max)
{
	del->max=max;
}

void Sheet::min_changed(int min)
{
	del->min=min;
}

void Sheet::setHeader(QStringList qsl)
{
	tabla->setHorizontalHeaderLabels(qsl);
}

void Sheet::setVisible(bool visible)
{
	tabla->setVisible(visible);
	BTNadd->setVisible(visible);
	BTNdel->setVisible(visible);
}

void Sheet::somethingChanged()
{
    emit modified();
}


void Sheet::BTNadd_clicked()
{
	int rows=tabla->rowCount();
	if (tabla->rowCount() < max_sel) {
		resizeRows(rows+1);
		tabla->insertRow(rows);
		resizeRows(rows+1);
		int col=tabla->columnCount();
		for (int i=0;i<col;i++)
			setcellvalue("",rows,i);
	}
}

void Sheet::BTNdel_clicked()
{
	int cr=tabla->currentRow();
	if (cr>=0){
		int rows=tabla->rowCount();
		if(rows > 0)
		{
			resizeRows(rows-1);
			tabla->removeRow(tabla->currentRow());
			resizeRows(rows-1);	
		}
	}
}

void Sheet::resizeRows(int rows)
{
	tabla->resizeRowsToContents();
	int rowheight=tabla->rowHeight(0);
	int height=0;
	QScrollBar *Hscrollbar=tabla->horizontalScrollBar();
	QHeaderView *Hheader=tabla->horizontalHeader();
#if QT_VERSION >= 0x050000
        Hheader->setSectionsClickable(false);
        if (Hheader->count()<4){
		Hheader->setSectionResizeMode(QHeaderView::Stretch);
	}else{
		height=Hscrollbar->height();
	}
#else
        Hheader->setClickable(false);
        if (Hheader->count()<4){
		Hheader->setResizeMode(QHeaderView::Stretch);
	}else{
		height=Hscrollbar->height();
	}
#endif
	if (rows==0){
		height=height+2*Hheader->height();
	}else if (rows>0 && rows<5){
		height=height+rows*rowheight+Hheader->height();
	}else{
		height=height+5*rowheight+Hheader->height();
	}
	tabla->setMinimumHeight(height+5);
	tabla->setMaximumHeight(height+5);
}

QString Sheet::getcellvalue(int row,int col)
{
	QTableWidgetItem *newItem = tabla->item(row,col);
	if (newItem != 0){
		QString qv = newItem->text();
		return qv;
	}else{
		return 0;
	}
}

void Sheet::setcellvalue(QString value,int row,int col)
{
	QTableWidgetItem *newItem = new QTableWidgetItem(row*col);
	if (value.toDouble() == 0)
        if (value.size() > 0)
            newItem->setText("0");
        else
            newItem->setText("");
	else
		newItem->setText(value);
	newItem->setTextAlignment(Qt::AlignCenter);
	tabla->setItem(row,col,newItem);
}

void Sheet::setcellfloatvalue(QString value,int row,int col)
{
	QTableWidgetItem *newItem = new QTableWidgetItem(row*col);
        newItem->setText(value);
	newItem->setTextAlignment(Qt::AlignCenter);
	tabla->setItem(row,col,newItem);
}

void Sheet::clear()
{
	tabla->setRowCount(0);
}

//*************************************************************************
//*****************************  DELEGATE  ********************************
//*************************************************************************
LineEditDelegate::LineEditDelegate(QObject *parent): QItemDelegate(parent)
{
}


QWidget *LineEditDelegate::createEditor(QWidget *parent,
    const QStyleOptionViewItem &/* option */,
    const QModelIndex &index) const
{
	QDoubleValidator *myDoubleValidator = new QDoubleValidator(0);
	myDoubleValidator->setLocale(QLocale::English);
	if (type){
		QLineEdit *editor = new QLineEdit(parent);
		editor->setValidator(myDoubleValidator);
		editor->installEventFilter(const_cast<LineEditDelegate*>(this));
		return editor;
	}else{
		QSpinBox *editor = new QSpinBox(parent);
		editor->setMinimum(min);
		editor->setMaximum(max);
		editor->installEventFilter(const_cast<LineEditDelegate*>(this));
		return editor;
	}	
}

void LineEditDelegate::setEditorData(QWidget *editor,const QModelIndex &index) const
{
    if (!type){
		int value = index.model()->data(index, Qt::DisplayRole).toInt();

		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->setValue(value);
		spinBox->setCorrectionMode(QAbstractSpinBox::CorrectToNearestValue); // To fill with zero when content is deleted 
	}else{
		QString value=index.model()->data(index,Qt::DisplayRole).toString();

		QLineEdit *le = static_cast<QLineEdit*>(editor);
		le->setText(value);
	}
}

void LineEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                   const QModelIndex &index) const
{
	if (!type){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		int value = spinBox->value();
		if (value != 0)
			model->setData(index, value);
		else
			model->setData(index,"");
	}else{
		QLineEdit *le = static_cast<QLineEdit*>(editor);

	    model->setData(index, le->text());
	}
}

void LineEditDelegate::updateEditorGeometry(QWidget *editor,
    const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
    editor->setGeometry(option.rect);
}
