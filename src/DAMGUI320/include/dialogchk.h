#ifndef DIALOGCHK_H
#define DIALOGCHK_H

#include <QObject>
#include <QWidget>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>

class dialogCHK : public QObject
{
    Q_OBJECT
public:
    explicit dialogCHK(QObject *parent = 0);

signals:

public slots:
    void create();
    void delet();

private:
    QGroupBox *FRMdialogCHK;
    QLabel *LBLdialogCHK;
    QPushButton *BTNcreate;
    QPushButton *BTNdelete;
};

#endif // DIALOGCHK_H
