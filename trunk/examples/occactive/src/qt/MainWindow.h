
#pragma once

#include <QtGui/QApplication>
#include <QtGui/QProgressBar>
#include <QtGui/QStandardItemModel>
#include <QtGui/QStandardItem>
#include <QtGui/QDockWidget>
#include <QtCore/QSettings>
#include <QtCore/QMutex>
#include <QtCore/QPointer>
#include <QtCore/QScopedPointer>
#include <QtCore/QSharedPointer>
#include <QtCore/QMap>
#include <QtGui/QShortcut>
#include <QtCore/QDebug>

#include "ui_MainWindow.h"

#include "Manager.h"

#include "occ/OccManager.h"

class QtMainWindow : public QMainWindow, protected Ui::MainWindow 
{
    Q_OBJECT

public:
    QtMainWindow(QApplication *app);
    ~QtMainWindow();

    void setManager(QSharedPointer<OccManager> aManager);

protected:
     void closeEvent(QCloseEvent *event);

private slots:

    void fileNew();
    void fileOpen();
    void fileSaveAs();
    void fileImport();
    void fileExport();
    void applicationExit();

    void editDelete();

    void viewWireframe();
    void viewShaded();

    void viewFront();
    void viewBack();
    void viewLeft();
    void viewRight();
    void viewTop();
    void viewBottom();
    void viewIsometric();

    void viewColor();

    void selectVertices();
    void selectEdges();
    void selectWires();
    void selectFaces();
    void selectSolids();
    void selectAll();
    
    void decomposeShape();
    void fuseShape();
    void offsetShape();
    void buildShape();

    void meshSurface();
    void meshVolume();
    void meshStitch();

public slots:
    void selectionChanged();

private:

    // Application manager
    QSharedPointer<OccManager> m_manager;

    QApplication* m_app;

    size_t m_interval;
    size_t m_maxElem;

	QShortcut m_shortCutDelete;

    QString m_curDir;

    QProgressBar m_progressBar;

    void createProgressBar();
    void createActions();

    void readSettings();
    void writeSettings();

    void showProperties();

};

