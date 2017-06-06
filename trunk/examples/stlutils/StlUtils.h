
#pragma once

#include <vector>

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>
#include <QFileDialog>
#include <QDebug>

class vtkPolyData;

class vtkEventQtSlotConnect;
class vtkObject;
class vtkSTLReader;
class vtkPolyDataMapper;
class vtkActor;
class vtkAxesActor;
class vtkRenderer;
class vtkFeatureEdges;

// Forward Qt class declarations
class Ui_StlUtils;

// Forward VTK class declarations
class vtkQtTableView;


class StlUtils : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    StlUtils(); 
    ~StlUtils();

    public slots:

        virtual void slotOpenFile();
        virtual void slotExit();

        virtual void slotKeyPressed(vtkObject *, unsigned long, void *, void *);

protected:

    protected slots:

private:

    void readStlFile();
    void drawStlFile();
    
    QString m_stlFileName;

    vtkSmartPointer<vtkFeatureEdges> m_boundaryEdges;

    vtkSmartPointer<vtkPolyData> m_points;

    vtkSmartPointer<vtkPolyDataMapper> m_stlMapper;
    vtkSmartPointer<vtkPolyDataMapper> m_edgeMapper;
    vtkSmartPointer<vtkPolyDataMapper> m_pointMapper;

    vtkSmartPointer<vtkActor> m_stlActor;
    vtkSmartPointer<vtkActor> m_edgeActor;
    vtkSmartPointer<vtkActor> m_pointActor;

    vtkSmartPointer<vtkAxesActor> m_axesActor;

    vtkSmartPointer<vtkRenderer> m_renderer;

    vtkSmartPointer<vtkSTLReader> m_stlReader;

    vtkSmartPointer<vtkQtTableView> m_tableView;
    vtkSmartPointer<vtkEventQtSlotConnect> m_connections;

    // Designer form
    Ui_StlUtils *ui;

};
