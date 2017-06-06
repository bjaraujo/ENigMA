
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
class vtkScalarBarActor;
class vtkRenderer;
class vtkFeatureEdges;
class vtkUnstructuredGrid;

// Forward Qt class declarations
class Ui_FemFlow;

// Forward VTK class declarations
class vtkQtTableView;

#include "MshMesh.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::pde;
using namespace ENigMA::post;

class FemFlow : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    FemFlow();
    ~FemFlow();

    void init();

    public slots:

    virtual void slotExit();
    virtual void slotSolve();

protected:

    protected slots :

private:

    void solve1();
    void solve2();

    vtkSmartPointer<vtkRenderer> m_renderer;

    vtkSmartPointer<vtkActor> m_resultsActor;
    vtkSmartPointer<vtkScalarBarActor> m_scalarBarActor;

    // Designer form
    Ui_FemFlow *ui;

};

