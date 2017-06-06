
#ifndef TetrahedronIntersection_H
#define TetrahedronIntersection_H

#include <vector>

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>

#include "GeoCoordinate.hpp"
#include "GeoLine.hpp"
#include "GeoTetrahedron.hpp"
#include "MshBasicMesher.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::mesh;

class vtkPolyData;

class vtkEventQtSlotConnect;
class vtkObject;
class vtkActor;
class vtkAxesActor;

// Forward Qt class declarations
class Ui_TetrahedronIntersection;

// Forward VTK class declarations
class vtkQtTableView;


class TetrahedronIntersection : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    TetrahedronIntersection(); 
    ~TetrahedronIntersection();

    public slots:

        virtual void slotOpenFile();
        virtual void slotExit();

        virtual void slotKeyPressed(vtkObject *, unsigned long, void *, void *);

protected:

    protected slots:

private:
    
    CMshMesh<double> m_mesh;
    CGeoCoordinate<double> m_point;

    void buildMesh();
    
    void drawMesh();
    void drawTetrahedron();
    void drawPoint();

    vtkSmartPointer<vtkQtTableView>  m_tableView;
    vtkSmartPointer<vtkEventQtSlotConnect> m_connections;

    // Designer form
    Ui_TetrahedronIntersection *ui;

    vtkSmartPointer<vtkActor> m_meshActor;
    vtkSmartPointer<vtkActor> m_pointActor;
    vtkSmartPointer<vtkActor> m_tetrahedronActor;
    vtkSmartPointer<vtkAxesActor> m_axesActor;

};

#endif // TetrahedronIntersection_H
