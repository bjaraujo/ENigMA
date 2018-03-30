
#ifndef Mesh2d_H
#define Mesh2d_H

#include <vector>

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>
#include <QtCore/QMutex>

#include "MshBasicMesher.hpp"
#include "MshTriangleMesher.hpp"
#include "MshQuadrilateralMesher.hpp"
#include "PosGmsh.hpp"
#include "PosQuickMesh.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::mesh;
using namespace ENigMA::post;

class vtkPolyData;

class vtkEventQtSlotConnect;
class vtkObject;

// Forward Qt class declarations
class Ui_Mesh2d;

// Forward VTK class declarations
class vtkQtTableView;

class Mesh2D : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    Mesh2D(); 
    ~Mesh2D();

    public slots:

        virtual void slotOpenFile();
        virtual void slotExit();

        virtual void slotKeyPressed(vtkObject *, unsigned long, void *, void *);

    int  updateData(int );

protected:

    protected slots:

private:
    
    //CMshTriangleMesher<double> m_mesher;
    CMshQuadrilateralMesher<double> m_mesher;

    bool generateMesh(const unsigned int nu, const unsigned int nv);
    void drawMesh();

    vtkSmartPointer<vtkQtTableView>  m_tableView;
    vtkSmartPointer<vtkEventQtSlotConnect> m_connections;

    // Designer form
    Ui_Mesh2d *ui;

    QMutex m_mutex;

};

#endif // Mesh2d_H

