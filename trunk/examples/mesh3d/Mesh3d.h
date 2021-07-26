
#ifndef Mesh3d_H
#define Mesh3d_H

#include <vector>

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>
#include <QtCore/QMutex>

#include "MshBasicMesher.hpp"
#include "MshTetrahedronMesher.hpp"
#include "PosGmsh.hpp"
#include "PosQuickMesh.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::mesh;
using namespace ENigMA::post;

class vtkPolyData;

class vtkEventQtSlotConnect;
class vtkObject;

// Forward Qt class declarations
class Ui_Mesh3d;

// Forward VTK class declarations
class vtkQtTableView;


class Mesh3D : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    Mesh3D(QString fileName); 
    ~Mesh3D();

    public slots:

        virtual void slotOpenFile();
        virtual void slotExit();

        virtual void slotKeyPressed(vtkObject *, unsigned long, void *, void *);

    int  update(bool updateData);

protected:

    protected slots:

private:

    QString m_fileName;

    CMshTetrahedronMesher<double> m_mesher;

    bool generateMesh(const unsigned int nu, const unsigned int nv, const unsigned int nw);
    bool loadMesh();
    void drawMesh();

    vtkSmartPointer<vtkQtTableView>  m_tableView;
    vtkSmartPointer<vtkEventQtSlotConnect> m_connections;

    // Designer form
    Ui_Mesh3d *ui;

    QMutex m_mutex;

};

#endif // Mesh3d_H

