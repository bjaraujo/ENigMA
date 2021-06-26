
#ifndef ClipPolyhedron_H
#define ClipPolyhedron_H

#include <vector>

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>

#include "GeoCoordinate.hpp"
#include "GeoPlane.hpp"
#include "GeoPolyhedron.hpp"

using namespace ENigMA::geometry;

class vtkPolyData;

class vtkEventQtSlotConnect;
class vtkObject;

// Forward Qt class declarations
class Ui_ClipPolyhedron;

// Forward VTK class declarations
class vtkQtTableView;


class ClipPolyhedron : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    ClipPolyhedron(); 
    ~ClipPolyhedron();

    public slots:

        virtual void slotExit();

        virtual void slotKeyPressed(vtkObject *, unsigned long, void *, void *);

protected:

    protected slots:

private:
    
    size_t m_errors;

    double m_fraction;

    double m_size;

    CGeoPlane<double> m_plane;
    CGeoCoordinate<double> m_position;

    void drawClippedCube();

    vtkSmartPointer<vtkQtTableView>  m_tableView;
    vtkSmartPointer<vtkEventQtSlotConnect> m_connections;

    // Designer form
    Ui_ClipPolyhedron *ui;

};

#endif // ClipPolyhedron_H
