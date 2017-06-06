
#ifndef TriangleIntersection_H
#define TriangleIntersection_H

#include <vector>

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>

#include "GeoCoordinate.hpp"
#include "GeoLine.hpp"
#include "GeoTriangle.hpp"

using namespace ENigMA::geometry;

class vtkPolyData;

class vtkEventQtSlotConnect;
class vtkObject;

// Forward Qt class declarations
class Ui_TriangleIntersection;

// Forward VTK class declarations
class vtkQtTableView;


class TriangleIntersection : public QMainWindow
{
    Q_OBJECT

public:

    // Constructor/Destructor
    TriangleIntersection(); 
    ~TriangleIntersection();

    public slots:

        virtual void slotOpenFile();
        virtual void slotExit();

        virtual void slotKeyPressed(vtkObject *, unsigned long, void *, void *);

protected:

    protected slots:

private:
    
    CGeoTriangle<double> m_triangle1;
    CGeoTriangle<double> m_triangle2;

    CGeoCoordinate<double> m_position;
    CGeoCoordinate<double> m_p1, m_p2, m_p3, m_p4, m_p5, m_p6;

    void drawTriangles();

    vtkSmartPointer<vtkQtTableView>  m_tableView;
    vtkSmartPointer<vtkEventQtSlotConnect> m_connections;

    // Designer form
    Ui_TriangleIntersection *ui;

};

#endif // TriangleIntersection_H
