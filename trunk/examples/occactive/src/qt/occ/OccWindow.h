#pragma once

#include <QtGui/QMdiSubWindow>
#include <QtCore/QList>

#include "OccWidget.h"

enum eMaterial {
    MAT_BRASS = 1,
    MAT_BRONZE = 2,
    MAT_COPPER = 3,
    MAT_GOLD = 4,
    MAT_PEWTER = 5,
    MAT_PLASTIC = 6,
    MAT_SILVER = 7,
    MAT_STEEL = 8
};

class OccWindow: public QMdiSubWindow
{
    Q_OBJECT
public:

    OccWindow();
    ~OccWindow();

    void eraseAll();

    void viewWireframe();
    void viewShaded();

    void viewFront();
    void viewBack();
    void viewLeft();
    void viewRight();
    void viewTop();
    void viewBottom();
    void viewIsometric();

    void setColor(QColor& aColor);
    void setMaterial(eMaterial aMaterial);

    void selectVertices();
    void selectEdges();
    void selectWires();
    void selectFaces();
    void selectSolids();
    void selectAll();

    void fitAll();

    Handle(AIS_InteractiveContext) context();

signals:
    void onSelectionChanged();

public slots:
    void selectionChanged();

private:

    Handle(Aspect_DisplayConnection) m_displayConnection;
    Handle(Graphic3d_GraphicDriver) m_graphicDriver;
    Handle(V3d_Viewer) m_viewer;
    Handle(AIS_InteractiveContext) m_context;

    OccWidget m_occWidget;

};


