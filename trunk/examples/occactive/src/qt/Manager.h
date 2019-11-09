#pragma once

#include <string>

#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtGui/QMdiSubWindow>
#include <QtGui/QTreeView>
#include <QtGui/QTableView>
#include <QtCore/QThread>

class IManager : public QObject
{   
public:

    virtual QMdiSubWindow* window() = 0;
    virtual QStandardItemModel* treeViewModel() = 0;

    virtual void reset() = 0;
    
    virtual void importGeometry(int nFormat, QString& aFileName) = 0;
    virtual void exportGeoemtry(int nFormat, QString& aFileName) = 0;

    virtual void importMesh(int nFormat, QString& aFileName) = 0;
    virtual void exportMesh(int nFormat, QString& aFileName) = 0;

    virtual void selectVertices() = 0;
    virtual void selectEdges() = 0;
    virtual void selectWires() = 0;
    virtual void selectFaces() = 0;
    virtual void selectSolids() = 0;
    virtual void selectAll() = 0;
    
    virtual void deleteSelected() = 0;
    
    virtual void viewWireframe() = 0;
    virtual void viewShaded() = 0;

    virtual void viewFront() = 0;
    virtual void viewBack() = 0;
    virtual void viewLeft() = 0;
    virtual void viewRight() = 0;
    virtual void viewTop() = 0;
    virtual void viewBottom() = 0;
    virtual void viewIsometric() = 0;

    virtual void setColor(QColor aColor) = 0;
    
    virtual void decomposeShape() = 0;
    virtual void fuseShape() = 0;
    virtual void offsetShape(double offsetValue) = 0;
    virtual void buildShape() = 0;

    virtual void meshSurface(double meshSize) = 0;
    virtual void meshVolume(double meshSize) = 0;

    virtual void getProperties(QList<QString>& sPropertyNames, QList<QVariant>& sPropertyValues) = 0;

};


