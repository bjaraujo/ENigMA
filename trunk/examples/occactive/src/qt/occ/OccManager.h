#pragma once

#include <QtCore/QObject>
#include <QtGui/QStandardItemModel>
#include <QtGui/QItemSelectionModel>

#include "OccWindow.h"
#include "OccStructureView.h"
#include "OccPropertyView.h"

#include "../../occ/OccTranslator.h"
#include "../../occ/OccMesher.h"

#include "../Manager.h"

class gp_Pnt;
class Poly_Triangle;

class OccManager : public IManager
{
    Q_OBJECT
public:

    OccManager();
    ~OccManager();

    static QSharedPointer<OccManager> instance()
    {

        if (m_instance == NULL)
        {
            m_instance = QSharedPointer<OccManager>(new OccManager());
        }

        return m_instance;

    }

    OccWindow* window();
    QStandardItemModel* treeViewModel();
    QItemSelectionModel* treeViewSelectionModel();

    void reset();

    void displayShape(const TopoDS_Shape& aShape);

    void importGeometry(int nFormat, QString& aFileName);
    void exportGeoemtry(int nFormat, QString& aFileName);

    void importMesh(int nFormat, QString& aFileName);
    void exportMesh(int nFormat, QString& aFileName);

    void selectVertices();
    void selectEdges();
    void selectWires();
    void selectFaces();
    void selectSolids();
    void selectAll();
    
    void deleteSelected();
    
    void selectShapes(QModelIndexList selected);

    void viewWireframe();
    void viewShaded();

    void viewFront();
    void viewBack();
    void viewLeft();
    void viewRight();
    void viewTop();
    void viewBottom();
    void viewIsometric();

    void setColor(QColor aColor);

    void decomposeShape();
    void fuseShape();
    void offsetShape(double offsetValue);
    void buildShape();

    void meshSurface(double meshSize);
    void meshVolume(double meshSize);
    void meshStitch();

    void getProperties(QList<QString>& sPropertyNames, QList<QVariant>& sPropertyValues);

public slots:
    void selectionChanged();
    void treeViewSelectionChanged(const QItemSelection& selected, const QItemSelection& deselected);

signals:
    void onSelectionChanged();

protected:
    static QSharedPointer<OccManager> m_instance;

private:

    OccWindow m_occWindow;

    QStandardItemModel m_treeViewModel;
    QItemSelectionModel m_selectionModel;

    QStandardItem* m_headerItem;
    QStandardItem* m_modelRoot;

    OccTranslator m_translate;
    OccMesher m_mesher;

    QMap<QModelIndex, TopoDS_Shape> m_shapeMap;

};
