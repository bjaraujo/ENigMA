#pragma once

#include <QtCore/QSharedPointer>
#include <QtCore/QList>
#include <QtGui/QStandardItemModel>
#include <QtGui/QStandardItem>
#include <QtGui/QTreeView>

#include "OccStructureModel.h"

class OccStructureView: public QTreeView
{
    Q_OBJECT
public:

    OccStructureView();
    ~OccStructureView();

private:

    OccStructureModel m_model;

};


