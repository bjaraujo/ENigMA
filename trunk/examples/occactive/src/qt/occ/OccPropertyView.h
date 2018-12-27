#pragma once

#include <QtCore/QSharedPointer>
#include <QtCore/QMap>
#include <QtGui/QStandardItemModel>
#include <QtGui/QStandardItem>
#include <QtGui/QTreeView>

#include "OccPropertyModel.h"

class OccPropertyView: public QTreeView
{
    Q_OBJECT
public:

    OccPropertyView();
    ~OccPropertyView();
    
private:

    OccPropertyModel m_model;

};


