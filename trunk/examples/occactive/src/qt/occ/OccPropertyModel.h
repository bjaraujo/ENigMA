#pragma once

#include <QtCore/QSharedPointer>
#include <QtCore/QMap>
#include <QtCore/QAbstractItemModel>
#include <QtGui/QStandardItemModel>
#include <QtGui/QStandardItem>

class OccPropertyModel: public QAbstractItemModel
{
    Q_OBJECT
public:

    OccPropertyModel();
    OccPropertyModel(QObject* parent);
    ~OccPropertyModel();

    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;

    Qt::ItemFlags flags(const QModelIndex &index) const;
    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &index) const;

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;

};


