
#include "OccPropertyModel.h"

OccPropertyModel::OccPropertyModel()
{

}

OccPropertyModel::OccPropertyModel(QObject* parent) : QAbstractItemModel(parent)
{

}

OccPropertyModel::~OccPropertyModel()
{

}

int OccPropertyModel::rowCount(const QModelIndex & /*parent*/) const
{
    return 3;
}

int OccPropertyModel::columnCount(const QModelIndex & /*parent*/) const
{
    return 2;
}

Qt::ItemFlags OccPropertyModel::flags(const QModelIndex &index) const
{
     if (!index.isValid())
         return 0;

     return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

QModelIndex OccPropertyModel::index(int row, int column, const QModelIndex &parent) const
{
 
    return QModelIndex();

}

QModelIndex OccPropertyModel::parent(const QModelIndex &index) const
{

    return QModelIndex();

}

QVariant OccPropertyModel::headerData(int section, Qt::Orientation orientation, int role) const
{

    if (role == Qt::DisplayRole)
    {
        return QString("Row%1, Column%2").arg("Property").arg("Value");
    }

    return QVariant();

}

QVariant OccPropertyModel::data(const QModelIndex &index, int role) const
{

    if (role == Qt::DisplayRole)
    {
        return QString("Row%1, Column%2").arg(index.row() + 1).arg(index.column() + 1);
    }

    return QVariant();

}

