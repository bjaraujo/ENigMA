
#include "OccStructureModel.h"

OccStructureModel::OccStructureModel()
{

}

OccStructureModel::OccStructureModel(QObject* parent) : QAbstractItemModel(parent)
{


}

OccStructureModel::~OccStructureModel()
{

}

int OccStructureModel::rowCount(const QModelIndex & /*parent*/) const
{
    return 3;
}

int OccStructureModel::columnCount(const QModelIndex & /*parent*/) const
{
    return 2;
}

Qt::ItemFlags OccStructureModel::flags(const QModelIndex &index) const
{
     if (!index.isValid())
         return 0;

     return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}


QModelIndex OccStructureModel::index(int row, int column, const QModelIndex &parent) const
{
 
    return QModelIndex();

}

QModelIndex OccStructureModel::parent(const QModelIndex &index) const
{

    return QModelIndex();

}

QVariant OccStructureModel::headerData(int section, Qt::Orientation orientation, int role) const
{

    if (role == Qt::DisplayRole)
    {
        return QString("Row%1, Column%2").arg("Property").arg("Value");
    }

    return QVariant();

}

QVariant OccStructureModel::data(const QModelIndex &index, int role) const
{

    if (role == Qt::DisplayRole)
    {
        return QString("Row%1, Column%2").arg(index.row() + 1).arg(index.column() + 1);
    }

    return QVariant();

}

