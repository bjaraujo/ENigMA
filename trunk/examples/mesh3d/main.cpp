
// QT includes
#include <QApplication>

#include "Mesh3d.h"

extern int qInitResources_icons();

int main(int argc, char** argv)
{

    // QT Stuff
    QApplication app(argc, argv);

    QString fileName;

    if (argc > 1)
        fileName = QString::fromStdString(argv[1]);

    qInitResources_icons();

    Mesh3D myMesh3D(fileName);
    myMesh3D.show();

    return app.exec();
}
