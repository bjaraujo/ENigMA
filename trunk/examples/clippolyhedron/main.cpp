// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

// QT includes
#include <QApplication>

#include "ClipPolyhedron.h"

int main(int argc, char** argv)
{

    // QT Stuff
    QApplication app(argc, argv);

    ClipPolyhedron myClipPolyhedron;
    myClipPolyhedron.show();

    return app.exec();
}
