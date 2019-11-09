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

#include "FvmChannel.h"

extern int qInitResources_icons();

int main(int argc, char** argv)
{

    // QT Stuff
    QApplication app(argc, argv);

    qInitResources_icons();

    FvmChannel myFvmChannel;
    myFvmChannel.show();

    myFvmChannel.init();

    return app.exec();

}
