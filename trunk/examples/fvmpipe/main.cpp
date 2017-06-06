// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

// QT includes
#include <QApplication>

#include "FvmPipe.h"

extern int qInitResources_icons();

int main(int argc, char** argv)
{

    // QT Stuff
    QApplication app(argc, argv);

    qInitResources_icons();

    FvmPipe myFvmPipe;
    myFvmPipe.show();

    myFvmPipe.init();

    return app.exec();

}
