// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include <QApplication>

#include "MainWindow.h"

int main(int argc, char *argv[])
{

     QApplication app(argc, argv);
     MainWindow window;

     window.show();

     return app.exec();
}

