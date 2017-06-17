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

#include "TetrahedronIntersection.h"

extern int qInitResources_icons();

int main(int argc, char** argv)
{
  
  // QT Stuff
  QApplication app(argc, argv);
  
  qInitResources_icons();
  
  TetrahedronIntersection myTetrahedron3d;
  myTetrahedron3d.show();
  
  return app.exec();
}
