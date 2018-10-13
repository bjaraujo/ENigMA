// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include <time.h>

#include "MshMesh.hpp"
#include "PosGmsh.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::pde;

int main(int argc, char *argv[])
{

    CPdeField<double> aPdeField1, aPdeField2;
    CPosGmsh<double> aPosGmsh;

    //aPosGmsh.load(aPdeField1, "C:\\Temp\\merge\\telemac_domain2.msh");
    aPosGmsh.load(aPdeField1, "C:\\Temp\\merge\\artemis_domain1.msh");
    aPosGmsh.load(aPdeField2, "C:\\Temp\\merge\\blf_patch.msh");

    //aPosGmsh.save(aPdeField1, "C:\\Temp\\merge\\telemac_domain2_renumbered.msh", "tris");
    aPosGmsh.save(aPdeField1, "C:\\Temp\\merge\\artemis_domain1_renumbered.msh", "tris");

    aPdeField1.mesh().addMesh(aPdeField2.mesh());
    
    aPdeField1.mesh().mergeNodes(1E-6);

    //aPosGmsh.save(aPdeField1, "C:\\Temp\\merge\\telemac_domain2_merged.msh", "tris");
    aPosGmsh.save(aPdeField1, "C:\\Temp\\merge\\artemis_domain1_merged.msh", "tris");

}
