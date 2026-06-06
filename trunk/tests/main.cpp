// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "gtest/gtest.h"
#include <algorithm>
#include <string>
#include <vector>

bool g_saveFiles = false;

int main(int argc, char** argv)
{
    // Strip --save-files before passing args to GTest
    std::vector<char*> filteredArgs;
    for (int i = 0; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--save-files")
            g_saveFiles = true;
        else
            filteredArgs.push_back(argv[i]);
    }

    int filteredArgc = static_cast<int>(filteredArgs.size());
    ::testing::InitGoogleTest(&filteredArgc, filteredArgs.data());
    return RUN_ALL_TESTS();
}
