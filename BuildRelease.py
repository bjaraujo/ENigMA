
import os
import shutil
   
# Read version info
fi = open('trunk/src/version.h', 'r')
strLines = fi.readlines()
fi.close()

intMajorNumber = 0
intMinorNumber = 0
intReleaseNumber = 0
intBuildNumber = 0

strNewVersion = str(intMajorNumber) + '.' + str(intMinorNumber) + '.' + str(intReleaseNumber) + '.' + str(intBuildNumber)

# Increment version
incVersionInput = input('Increment version number [Y/n]?')

incVersion = True 
if incVersionInput == 'n':
    incVersion = False

fo = open('trunk/src/version.h', 'w')

for strLine in strLines:

    if 'APP_VERSION_INFO' in strLine.upper(): 
        strPrevVersion = strLine.split(' ')[-1].rstrip().replace('"', '')
                
        intMajorNumber = int(strPrevVersion.split('.')[0])
        intMinorNumber = int(strPrevVersion.split('.')[1])
        intReleaseNumber = int(strPrevVersion.split('.')[2])
        intBuildNumber = int(strPrevVersion.split('.')[3])
        
        if incVersion:
            intReleaseNumber = intReleaseNumber + 1
            
            if intReleaseNumber > 9:
                intMinorNumber = intMinorNumber + 1
                intReleaseNumber = 0

            if intMinorNumber > 9:
                intMajorNumber = intMajorNumber + 1
                intMinorNumber = 0
            
        strNewVersion = str(intMajorNumber) + '.' + str(intMinorNumber) + '.' + str(intReleaseNumber) + '.' + str(intBuildNumber)
        
        strLine = strLine.replace(strPrevVersion, strNewVersion)
                
    fo.write(strLine)
    
fo.close()

# Build packages
print('------------ Building release ------------')
print('version: ' + strNewVersion)

os.system('pause')

if not os.path.exists('build'):
    os.mkdir('build')
    
os.system('cmake -Bbuild trunk -G "Visual Studio 16 2019" -A Win32 -DENIGMA_BUILD_UNIT_TESTS:BOOL=OFF -DENIGMA_BUILD_WRAPPERS_SWIG:BOOL=ON -DWRAP_SWIG_PYTHON:BOOL=ON -DWRAP_SWIG_CSHARP:BOOL=ON -DSWIG_EXECUTABLE=D:/Libraries/Swig/swigwin-4.0.1/swig.exe')
os.system('cmake --build build --config Release')

if not os.path.exists('release'):
    os.mkdir('release')

pyReleaseFolder = 'release/ENigMApy_' + strNewVersion + '_win32'
if not os.path.isdir(pyReleaseFolder):
    os.mkdir(pyReleaseFolder)

# Copy files
setup = ['import setuptools', \
         '', \
         'with open("README.md", "r") as fh:', \
         '    long_description = fh.read()', \
         '', \
         'setuptools.setup(', \
         '    name="ENigMApy",', \
         '    version="' + strNewVersion + '",', \
         '    author="baraujo",', \
         '    author_email="",', \
         '    description="ENigMA - Extended Numerical Multiphysics Analysis",', \
         '    long_description=long_description,', \
         '    long_description_content_type="text/markdown",', \
         '    url="https://github.com/bjaraujo/ENigMA",', \
         '    packages=["ENigMA"],', \
         '    include_package_data=True,', \
         '    classifiers=[', \
         '        "Programming Language :: Python :: 3",', \
         '        "License :: OSI Approved :: MIT License",', \
         '        "Operating System :: OS Independent",', \
         '    ],', \
         '    python_requires=">=3.6",', \
         ')']

with open(pyReleaseFolder + '/setup.py', 'w') as f:
    for l in setup:
        f.write(l + '\n')

with open(pyReleaseFolder + '/MANIFEST.in', 'w') as f:
    f.write('recursive-include ENigMA *.pyd\n')

if not os.path.isdir(pyReleaseFolder + '/ENigMA'):
    os.mkdir(pyReleaseFolder + '/ENigMA')
        
shutil.copy2('build/bin/_ENigMA.pyd', pyReleaseFolder + '/ENigMA/_ENigMA.pyd')
shutil.copy2('build/bin/ENigMA.py', pyReleaseFolder + '/ENigMA/_ENigMA.py')

# Copy license
shutil.copy2('LICENSE.txt', pyReleaseFolder + '/LICENSE.txt')
shutil.copy2('README.md', pyReleaseFolder + '/README.md')

# Create tag
#os.system('git commit -a -m v' + strNewVersion)
#os.system('git tag v' + strNewVersion)
#os.system('git push --tags')

