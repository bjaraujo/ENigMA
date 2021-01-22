
import os
import sys
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
incVersionInput = input('Increment version number [y/N]? ') or 'N'

incVersion = False
if incVersionInput.upper() == 'Y':
    incVersion = True

print('1 - 32bit')
print('2 - 64bit')
version32or64Bit = input('Which version [2]? ') or '2'

version64 = True
if version32or64Bit == '1':
    version64 = False
elif version32or64Bit == '2':
    version64 = True

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

print(sys.executable)

# Build packages
print('------------ Building release ------------')
print('version: ' + strNewVersion)
if version64:
    print('64 bit')
else:
    print('32 bit')

os.system('pause')

if not os.path.exists('build'):
    os.mkdir('build')

os.chdir('build')

swigExecuable = 'D:/Libraries/Swig/swigwin-4.0.2/swig.exe'
unitTests = 'OFF'

if sys.platform == 'win32':
    if version64:
        os.system('cmake ../trunk -G "Visual Studio 16 2019" -A x64 -DENIGMA_BUILD_UNIT_TESTS:BOOL=' + unitTests + ' -DENIGMA_BUILD_WRAPPERS_SWIG:BOOL=ON -DWRAP_SWIG_PYTHON:BOOL=ON -DWRAP_SWIG_CSHARP:BOOL=OFF -DSWIG_EXECUTABLE=' + swigExecuable)
    else:
        os.system('cmake ../trunk -G "Visual Studio 16 2019" -A Win32 -DENIGMA_BUILD_UNIT_TESTS:BOOL=' + unitTests + ' -DENIGMA_BUILD_WRAPPERS_SWIG:BOOL=ON -DWRAP_SWIG_PYTHON:BOOL=ON -DWRAP_SWIG_CSHARP:BOOL=OFF -DSWIG_EXECUTABLE=' + swigExecuable)
else:
    os.system('cmake ../trunk -G "Ninja" -DENIGMA_BUILD_UNIT_TESTS:BOOL=' + unitTests + ' -DENIGMA_BUILD_WRAPPERS_SWIG:BOOL=ON -DWRAP_SWIG_PYTHON:BOOL=ON')

os.system('cmake --build . --config Release')

# https://packaging.python.org/tutorials/packaging-projects/ 
pyReleaseFolderTemp = 'bin/ENigMApy'
if not os.path.isdir(pyReleaseFolderTemp):
    os.mkdir(pyReleaseFolderTemp)

# Copy files
setup = ['import setuptools', \
         '', \
         'with open("README.md", "r") as fh:', \
         '    long_description = fh.read()', \
         '', \
         'setuptools.setup(', \
         '    name="ENigMApy",', \
         '    version="' + strNewVersion + '",', \
         '    author="bjaraujo",', \
         '    author_email="",', \
         '    description="ENigMA - Extended Numerical Multiphysics Analysis",', \
         '    long_description=long_description,', \
         '    long_description_content_type="text/markdown",', \
         '    url="https://github.com/bjaraujo/ENigMA",', \
         '    packages=["ENigMA"],', \
         '    include_package_data=True,', \
         '    classifiers=[', \
         '        "Programming Language :: Python :: 3",', \
         '        "License :: OSI Approved :: GNU General Public License (GPL)",', \
         '        "Operating System :: Microsoft :: Windows",' if sys.platform == 'win32' else '        "Operating System :: POSIX :: Linux",', \
         '    ],', \
         '    python_requires=">=3.6",', \
         ')']

with open(pyReleaseFolderTemp + '/setup.py', 'w') as f:
    for l in setup:
        f.write(l + '\n')

if sys.platform == 'win32':
    with open(pyReleaseFolderTemp + '/MANIFEST.in', 'w') as f:
        f.write('recursive-include ENigMA *.pyd\n')
else:
    with open(pyReleaseFolderTemp + '/MANIFEST.in', 'w') as f:
        f.write('recursive-include ENigMA *.so\n')

if not os.path.isdir(pyReleaseFolderTemp + '/ENigMA'):
    os.mkdir(pyReleaseFolderTemp + '/ENigMA')

if sys.platform == 'win32':
    open(pyReleaseFolderTemp + '/ENigMA/__init__.py', 'w').close()
    shutil.copy2('bin/_ENigMA.pyd', pyReleaseFolderTemp + '/ENigMA/_ENigMA.pyd')
    shutil.copy2('bin/ENigMA.py', pyReleaseFolderTemp + '/ENigMA/ENigMA.py')
else:
    open(pyReleaseFolderTemp + '/ENigMA/__init__.py', 'w').close()
    shutil.copy2('bin/ENigMA.so', pyReleaseFolderTemp + '/ENigMA/_ENigMA.so')
    shutil.copy2('bin/ENigMA.py', pyReleaseFolderTemp + '/ENigMA/ENigMA.py')

shutil.copy2('../LICENSE.txt', pyReleaseFolderTemp + '/LICENSE.txt')
shutil.copy2('../README.md', pyReleaseFolderTemp + '/README.md')

os.chdir(pyReleaseFolderTemp)
os.system(sys.executable + ' setup.py bdist_wheel')
os.chdir('../..')

if not os.path.exists('../release'):
    os.mkdir('../release')

# Python module
if sys.platform == 'win32':
    if version64:
        pyReleaseFolder = '../release/ENigMApy_' + strNewVersion + '_x64'
    else:
        pyReleaseFolder = '../release/ENigMApy_' + strNewVersion + '_win32'
        
    if not os.path.isdir(pyReleaseFolder):
        os.mkdir(pyReleaseFolder)

    if version64:
        shutil.copy2('bin/ENigMApy/dist/ENigMApy-' + strNewVersion + '-py3-none-any.whl', pyReleaseFolder + '/ENigMApy-' + strNewVersion + '-py3-none-win_amd64.whl')
    else:
        shutil.copy2('bin/ENigMApy/dist/ENigMApy-' + strNewVersion + '-py3-none-any.whl', pyReleaseFolder + '/ENigMApy-' + strNewVersion + '-py3-none-win32.whl')
else:
    pyReleaseFolder = '../release/ENigMApy_' + strNewVersion + '_linux'
    if not os.path.isdir(pyReleaseFolder):
        os.mkdir(pyReleaseFolder)

    shutil.copy2('bin/ENigMApy/dist/ENigMApy-' + strNewVersion + '-py3-none-any.whl', pyReleaseFolder + '/ENigMApy-' + strNewVersion + '-py3-none-linux_x86_64.whl')

# CSharp module
if sys.platform == 'win32':
    if version64:
        csReleaseFolder = '../release/ENigMAcs_' + strNewVersion + '_x64'
    else:
        csReleaseFolder = '../release/ENigMAcs_' + strNewVersion + '_win32'
        
    if not os.path.isdir(csReleaseFolder):
        os.mkdir(csReleaseFolder)

    shutil.copy2('bin/ENigMAcs.dll', csReleaseFolder + '/ENigMAcs.dll')
    shutil.copy2('../LICENSE.txt', csReleaseFolder + '/LICENSE.txt')
    shutil.copy2('../README.md', csReleaseFolder + '/README.md')
    
# Create tag

if incVersion:
    os.system('git commit -a -m v' + strNewVersion)
    os.system('git tag v' + strNewVersion)
    os.system('git push --tags')


