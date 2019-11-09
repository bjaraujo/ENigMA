
import os
import shutil

# Build platform
build = 'msvc15-win64'

configuration = 'Release'
    
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

# Build
print('------------ Building release ------------')
print('version: ' + strNewVersion)
print('build: ' + build)

os.system("pause")
os.system('"C:\Program Files (x86)/Microsoft Visual Studio/2017/BuildTools/MSBuild/15.0/Bin/MSBuild.exe" build/' + build + '/wrappers/swig/python/_ENigMA.vcxproj /p:Configuration=' + configuration + ' /t:Rebuild')

if not os.path.exists('releases'):
    os.mkdir('releases')

strNewFolder = 'releases/' + build
if not os.path.isdir(strNewFolder):
	os.mkdir(strNewFolder)

strNewFolder = 'releases/' + build + '/ENigMA_python3_64bit_' + strNewVersion
if not os.path.isdir(strNewFolder):
	os.mkdir(strNewFolder)

# Copy files
shutil.copy2('build/' + build + '/wrappers/swig/python/' + configuration + '/_ENigMA.pyd', strNewFolder + '/_ENigMA.pyd')
shutil.copy2('build/' + build + '/wrappers/swig/python/' + configuration + '/ENigMA.py', strNewFolder + '/ENigMA.py')
    
# Copy license
shutil.copy2('LICENSE.txt', strNewFolder + '/LICENSE.txt')

# Create tag
os.system('git commit -a -m v' + strNewVersion)
os.system('git tag v' + strNewVersion)
os.system('git push --tags')




