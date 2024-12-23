import sys
import os
import glob
from zipfile import ZipFile

with open('version.h', 'r') as fh:
    version = fh.read().split('"')[1]

if not os.path.exists('dist'):
    os.mkdir('dist')

if len(sys.argv) > 1:
    arch = sys.argv[1]

file = 'dist/ENigMAcs-' + version + '_' + arch + '.zip'
print(file)

pattern = 'ENigMA/csharp/*.cs' 
cs_files = glob.glob(pattern)

with ZipFile(file, 'w') as fh:
    fh.write('ENigMAcs.dll')
    fh.write('README.md')
    fh.write('LICENSE.txt')
    for cs in cs_files:
        fh.write(cs)
