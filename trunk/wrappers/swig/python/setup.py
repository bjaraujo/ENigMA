import sys
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

with open('version.h', 'r') as fh:
    version = fh.read().split('"')[1]

setuptools.setup(
    name='ENigMApy',
    version=version,
    author='bjaraujo',
    author_email='',
    description='ENigMA - Extended Numerical Multiphysics Analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bjaraujo/ENigMA',
    packages=['ENigMA'],
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: Microsoft :: Windows' if sys.platform == 'win32' else 'Operating System :: POSIX :: Linux'
    ],
    python_requires='>=3.5'
)
