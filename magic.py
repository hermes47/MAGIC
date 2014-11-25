#!/usr/bin/env python3

import sys, getopt, os
from buildlipid import MAGIC
from config import ROOT_DIRECTORY, MTB_PDB_ROOT_NAMES, NAME_LENGTH, RANDOM_NAME_CHARS

def main(argv):
    ifpfile = '54A7.ifp'
    #pString = 'P8EO(B16HMCW5MCW5)'
    #pString = 'GGKT(7IRE)'
    #pString = 'P8EO(B16HMCW5FRFY)'
    #pString = 'P8EO(B16HFRFYFRFY)'
    pString = 'GLYM(PCHMMONMOLEM)'
    # generate the list of unique atom names in the various data directories
    # store in the config.MTB_PDB_ROOT_NAMES data point
    fileNames = []
    for folder in ('/Input/PDBFiles','/Input/MTBFiles','/Input/OVLFiles','/OutputFiles','/PickledFiles'):
        fileNames += os.listdir(path=ROOT_DIRECTORY+folder) 
    for file in fileNames:
        if len(file) < NAME_LENGTH: #short names are pointless
            continue
        if len(file) > NAME_LENGTH and file[4] != ".": #if not dot, is something other than an extension to the root name
            continue
        if file[0] not in RANDOM_NAME_CHARS:
            continue
        if file[1] not in RANDOM_NAME_CHARS:
            continue
        if file[2] not in RANDOM_NAME_CHARS:
            continue
        if file[3] not in RANDOM_NAME_CHARS:
            continue
        if file[0:4] not in MTB_PDB_ROOT_NAMES:
            MTB_PDB_ROOT_NAMES.append(file[0:4])
    try:
        opts, args = getopt.getopt(argv, "i:s:",["ifpfile=","parsestring="])
    except getopt.GetoptError:
        print("magic.py -s <StringToParse> -i <IFPFile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifpfile"):
            ifpfile = arg
        elif opt in ("-s", "--parsestring"):
            pString = arg
    MAGIC(ifpfile, pString)
            

if __name__ == "__main__":
    main(sys.argv[1:])
'''
Installing numpy for 3.4 is a requirement.
install distribute
- download 'distribute_setup.py' from http://python-distribute.org
- download https://pypi.python.org/packages/source/d/distribute/distribute-0.6.49.tar.gz
- make sure both in same directory
- run: python3 distribute_setup.py

install pip
- c&p https://raw.github.com/pypa/pip/master/contrib/get-pip.py into a file
- run: python3 get-pip.py

make sure pip is in your system path
- sudo ln -s /Library/Frameworks/Python.framework/Versions/3.4/bin/pip /usr/local/bin

install nose
- pip install nose

download NumPy from http://sourceforge.net/projects/numpy/files/NumPy/    most recent version 1.8.2
unarchive, go into folder
- export CC=clang
- export CXX=clang
- export FFLAGS=-ff2c
where X is most recent SDK in the folder
- export LDSHARED='clang -bundle -undefined dynamic_lookup -arch i386 -arch x86_64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.X.sdk -g'
- python3 setup.py build
- python3 setup.py install
'''