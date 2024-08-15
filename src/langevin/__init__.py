import os, sys
#Taken from https://stackoverflow.com/a/49375740
dirName = os.path.dirname(os.path.realpath(__file__))
if dirName not in sys.path:
    sys.path.append(dirName)

#Interesting to note that this script fails, but this is correct for loading
#the package
from .langevin import *