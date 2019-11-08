#! /usr/bin/python3
import os, sys
thisdir = os.path.dirname(__file__)
libdir = os.path.join(thisdir, '../modules')
print(thisdir)
print(libdir)
if libdir not in sys.path:
	print(f'\tadding {libdir} to sys.path')
	sys.path.insert(0, libdir)