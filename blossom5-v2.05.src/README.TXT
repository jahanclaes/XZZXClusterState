MINIMAL PYTHON CTYPES WRAPPER OF BLOSSOM V
##########################################

ORIGINAL:
- BLOSSOM V - implementation of Edmonds' algorithm for computing a minimum cost perfect matching in a graph
  Version 2.05
  http://pub.ist.ac.at/~vnk/software.html

CHANGES:
- Add pypm.cpp (to expose basic API to ctypes)
- Add pypm_example.py (to show how to access API via Python)
- Rename README.TXT to README.TXT.ORIG
- Add README.TXT (this file)
- Rename Makefile to Makefile.ORIG
- Add Makefile (to build shared library)

COMPILATION:
- Type "make".

USAGE:
- See pypm_example.py for example.

FURTHER INFORMATION:
- See README.TXT.ORIG for original readme file.
