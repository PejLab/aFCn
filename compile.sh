###############################################################################
#
# This script just compiles cython modules for the user 
#
###############################################################################

#compile the extensions
cd src && python3 setup.py build_ext --inplace 
