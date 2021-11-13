# *******************************************************************
# Main Makefile
# Christian Wolf, christian.wolf@insa-lyon.fr
#
# Possible values for the global CFLAGS settings:
#
# ************ General 
#
# CHECK_CODE			Include checks for data integrity
# PEDANTIC_CHECK_CODE	Include additional checks (very slow)
# PEDANTIC_MATRIX_FINITE_CHECK
#                       Check whether values in the Matrix<T>
#                       class are finite or note (NaN, +/-Inf)
# PEDANTIC_MATRIX_FINITE_CHECK_ISFLOAT
#                       Supposes that the matrix values are floats
#                       if not defined ---> doubles
# WRITE_DEBUG_IMAGES	Write intermediate debug images to disk
# HAVE_VISUALIZE		Do we have the visualization module? 
# UNIX_ARCHITECTURE 	Enable UNIX goodies
# DEBUG
# ENABLE_PROFILING
# DO_FFT				(Add the dependency to the FFT group)
# VISUALBUG				Enables workarounds for Bugs in Visual C++
# TELESUN				Disables some code which does not work on HP-UX
#
# ********************************************************************

MACROS		 = -DCHECK_CODE -DUSE_LIBJPEG -DUNIX_ARCHITECTURE -DDEBUG -DHAVE_VISUALIZE -DTEXTDETECT_EVALUATION -DVISUALIZE  -DDBLMRF_WITH_POTTS_AND_GRAPHCUT -DHAVE_LIBXML2 

# ---- Uncomment this for full optimization
OPT			= -O3 -fsingle-precision-constant
LOPT        	=

# ---- Compiler flags
CFLAGS_GLOBAL 	= -Wall $(OPT) -I../PortVisualC++ -I/usr/include/qt4 -I/usr/include/libxml2 $(MACROS)

# ---- Linker flags
LFLAGS_GLOBAL	= $(LOPT) -L/usr/lib/qt/lib -lm -lstdc++ -ljpeg 

INSTALLDIR		= ../bin

all:
			../tools/makeexecutable exmrf $(INSTALLDIR)/mrf $(CFLAGS_GLOBAL) , $(LFLAGS_GLOBAL)

addpath:			
			echo $(PATH)

clean:
			find . -name "Makefile_[0-9]*.tmp" -exec rm -f {} \; -print
			find . -name "core.[0-9]*" -exec rm -f {} \; -print
			find . -name "*.cpp~" -exec rm -f {} \; -print
			find . -name "*.h~" -exec rm -f {} \; -print
			find . -name Makefile.makedefault -exec rm -f {} \; -print
			find . -name *.o -exec rm -f {} \; -print
			for i in $(EX_MF) ; do make -C ex$$i  clean ; done


