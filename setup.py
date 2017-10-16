import sys
from numpy.distutils.core import Extension, setup

__author__ = "Lars Andersen Bratholm"
__copyright__ = "Copyright 2017"
__credits__ = ["Lars Andersen Bratholm (2017) https://github.com/larsbratholm/fns"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Lars Andersen Bratholm"
__email__ = "lars.bratholm@bristol.ac.uk"
__status__ = "Alpha"
__description__ = "Furthest Neighbour Search"
__url__ = "https://github.com/larsbratholm/fns"


FORTRAN = "f90"

# GNU (default)
COMPILER_FLAGS = ["-fopenmp", "-m64", "-march=native", "-fPIC", "-Ofast", "-ffast-math", "-funroll-loops",
                    "-Wno-maybe-uninitialized", "-Wno-unused-function", "-Wno-cpp"]#, "-fcheck=all"]
#LINKER_FLAGS = ["-lgomp"]
#MATH_LINKER_FLAGS = ["-lblas", "-llapack"]

# For clang without OpenMP: (i.e. most Apple/mac system)
if sys.platform == "darwin" and all(["gnu" not in arg for arg in sys.argv]):
    COMPILER_FLAGS = ["-O3", "-m64", "-march=native", "-fPIC"]
    LINKER_FLAGS = []
    MATH_LINKER_FLAGS = ["-lblas", "-llapack"]


# Intel
if any(["intelem" in arg for arg in sys.argv]):
    COMPILER_FLAGS = ["-xHost", "-O3", "-axAVX", "-qopenmp"]
    LINKER_FLAGS = ["-liomp5", " -lpthread", "-lm", "-ldl"]
    MATH_LINKER_FLAGS = ["-L${MKLROOT}/lib/intel64", "-lmkl_rt"]


# UNCOMMENT TO FORCE LINKING TO MKL with GNU compilers:
LINKER_FLAGS = ["-lgomp", " -lpthread", "-lm", "-ldl"]
MATH_LINKER_FLAGS = ["-L${MKLROOT}/lib/intel64", "-lmkl_rt"]


#ext_frepresentations = Extension(name = 'frepresentations',
#                          sources = ['qml/frepresentations.f90'],
#                          extra_f90_compile_args = COMPILER_FLAGS,
#                          extra_f77_compile_args = COMPILER_FLAGS,
#                          extra_compile_args = COMPILER_FLAGS,
#                          extra_link_args = MATH_LINKER_FLAGS + LINKER_FLAGS,
#                          language = FORTRAN,
#                          f2py_options=['--quiet'])

# use README.md as long description
def readme():
    with open('README.md') as f:
        return f.read()

def setup_pepytools():

    setup(

        name="fns",
        packages=['fns'],

        # metadata
        version=__version__,
        author=__author__,
        author_email=__email__,
        platforms = 'Any',
        description = __description__,
        long_description = readme(),
        keywords = ['Furthest Neighbour'],
        classifiers = [],
        url = __url__,

        # set up package contents

        ext_package = 'fns',
        #ext_modules = [
        #      ext_frepresentations,
        #],
)

if __name__ == '__main__':

    setup_pepytools()
