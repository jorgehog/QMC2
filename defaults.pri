CONFIG -= app_bundle
CONFIG -= qt

ABEL {
    ARMAPATH = /usit/abel/u1/jorgehog/libs/armadillo-3.920.1

    INCLUDEPATH += $$ARMAPATH/usr/include \
                /cluster/software/VERSIONS/intel-2013.2/mkl/include
                .

    LIBS +=  -L/cluster/software/VERSIONS/intel-2013.2/mkl/lib/intel64 \
             -L$$ARMAPATH/usr/lib64 \
             -lpthread \
             -liomp5

    DEFINES += MKL_LP64
    QMAKE_CXXFLAGS += -openmp

    QMAKE_LFLAGS -= -lm -O1

}

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK  #--enable-shared --enable-static


#C++11 features
QMAKE_CXXFLAGS += -std=c++0x -DARMA_USE_CXX11

QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS -g

QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS -O3 -DARMA_NO_DEBUG


DIRS = scratch/QMC_SCRATCH/walker_positions

for(DIR, DIRS) {
     mkcommands += $$TOP_OUT_PWD/$$DIR
}

#createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) #createDirs
export(first.depends)
#export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first #createDirs

