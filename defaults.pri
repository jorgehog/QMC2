CONFIG -= app_bundle
CONFIG -= qt

ABEL {
    ARMAPATH = /usit/abel/u1/jorgehog/libs/armadillo-3.920.1

    INCLUDEPATH += $$ARMAPATH/usr/include \
                /cluster/software/VERSIONS/intel-2013.2/mkl/include

    LIBS +=  -L/cluster/software/VERSIONS/intel-2013.2/mkl/lib/intel64 \
             -L$$ARMAPATH/usr/lib64 \
             -lpthread \
             -liomp5

    DEFINES += MKL_LP64 LINED_OUTPUT
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
QMAKE_LFLAGS += $$system(mpic++ --showme:link)
QMAKE_CXXFLAGS += $$system(mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK


#C++11 features
QMAKE_CXXFLAGS += -std=c++0x

QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS -g

QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS -O3 -DARMA_NO_DEBUG

TOP_PWD=$$PWD
TOP_OUT_PWD=$$shadowed($$PWD)
