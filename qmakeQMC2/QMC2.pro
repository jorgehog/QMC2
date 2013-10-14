
TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += CXX11

LIBS += -llapack -lblas -larmadillo

SOURCES += \
    ../src/BasisFunctions/HarmonicOscillator/HarmonicOscillator.cpp \
    ../src/BasisFunctions/HarmonicOscillator3D/HarmonicOscillator3D.cpp \
    ../src/BasisFunctions/hydrogenic/hydrogenic.cpp \
    ../src/Diffusion/Diffusion.cpp \
    ../src/Diffusion/Fokker_Planck/Fokker_Planck.cpp \
    ../src/Diffusion/RNGs/gaussian_deviate.cpp \
    ../src/Diffusion/RNGs/ran2.cpp \
    ../src/Diffusion/RNGs/zignor.cpp \
    ../src/Diffusion/RNGs/zigrandom.cpp \
    ../src/Diffusion/Simple/Simple.cpp \
    ../src/ErrorEstimator/Blocking/Blocking.cpp \
    ../src/ErrorEstimator/ErrorEstimator.cpp \
    ../src/ErrorEstimator/SimpleVar/SimpleVar.cpp \
    ../src/Jastrow/Jastrow.cpp \
    ../src/Jastrow/No_Jastrow/No_Jastrow.cpp \
    ../src/Jastrow/Pade_Jastrow/Pade_Jastrow.cpp \
    ../src/Minimizer/ASGD/ASGD.cpp \
    ../src/Minimizer/Minimizer.cpp \
    ../src/Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.cpp \
    ../src/Orbitals/DiTransform/DiTransform.cpp \
    ../src/Orbitals/ExpandedBasis/ExpandedBasis.cpp \
    ../src/Orbitals/Orbitals.cpp \
    ../src/Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.cpp \
    ../src/OutputHandler/Distribution/Distribution.cpp \
    ../src/OutputHandler/OutputHandler.cpp \
    ../src/OutputHandler/stdoutASGD/stdoutASGD.cpp \
    ../src/OutputHandler/stdoutDMC/stdoutDMC.cpp \
    ../src/Potential/AtomCore/AtomCore.cpp \
    ../src/Potential/Coulomb/Coulomb.cpp \
    ../src/Potential/DiAtomCore/DiAtomCore.cpp \
    ../src/Potential/DoubleWell/DoubleWell.cpp \
    ../src/Potential/Harmonic_osc/Harmonic_osc.cpp \
    ../src/Potential/Potential.cpp \
    ../src/QMC/DMC/DMC.cpp \
    ../src/QMC/QMC.cpp \
    ../src/QMC/VMC/VMC.cpp \
    ../src/BasisFunctions/BasisFunctions.cpp \
    ../src/Sampling/Brute_Force/Brute_Force.cpp \
    ../src/Sampling/Importance/Importance.cpp \
    ../src/Sampling/Sampling.cpp \
    ../src/System/Fermions/Fermions.cpp \
    ../src/System/System.cpp \
    ../src/Walker/Walker.cpp \
    ../src/Sampler/sampleMethods/SampleForce.cpp \
    ../src/forcesMain.cpp \
    ../src/Orbitals/NBodyTransform/nbodytransform.cpp \
    ../src/Potential/MolecularCoulomb/molecularcoulomb.cpp
    

HEADERS += \
    ../src/structs.h \
    ../src/BasisFunctions/BasisFunctions.h \
    ../src/BasisFunctions/HarmonicOscillator/HarmonicOscillator.h \
    ../src/BasisFunctions/HarmonicOscillator3D/HarmonicOscillator3D.h \
    ../src/BasisFunctions/hydrogenic/hydrogenic.h \
    ../src/Diffusion/Diffusion.h \
    ../src/Diffusion/Fokker_Planck/Fokker_Planck.h \
    ../src/Diffusion/RNGs/gaussian_deviate.h \
    ../src/Diffusion/RNGs/ran2.h \
    ../src/Diffusion/RNGs/zignor.h \
    ../src/Diffusion/RNGs/zigrandom.h \
    ../src/Diffusion/Simple/Simple.h \
    ../src/ErrorEstimator/Blocking/Blocking.h \
    ../src/ErrorEstimator/ErrorEstimator.h \
    ../src/ErrorEstimator/SimpleVar/SimpleVar.h \
    ../src/Jastrow/Jastrow.h \
    ../src/Jastrow/No_Jastrow/No_Jastrow.h \
    ../src/Jastrow/Pade_Jastrow/Pade_Jastrow.h \
    ../src/Minimizer/ASGD/ASGD.h \
    ../src/Minimizer/Minimizer.h \
    ../src/Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h \
    ../src/Orbitals/DiTransform/DiTransform.h \
    ../src/Orbitals/ExpandedBasis/ExpandedBasis.h \
    ../src/Orbitals/Orbitals.h \
    ../src/Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.h \
    ../src/OutputHandler/Distribution/Distribution.h \
    ../src/OutputHandler/OutputHandler.h \
    ../src/OutputHandler/stdoutASGD/stdoutASGD.h \
    ../src/OutputHandler/stdoutDMC/stdoutDMC.h \
    ../src/Potential/AtomCore/AtomCore.h \
    ../src/Potential/Coulomb/Coulomb.h \
    ../src/Potential/DiAtomCore/DiAtomCore.h \
    ../src/Potential/DoubleWell/DoubleWell.h \
    ../src/Potential/Harmonic_osc/Harmonic_osc.h \
    ../src/Potential/Potential.h \
    ../src/QMC/DMC/DMC.h \
    ../src/QMC/QMC.h \
    ../src/QMC/VMC/VMC.h \
    ../src/Sampler/Sampler.h \
    ../src/Sampling/Brute_Force/Brute_Force.h \
    ../src/Sampling/Importance/Importance.h \
    ../src/Sampling/Sampling.h \
    ../src/System/Fermions/Fermions.h \
    ../src/System/System.h \
    ../src/Walker/Walker.h \
    ../src/defines.h \
    ../src/Sampler/sampleMethods/SampleForce.h \
    ../src/Orbitals/NBodyTransform/nbodytransform.h \
    ../src/Potential/MolecularCoulomb/molecularcoulomb.h

ABEL {
    ARMAPATH = /usit/abel/u1/jorgehog/libs/armadillo-3.920.1

    INCLUDEPATH += $$ARMAPATH/usr/include \
                /cluster/software/VERSIONS/intel-2013.2/mkl/include \
                .

    LIBS -= -llapack -lblas
    LIBS +=  -L/cluster/software/VERSIONS/intel-2013.2/mkl/lib/intel64 \
             -L$$ARMAPATH/usr/lib64 \
             -lpthread \
             -liomp5
    DEFINES += MKL_LP64
    QMAKE_CXXFLAGS += -openmp
}

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc
 
QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK


#C++11 features
CXX11 {
    QMAKE_CXXFLAGS += -std=c++0x -DARMA_USE_CXX11
}

QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS -g

QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS -O3 -DARMA_NO_DEBUG

!ABEL {
    #remove this line if you don't use ccache
    QMAKE_CXX = ccache $$QMAKE_CXX
}

DEFINES += LINED_OUTPUT
