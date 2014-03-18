include(../defaults.pri)

TEMPLATE = lib
TARGET = ../lib/QMC2

LIBS += -larmadillo

QMAKE_LFLAGS += -g

SOURCES += \
    BasisFunctions/HarmonicOscillator/HarmonicOscillator.cpp \
    BasisFunctions/HarmonicOscillator3D/HarmonicOscillator3D.cpp \
    BasisFunctions/hydrogenic/hydrogenic.cpp \
    Diffusion/Diffusion.cpp \
    Diffusion/Fokker_Planck/Fokker_Planck.cpp \
    Diffusion/RNGs/zignor.cpp \
    Diffusion/RNGs/zigrandom.cpp \
    Diffusion/Simple/Simple.cpp \
    ErrorEstimator/Blocking/Blocking.cpp \
    ErrorEstimator/ErrorEstimator.cpp \
    ErrorEstimator/SimpleVar/SimpleVar.cpp \
    Jastrow/Jastrow.cpp \
    Jastrow/No_Jastrow/No_Jastrow.cpp \
    Jastrow/Pade_Jastrow/Pade_Jastrow.cpp \
    Minimizer/ASGD/ASGD.cpp \
    Minimizer/Minimizer.cpp \
    Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.cpp \
    Orbitals/Orbitals.cpp \
    Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.cpp \
    OutputHandler/Distribution/Distribution.cpp \
    OutputHandler/OutputHandler.cpp \
    OutputHandler/stdoutASGD/stdoutASGD.cpp \
    OutputHandler/stdoutDMC/stdoutDMC.cpp \
    Potential/AtomCore/AtomCore.cpp \
    Potential/Coulomb/Coulomb.cpp \
    Potential/DoubleWell/DoubleWell.cpp \
    Potential/Harmonic_osc/Harmonic_osc.cpp \
    Potential/Potential.cpp \
    QMC/DMC/DMC.cpp \
    QMC/QMC.cpp \
    QMC/VMC/VMC.cpp \
    Sampling/Brute_Force/Brute_Force.cpp \
    Sampling/Importance/Importance.cpp \
    Sampling/Sampling.cpp \
    System/Fermions/Fermions.cpp \
    System/System.cpp \
    Walker/Walker.cpp \
    Orbitals/NBodyTransform/nbodytransform.cpp \
    Potential/MolecularCoulomb/molecularcoulomb.cpp \
    Sampler/Sampler.cpp

HEADERS += \
    structs.h \
    defines.h \
    BasisFunctions/BasisFunctions.h \
    BasisFunctions/HarmonicOscillator/HarmonicOscillator.h \
    BasisFunctions/HarmonicOscillator3D/HarmonicOscillator3D.h \
    BasisFunctions/hydrogenic/hydrogenic.h \
    Diffusion/Diffusion.h \
    Diffusion/Fokker_Planck/Fokker_Planck.h \
    Diffusion/RNGs/gaussian_deviate.h \
    Diffusion/RNGs/ran2.h \
    Diffusion/RNGs/zignor.h \
    Diffusion/RNGs/zigrandom.h \
    Diffusion/Simple/Simple.h \
    ErrorEstimator/Blocking/Blocking.h \
    ErrorEstimator/ErrorEstimator.h \
    ErrorEstimator/SimpleVar/SimpleVar.h \
    Jastrow/Jastrow.h \
    Jastrow/No_Jastrow/No_Jastrow.h \
    Jastrow/Pade_Jastrow/Pade_Jastrow.h \
    Minimizer/ASGD/ASGD.h \
    Minimizer/Minimizer.h \
    Orbitals/AlphaHarmonicOscillator/AlphaHarmonicOscillator.h \
    Orbitals/ExpandedBasis/ExpandedBasis.h \
    Orbitals/Orbitals.h \
    Orbitals/OrbitalsFactory.h \
    Orbitals/hydrogenicOrbitals/hydrogenicOrbitals.h \
    OutputHandler/Distribution/Distribution.h \
    OutputHandler/OutputHandler.h \
    OutputHandler/stdoutASGD/stdoutASGD.h \
    OutputHandler/stdoutDMC/stdoutDMC.h \
    Potential/AtomCore/AtomCore.h \
    Potential/Coulomb/Coulomb.h \
    Potential/DoubleWell/DoubleWell.h \
    Potential/Harmonic_osc/Harmonic_osc.h \
    Potential/Potential.h \
    QMC/DMC/DMC.h \
    QMC/QMC.h \
    QMC/VMC/VMC.h \
    Sampler/Sampler.h \
    Sampling/Brute_Force/Brute_Force.h \
    Sampling/Importance/Importance.h \
    Sampling/Sampling.h \
    System/Fermions/Fermions.h \
    System/System.h \
    Walker/Walker.h \
    Orbitals/NBodyTransform/nbodytransform.h \
    Potential/MolecularCoulomb/molecularcoulomb.h \
    Sampler/standardsamplers.h \
    misc.h

OTHER_FILES += ../include/QMC2.h ../install/include/QMC2.h

NO_MPI {
    DEFINES += QMC2_NO_MPI
}

!equals($$OUT_PWD/.., $$TOP_PWD) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}

#target.path  = /usr/lib/

#all_headers.files  = $$HEADERS
#all_headers.path   = /usr/include/QMC2_bits

#main_include.files = $$TOP_PWD/install/include/*
#main_include.path  = /usr/include/

#INSTALLS += target main_include all_headers
