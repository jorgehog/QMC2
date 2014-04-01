include(../apps_defaults.pri)

INCLUDEPATH += $$PWD/HF/include

LIBS += -L$$PWD/HF/lib/ -lhartree-fock
LIBS += -lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
LIBS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)

SOURCES = hfqmcmain.cpp \
    hforbitals.cpp

#Copy infiles
copydata.commands = $(COPY_DIR) $$PWD/HF/infiles $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata

HEADERS += \
    hforbitals.h
