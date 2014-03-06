include(../../defaults.pri) 
include(../apps_defaults.pri)

INCLUDEPATH += $$PWD/HF/src/libs

LIBS += -L$$PWD/HF/src/libs/ -lhartree-fock

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
