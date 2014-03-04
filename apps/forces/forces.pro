include(../../defaults.pri) 

TEMPLATE = app 


SOURCES = forcesmain.cpp \
          sampleforce.cpp
HEADERS = sampleforce.h


LIBS += -L$$TOP_OUT_PWD/lib -lQMC2

INCLUDEPATH += $$TOP_PWD/include

