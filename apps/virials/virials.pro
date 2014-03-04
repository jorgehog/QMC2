include(../../defaults.pri) 

TEMPLATE = app 


SOURCES = virialsmain.cpp


LIBS += -L$$TOP_OUT_PWD/lib -lQMC2

INCLUDEPATH += $$TOP_PWD/include

