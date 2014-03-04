include(../../defaults.pri) 

TEMPLATE = app 

TARGET = QMC2_default


SOURCES = defaultmain.cpp


LIBS += -L$$TOP_OUT_PWD/lib -lQMC2

INCLUDEPATH += $$TOP_PWD/include

