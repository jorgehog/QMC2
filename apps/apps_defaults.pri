include(../defaults.pri)

TEMPLATE = app

QMAKE_LIBDIR += $$TOP_PWD/lib

LIBS += -L$$TOP_PWD/lib -lQMC2

INCLUDEPATH += $$TOP_PWD/include

