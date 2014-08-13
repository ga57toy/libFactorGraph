#-------------------------------------------------
#
# Project created by QtCreator 2014-08-13T16:29:55
#
#-------------------------------------------------

QT       -= gui

TARGET = libFactorGraph
TEMPLATE = lib

DEFINES += LIBFACTORGRAPH_LIBRARY

SOURCES += libfactorgraph.cpp

HEADERS += libfactorgraph.h\
        libfactorgraph_global.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}

OTHER_FILES += \
    README.txt
