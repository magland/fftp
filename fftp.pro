QT += core network
QT -= gui

CONFIG += c++11

CONFIG -= app_bundle #Please apple, don't make a bundle

DESTDIR = bin
OBJECTS_DIR = build
MOC_DIR = build
TARGET = fftp
TEMPLATE = app

HEADERS += \
    fftp.h

SOURCES += fftpmain.cpp \
    fftp.cpp

