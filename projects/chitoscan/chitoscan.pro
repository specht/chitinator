include(../base.pro)

TARGET = chitoscan
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

HEADERS += \
	../../src/ChitoScanner.h \

SOURCES += \
	../../src/ChitoScanner.cpp \
	../../src/chitoscan.cpp \
