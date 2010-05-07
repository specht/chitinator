include(../base.pro)

TARGET = expmasses
CONFIG(debug, debug|release) {
    TARGET = $$join(TARGET,,,_debug)
}

HEADERS += \
    ../../src/Digestion.h \
    ../../src/Enzyme.h \
    ../../src/Polymer.h \

SOURCES += \
    ../../src/expmasses.cpp \
    ../../src/Digestion.cpp \
    ../../src/Enzyme.cpp \
    ../../src/Polymer.cpp \
