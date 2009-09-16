TEMPLATE = app

CONFIG += debug_and_release console

DEPENDPATH += .
INCLUDEPATH += .

LIBS += -lptb -lbz2 -lquazip

macx {
	CONFIG -= app_bundle
	CONFIG += x86
}

win32 {
	LIBS += -lzdll
} else {
	LIBS += -lz
}

CONFIG(debug, debug|release) {
	OBJECTS_DIR = ../../obj/debug/
	MOC_DIR = ../../obj/debug/
	RCC_DIR = ../../obj/debug/
}
else {
	OBJECTS_DIR = ../../obj/release/
	MOC_DIR = ../../obj/release/
	RCC_DIR = ../../obj/release/
}

DESTDIR = ../../

QT = core xml

HEADERS += \
	../../src/RefPtr.h \

SOURCES += \
