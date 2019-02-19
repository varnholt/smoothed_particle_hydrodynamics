TARGET = sph
TEMPLATE = app

QT += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

DEFINES += _USE_MATH_DEFINES

# openmp support
msvc:QMAKE_CXXFLAGS_RELEASE += /openmp
QMAKE_CXXFLAGS += -openmp

win32
{
   LIBS += -lglu32
   LIBS += -lopengl32
   LIBS += -lvcomp
}

unix {
   LIBS += -lGLU
   LIBS += -lgomp
}

SOURCES += \
   src/main.cpp\
   src/sph.cpp \
   src/particle.cpp \
   src/visualization.cpp \
   src/widget.cpp \
   src/sphconfig.cpp \
   src/vec3.cpp

HEADERS += \
   src/sph.h \
   src/particle.h \
   src/visualization.h \
   src/widget.h \
   src/sphconfig.h \
   src/vec3.h

FORMS += src/widget.ui

INCLUDEPATH += src
