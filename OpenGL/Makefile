# Usage:
# For compiling for linux:
#    make OS="linux"
# For compiling for windows:
#    make OS="win32"
# or
#    make OS="win64"
# or
#    make OS="osx"
# For cross-compiling for windows:
#    make OS="cross"

OS = linux
ARCH = $(shell getconf LONG_BIT)
CC=c++
cc=gcc
WARCH=w64

ifeq ($(OS),linux)
   OGLLIBS = -lGLU -lGL -lglut -lm -lIL
   OGLINC = -I/usr/include/eigen3
   EXT = "-linux$(ARCH)"
else
   ifeq ($(OS),win32)
       OGLLIBS = -lglut32 -lglu32 -lopengl32
       EXT=".exe"
   else
        ifeq ($(OS),osx)
          OGLINC = -I/usr/X11/include/ -I/opt/local/include/ -I/opt/local/include/eigen3/
          OGLPATH = -L/usr/X11/lib/ -L/opt/local/lib
          OGLLIBS = -lGLU -lGL -lglut -lIL
          EXT=".osx"
        else
          CC=i686-$(WARCH)-mingw32-g++
          cc=i686-$(WARCH)-mingw32-gcc
          OGLLIBS = -lglut -lglu32 -lopengl32
          ifeq ($(WARCH),w64)
             EXT="64.exe"
          else
             EXT=".exe"
          endif
          OPTIONS=-static-libgcc -static-libstdc++
       endif
   endif
endif

FLAGS = -I. $(OGLINC) -D__KEEP_ASPECT__ -Wno-deprecated-declarations

all: drawBird$(EXT) housePerspective$(EXT) testeGL$(EXT) planet$(EXT) torus$(EXT) Arcball/Arcball.o
	@echo	Compilation completed.

drawBird$(EXT): drawBird.c
	$(cc) $(FLAGS) -o drawBird$(EXT) drawBird.c $(OGLPATH) $(OGLLIBS) $(OPTIONS)

anim$(EXT): anim.c
	$(cc) $(FLAGS) -o anim$(EXT) anim.c $(OGLPATH) $(OGLLIBS) $(OPTIONS)

housePerspective$(EXT): housePerspective.c
	$(cc) $(FLAGS) -o housePerspective$(EXT) housePerspective.c $(OGLPATH) $(OGLLIBS) $(OPTIONS)

Arcball/Arcball.o: Arcball/Arcball.cpp
	$(CC) $(FLAGS) -c Arcball/Arcball.cpp -o Arcball/Arcball.o

torus$(EXT): torus.cpp Arcball/Arcball.o
	$(CC) $(FLAGS) -o torus$(EXT) torus.cpp $(OGLPATH) $(OGLLIBS) $(OPTIONS) Arcball/Arcball.o

testeGL$(EXT): testeGL.c
	$(cc) $(FLAGS) -o testeGL$(EXT) testeGL.c $(OGLPATH) $(OGLLIBS) $(OPTIONS)

planet$(EXT): planet.c
	$(cc) $(FLAGS) -o planet$(EXT) planet.c $(OGLPATH) $(OGLLIBS) $(OPTIONS)

clip$(EXT): clip.c
	$(cc) $(FLAGS) -o clip$(EXT) clip.c $(OGLPATH) $(OGLLIBS) $(OPTIONS)

