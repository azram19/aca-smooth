CXXFLAGS = -v -std=c++0x -pthread -Ofast -I./vtk/include/vtk-5.10 -Wall -Wno-deprecated -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -mavx -funsafe-math-optimizations -march=corei7-avx -mtune=corei7-avx -mfpmath=sse -mveclibabi=svml -fomit-frame-pointer -funroll-loops

CXX = g++

OBJS = ACA2-2013.o Mesh.o Smooth.o SVD2x2.o

LIBS = -L./vtk/lib/vtk-5.10 -lvtkIO -lvtkFiltering -lvtkCommon -lvtkzlib -lvtkexpat -lvtksys -ldl -lpthread

TARGET = ACA2-2013

$(TARGET):	$(OBJS)
	@$(CXX) -v -o $(TARGET) $(OBJS) $(LIBS)

all: @$(TARGET)

clean:
	@rm -f $(OBJS) $(TARGET)
