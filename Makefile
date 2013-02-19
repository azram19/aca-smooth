CXXFLAGS = -std=c++0x -g -pthread -O3 -mrecip=all -ffast-math -Wall -Wno-deprecated -I./vtk/include/vtk-5.10 -funsafe-math-optimizations -march=core2 -mtune=core2 -mfpmath=sse -mveclibabi=svml -fomit-frame-pointer -funroll-loops -ftree-vectorize

OBJS = ACA2-2013.o Mesh.o Smooth.o SVD2x2.o Q.o

LIBS = -L./vtk/lib/vtk-5.10 -lvtkIO -lvtkFiltering -lvtkCommon -lvtkzlib -lvtkexpat -lvtksys -ldl -lpthread

TARGET = ACA2-2013

$(TARGET):	$(OBJS)
	@$(CXX) -S $(TARGET) $(OBJS) $(LIBS)
	@$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	@$(TARGET)

clean:
	@rm -f $(OBJS) $(TARGET)
