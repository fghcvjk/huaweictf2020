CXX=g++
CXXFLAGS=-Iinclude/ -no-pie -g
LDFLAGS=lib/libpng16.a lib/libz.a -no-pie

OBJECTS=src/main.o src/image.o src/stencil.o

stencil: $(OBJECTS)
	$(CXX) -o bin/stencil $(OBJECTS) $(LDFLAGS)

all:	stencil

run:	all
	bin/stencil test-image.png

clean:
	rm -f $(OBJECTS) bin/stencil output.png src/*~ *~
