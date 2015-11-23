dbwtSrc=dbwt/dbwt.c dbwt/dbwt_queue.c dbwt/dbwt_utils.c dbwt/sais.c
lib_paths=-L lib/sdsl/build/lib -L lib/sdsl/build/external/libdivsufsort/lib/ -L lib
includes=-I lib/sdsl/include -I lib/sdsl/build/external/libdivsufsort/include -I src -I dbwt
cppflags = -std=c++11 -O3 -MMD
sources = main.cpp bwt.cpp tools.cpp
link = -pthread -fopenmp -ldbwt -ldivsufsort64 -lsdsl

# Search directory for source files
VPATH=src

# Substitute .cpp with .o
objects = $(patsubst %.cpp,build/%.o,$(sources)) 

.PHONY: directories clean

# Include header dependencies
-include $(objects:%.o=%.d)

all: tests

# Map object files to source files. " | build " is an order-only prerequisite meaning
# that the build directory must exist before building objects, but
# the objects should must not be rebuilt if the timestamp of the build
# directory changes, which happens every time a file is modified in the directory.
build/%.o : %.cpp
	$(CXX) -c $< $(cppFlags) $(includes) -o $@

directories:
	@mkdir build -p
	@mkdir include -p
	@mkdir lib -p

tests: $(objects)
	$(CXX)  $(cppFlags) $(objects) tests.cpp $(link) -o tests
	
dbwt:
	gcc -std=c99 -O3 -g $(dbwtSrc) -I ./dbwt -c
	ar rcs lib/libdbwt.a *.o
	rm *.o

clean:
	rm build/*
