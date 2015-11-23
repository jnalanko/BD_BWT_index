sources = main.cpp bwt.cpp tools.cpp
dbwt_sources=dbwt/dbwt.c dbwt/dbwt_queue.c dbwt/dbwt_utils.c dbwt/sais.c
lib_paths=-L sdsl-lite/build/lib -L sdsl-lite/build/external/libdivsufsort/lib/ -L lib
includes=-I sdsl-lite/include -I sdsl-lite/build/external/libdivsufsort/include -I include -I dbwt
cxxflags = -std=c++11 -O3 -MMD
ccflags = -std=c99 -O3 -g -MMD
link = -pthread -fopenmp -ldivsufsort64 -lsdsl

# Search directory for source files
VPATH=src:dbwt

# Substitute .cpp with .o
objects = $(patsubst %.cpp, build/%.o, $(sources)) 

# Substitute .c with .o
dbwt_objects = $(patsubst %.c, %.o, $(dbwt_sources)) 

.PHONY: clean

# Include header dependencies
-include $(objects:%.o=%.d)

all: tests

# Map object files to source files. " | build " is an order-only prerequisite meaning
# that the build directory must exist before building objects, but
# the objects should must not be rebuilt if the timestamp of the build
# directory changes, which happens every time a file is modified in the directory.
build/%.o : %.cpp | build
	$(CXX) -c $< $(cxxflags) $(includes) -o $@
	
dbwt/%.o: dbwt/%.c
	$(CC) $(ccflags) -I dbwt -c $< -o $@

build:
	@mkdir build -p
	
tests: $(objects) $(dbwt_objects)
	$(CXX) $(objects) $(dbwt_objects) $(lib_paths) $(link) -o tests
	
clean:
	rm build/*.o
	rm build/*.d
	rm dbwt/*.o
	rm dbwt/*.d
