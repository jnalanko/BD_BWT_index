sources = bwt.cpp io_tools.cpp
dbwt_sources=dbwt/dbwt.c dbwt/dbwt_queue.c dbwt/dbwt_utils.c dbwt/sais.c
lib_path= -L lib
includes= -I include 
cxxflags = -std=c++11 -O3 -g -MMD
ccflags = -std=c99 -O3 -g -MMD
link = -lbdbwt -ldbwt -ldivsufsort64 -lsdsl

# Search directory for source files
VPATH=src:dbwt

# Substitute .cpp with .o
objects = $(patsubst %.cpp, build/%.o, $(sources)) 

# Substitute .c with .o
dbwt_objects = $(patsubst %.c, %.o, $(dbwt_sources)) 

.PHONY: clean install

install: include/divsufsort.h include/divsufsort64.h include/dbwt.h include/dbwt_queue.h include/dbwt_utils.h lib/libbdbwt.a lib/libdbwt.a lib/libsdsl.a lib/libdivsufsort.a include/sdsl-lite 

all: install tests

include/dbwt.h:
	cp dbwt/dbwt.h include/dbwt.h

include/dbwt_queue.h:
	cp dbwt/dbwt_queue.h include/dbwt_queue.h

include/dbwt_utils.h:
	cp dbwt/dbwt_utils.h include/dbwt_utils.h

lib/libbdbwt.a: $(objects) | lib
	ar rcs lib/libbdbwt.a $(objects)

lib/libdbwt.a: $(dbwt_objects) | lib
	ar rcs lib/libdbwt.a $(dbwt_objects)

lib/libsdsl.a: | lib
	cp sdsl-lite/build/lib/libsdsl.a lib
    
lib/libdivsufsort.a: | lib
	cp sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a lib

include/sdsl-lite:
	cp -r sdsl-lite/include/sdsl include

include/divsufsort.h:
	cp sdsl-lite/build/external/libdivsufsort/include/divsufsort.h include

include/divsufsort64.h:
	cp sdsl-lite/build/external/libdivsufsort/include/divsufsort64.h include

# Include header dependencies (todo: main.o missing from here)
-include $(objects:%.o=%.d)

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
	
lib:
	@mkdir lib -p
	
tests: install
	$(CXX) src/main.cpp $(lib_path) $(link) $(includes) $(cxxflags) -o tests
	
clean:
	rm build/*.o
	rm build/*.d
	rm dbwt/*.o
	rm dbwt/*.d
