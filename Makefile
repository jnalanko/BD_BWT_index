CXX = g++
CXXFLAGS = -O3 -g -std=c++14

CC = gcc
CFLAGS = -O3 -g

#DBWT_OBJECTS = build/dbwt.o build/dbwt_queue.o build/dbwt_utils.o build/sais.o
#OUR_OBJECTS = build/io_tools.o

#build/dbwt.o:
	#$(CC) dbwt/dbwt.c $(CFLAGS) -c -o build/dbwt.o

#build/dbwt_queue.o:
#	$(CC) dbwt/dbwt_queue.c $(CFLAGS) -c -o build/dbwt_queue.o

#build/dbwt_utils.o:
#	$(CC) dbwt/dbwt_utils.c $(CFLAGS) -c -o build/dbwt_utils.o

#build/sais.o:
#	$(CC) dbwt/sais.c $(CFLAGS) -c -o build/sais.o

#build/io_tools.o:
#$(CXX) src/io_tools.cpp $(CXXFLAGS) -c -o build/io_tools.o -I include

all:
	cp sdsl-lite/build/lib/libsdsl.a build/lib
	cp sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a build/lib
	#ar -rs build/lib/libbdbwt.a $(OUR_OBJECTS) 

