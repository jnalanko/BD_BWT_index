#include "bwt.hh"
#include <cstring>
#include <iostream>
#include "tools.hh"
#include <fstream>
#include <string>

extern "C" { 
#include "dbwt.h"
}

using namespace std;

char* bwt_dbwt(char* text){

    unsigned int last;
    long n = strlen(text);
    char* d = (char*)dbwt_bwt((uchar*)text, n, &last, 0);
    d = (char*)realloc(d, (n + 2) * sizeof(char));
    d[last] = '$';
    d[n + 1] = 0;

    return d;
}
