The code in this repository implementes a bidirectional Burrows-Wheeler index. The implementation is loosely based on the paper "Versatile succinct representations of the bidirectional burrows-wheeler transform.", by Belazzougui et al. at the European Symposium on Algorithms 2013.

The code is set of lightweight headers built on top of the sdsl-lite and divsufsort libraries. Both of the
dependencies are included inside in the sdsl-lite subdirectory, and they can both be built as follows:

```
cd sdsl-lite
sh ./install.sh
```

The library can be used by including the BD_BWT_index.hh header. For example, below is code that 
iterates all right-maximal nodes of the suffix link tree of the string mississippi:

```
#include <iostream>
#include <string>
#include <set>
#include "include/BD_BWT_index.hh"

using namespace std;

set<char> alphabet;

void traverse(BD_BWT_index<>& index, Interval_pair I, string rev_label){
    for(char c : alphabet){
        Interval_pair I2 = index.left_extend(I,c);
        if(I2.forward.size() != 0 && index.is_right_maximal(I2)){
            rev_label += c;

            // Print the label of the current string
            string label(rev_label.rbegin(), rev_label.rend());
            cout << label << endl;

            // Recurse
            traverse(index, I2, rev_label);
            rev_label.pop_back(); // Undo extension
        }
    }
}

int main(){
    string text = "mississippi";
    for(auto c : text) alphabet.insert(c);
    BD_BWT_index<> index((uint8_t*)text.c_str());
    traverse(index, Interval_pair(0, text.size(), 0, text.size()), "");
}
```

This code is at example.cpp. The build the example, you need to include the ./include directory, and link against libsdsl, libdivsufsort and libdivsufsort64. This can be done for example as follows:

```
g++ example.cpp -I include -I sdsl-lite/include -I ./sdsl-lite/build/external/libdivsufsort/include -L sdsl-lite/build/lib -lsdsl -L sdsl-lite/build/external/libdivsufsort/lib/ -ldivsufsort -ldivsufsort64 -o example
```

Documentation is poor at the moment, but the header include/BD_BWT_index.hh contains some documentation.

