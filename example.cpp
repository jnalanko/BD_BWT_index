#include <iostream>
#include "include/BD_BWT_index.hh"
#include <string>
#include <set>

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
