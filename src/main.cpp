#include <iostream>
#include <string>
#include "BD_BWT_index.hh"
#include "Iterators.hh"
#include <cassert>

using namespace std;

int main(int argc, char** argv){
    
    string s = "abracadabra";
    cout << s << endl;
    BD_BWT_index<sdsl::bit_vector> index(s);
    int dollar = -1;
    for(int i = 0; i < index.size(); i++){
        cout << index.forward_bwt_at(i);
        if(index.forward_bwt_at(i) == '$') dollar = i;
    }
    cout << endl;
    
    assert(dollar != -1);
    int pos = dollar;
    for(int i = 0; i < 2*(s.size() + 1); i++){
        pos = index.backward_step(pos);
        cout << index.forward_bwt_at(i % index.size());
    }
    cout << endl;
    
    BD_BWT_index_iterator<sdsl::bit_vector> it(&index);
    while(it.next()){
        cout << it.label << endl;
    }

    
    //BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index) : index(index), char_counts(256) {
    //void initialize_index_internal_memory(BD_BWT_index<t_bitvector>& index, std::string& input);
}