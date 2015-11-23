#include <iostream>
#include <string>
#include "BD_BWT_index.hh"
#include "Iterators.hh"
#include <cassert>
#include <set>

using namespace std;

set<string> get_right_maximal_substrings(const string& s){
    int n = s.size();
    set<string> result;
    set<char> alphabet(s.begin(),s.end());
    set<string> substrings;
    substrings.insert("");
    for(int k = 1; k <= n; k++){
        for(int i = 0; i <= n - k; i++){
            substrings.insert(s.substr(i,k));
        }
    }
    for(const string& x : substrings){
        int nExtensions = 0; // Number of ways to extend x right. End of string counts a distinct extension.
        for(char c : alphabet){
            if(substrings.count(x + c) == 1)
                nExtensions++;
        }
        if(s.substr(n-x.size()) == x) nExtensions++; // x is at the right end of s
        // Note subtr: "If this is equal to the string length, the function returns an empty string."
        
        if(nExtensions >= 2)
            result.insert(x);
    }
    return result;
}

vector<string> all_binary_strings_up_to(int k){
    vector<string> ans;
    ans.push_back("");
    for(int length = 1; length <= k; length++){
        for(int mask = 0; mask < (1 << length); mask++){
            string s = "";
            for(int i = 0; i < length; i++){
                if(mask & (1 << i)) s += 'a';
                else s += 'b';
            }
            ans.push_back(s);
        }
    }
    return ans;
}


bool test_suffix_link_tree_iteration(string& s){
    BD_BWT_index<sdsl::bit_vector> index(s);
    BD_BWT_index_iterator<sdsl::bit_vector> it(&index);
    set<string> labels;
    while(it.next()){
        string x(it.label.rbegin(), it.label.rend());
        labels.insert(x);
    }
    
    return (labels == get_right_maximal_substrings(s));       
}

int main(int argc, char** argv){
    
    vector<string> test_set = all_binary_strings_up_to(10);
    for(auto& s : test_set){
        cout << s << endl;
        assert(test_suffix_link_tree_iteration(s));
    }
    
    cerr << "All tests OK" << endl;
    
    //BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index) : index(index), char_counts(256) {
    //void initialize_index_internal_memory(BD_BWT_index<t_bitvector>& index, std::string& input);
}