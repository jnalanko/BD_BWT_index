#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <set>
#include <deque>
#include <stack>
#include <map>
#include "bwt.hh"
#include "malloc.h"
#include "tools.hh"

using namespace sdsl;
using namespace std;

template<class t_bitvector>
int64_t compute_cumulative_char_rank_in_interval(const wt_huff<t_bitvector>& wt, vector<uint8_t> alphabet, char c, Interval I){
    int64_t ans = 0;
    if(I.size() == 0) return 0;
    
    // Sum of ranks of all characters that are lexicographically smaller than c
    for(char d : alphabet){
        if(d == c) break;
        ans += wt.rank(I.right + 1,d) - wt.rank(I.left,d);
    }
    return ans;
}


template<class t_bitvector>
vector<uint8_t> get_interval_symbols(const wt_huff<t_bitvector>& wt, Interval I){
    if(I.size() == 0){
        vector<uint8_t> empty;
        return empty;
    }
    
    int_vector_size_type nExtensions;
    vector<uint8_t> symbols(wt.sigma);
    vector<uint64_t> ranks_i(wt.sigma);
    vector<uint64_t> ranks_j(wt.sigma);
    wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
    while(symbols.size() > nExtensions) symbols.pop_back();
    return symbols;
}

// Another version of get_interval_symbols when the allocation of a new vector is too slow.
// Counts the number of distinct symbols into nExtensions, and puts the symbols into the indices
// [0,nExtensions[ of symbols. Also stores ranks of the symbols at the endpoints of the interval I
// to ranks_i and ranks_j. Important: All the parameter vectors must have length at least equal to the size
// of the alphabet of the given wavelet tree. Symbols is not sorted
template<class t_bitvector>
void get_interval_symbols(const wt_huff<t_bitvector>& wt, Interval I, int_vector_size_type& nExtensions, vector<uint8_t>& symbols,
 vector<uint64_t>& ranks_i, vector<uint64_t>& ranks_j){
    if(I.size() == 0){
        nExtensions = 0;
        return;
    }
    
    wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
}


// Takes a backward step in the forward bwt
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::backward_step(int64_t pos){
    char c = forward_bwt[pos];
    return cumulative_char_count[c] + forward_bwt.rank(pos, c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::left_extend(Interval_pair intervals, char c){
    int64_t cumul_rank_c = compute_cumulative_char_rank_in_interval(forward_bwt, alphabet, c, intervals.forward);
    return left_extend(intervals,c,cumul_rank_c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::right_extend(Interval_pair intervals, char c){
    int64_t cumul_rank_c = compute_cumulative_char_rank_in_interval(reverse_bwt, alphabet, c, intervals.reverse);
    return right_extend(intervals,c,cumul_rank_c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::left_extend(Interval_pair intervals, char c, int64_t cumul_rank_c){
    if(intervals.forward.size() == 0)
        return Interval_pair(-1,-2,-1,-2);
    
    Interval forward = intervals.forward;
    Interval reverse = intervals.reverse;
    
    // Compute the new forward interval
    int64_t num_c_in_interval = forward_bwt.rank(forward.right + 1,c) - forward_bwt.rank(forward.left,c);
    int64_t start_f_new = cumulative_char_count[c] + forward_bwt.rank(forward.left, c); // Start in forward
    int64_t end_f_new = start_f_new + num_c_in_interval - 1; // End in forward

    if(start_f_new > end_f_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
    
    // Compute the new reverse interval
    int64_t start_r_new = reverse.left + cumul_rank_c;
    int64_t end_r_new = start_r_new + (end_f_new - start_f_new); // The forward and reverse intervals must have same length
    
    return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::right_extend(Interval_pair intervals, char c, int64_t cumul_rank_c){
    if(intervals.forward.size() == 0)
        return Interval_pair(-1,-2,-1,-2);
    
    Interval forward = intervals.forward;
    Interval reverse = intervals.reverse;
    
    // Compute the new reverse interval
    int64_t num_c_in_interval = reverse_bwt.rank(reverse.right + 1,c) - reverse_bwt.rank(reverse.left,c);
    int64_t start_r_new = cumulative_char_count[c] + reverse_bwt.rank(reverse.left, c); // Start in reverse
    int64_t end_r_new = start_r_new + num_c_in_interval - 1; // End in reverse

    if(start_r_new > end_r_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
    
    // Compute the new forward interval
    int64_t cumul_rank = cumul_rank_c;
    int64_t start_f_new = forward.left + cumul_rank;
    int64_t end_f_new = start_f_new + (end_r_new - start_r_new); // The forward and reverse intervals must have same length
    
    return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
    
}


// Compute the cumulative sum of character counts in lexicographical order
// Assumes alphabet is sorted
// Counts = vector with 256 elements
template<class t_bitvector>
void count_smaller_chars(const wt_huff<t_bitvector>& bwt, vector<uint8_t>& alphabet, vector<int64_t>& counts, Interval I){
    assert(alphabet.size() != 0);
    counts[alphabet[0]] = 0;
    if(I.size() == 0){
        for(int64_t i = 1; i < alphabet.size(); i++)
            counts[alphabet[i]] = 0;
        return;
    }
    
    for(int64_t i = 1; i < alphabet.size(); i++){
        int64_t count_prev = bwt.rank(I.right+1, alphabet[i-1]) - bwt.rank(I.left, alphabet[i-1]) ;
        counts[alphabet[i]] = counts[alphabet[i-1]] + count_prev;
    }
}

template<class t_bitvector>
bool BD_BWT_index<t_bitvector>::is_right_maximal(Interval_pair I){
    
    // An interval is right-maximal iff it has more than one possible right extension
    vector<uint8_t> symbols = get_interval_symbols(reverse_bwt, I.reverse);
    return (symbols.size() >= 2 || (I.reverse.size() >= 2 && symbols[0] == '$')); // TODO: don't hardcode the dollar
}

template<class t_bitvector>
bool BD_BWT_index<t_bitvector>::is_left_maximal(Interval_pair I){ // UNTESTED
    
    // An interval is left-maximal iff it has more than one possible left extension
    vector<uint8_t> symbols = get_interval_symbols(forward_bwt, I.forward);
    return (symbols.size() >= 2 || (I.forward.size() >= 2 && symbols[0] == '$')); // TODO: don't hardcode the dollar
}

// Returns the alphabet sorted order
vector<uint8_t> get_string_alphabet(const string& s){
    
    vector<bool> found(256,false);
    for(uint8_t c : s) found[c] = true;
    vector<uint8_t> alphabet;
    for(int i = 0; i < 256; i++){
        if(found[i]) alphabet.push_back((char)i);
    }
    
    return alphabet;
}

template<class t_bitvector>
BD_BWT_index<t_bitvector>::BD_BWT_index(const string& input){
    if(input == "") throw std::runtime_error("Tried to construct BD_BWT_index for an empty string");
    this->alphabet = get_string_alphabet(input);
    
    // If the dollar is not in the input already, add it to the alphabet
    if(find(alphabet.begin(), alphabet.end(), '$') == alphabet.end()){
        alphabet.push_back('$');
        sort(alphabet.begin(), alphabet.end());
    }
    
    // Build the two bwts
    char* forward = (char*) malloc(input.size() + 1);
    char* backward = (char*) malloc(input.size() + 1);
    for(int i = 0; i < input.size(); i++){
        forward[i] = input[i];
        backward[input.size()-1-i] = input[i];
    }
    forward[input.size()] = 0;
    backward[input.size()] = 0;

    char* forward_transform = bwt_dbwt(forward);
    char* backward_transform = bwt_dbwt(backward);
    free(forward);
    free(backward);
    
    // Build wavelet trees
    construct_im(this->forward_bwt, string(forward_transform).c_str(), 1); // Not casting through std::string segfaults
    construct_im(this->reverse_bwt, string(backward_transform).c_str(), 1); // Not casting through std::string segfaults
    free(forward_transform);
    free(backward_transform);
    
    // Compute cumulative character counts
    this->cumulative_char_count.resize(256);
    count_smaller_chars(forward_bwt,alphabet,cumulative_char_count,Interval(0,forward_bwt.size()-1));
}

template<class t_bitvector>
bool allDistinct(const wt_huff<t_bitvector>& bwt, Interval I){
    return get_interval_symbols<t_bitvector>(bwt,I).size() == I.size();
}

template<class t_bitvector>
bool BD_BWT_index<t_bitvector>::interval_is_supermaximal(Interval_pair I){
    return allDistinct(forward_bwt,I.forward) && allDistinct(reverse_bwt,I.reverse);
}

