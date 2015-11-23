#ifndef SDSL_ITERATE_HH
#define SDSL_ITERATE_HH

#include <sdsl/construct.hpp>
#include <vector>
#include <utility>
#include <string>
#include "bwt.hh"
#include "Interval.hh"

/*
 * Implements a bidictional BWT index.
 * All indices and ranks are indexed starting from zero.
 * All intervals are inclusive, i.e. the interval from i to j includes both i and j.
 */

template<class t_bitvector>
class BD_BWT_index{
    
private:
        
    sdsl::wt_huff<t_bitvector> forward_bwt;
    sdsl::wt_huff<t_bitvector> reverse_bwt;
    
    // cumulative_char_count[i] = total number of characters in the string with ascii value strictly less than i
    std::vector<int64_t> cumulative_char_count;
    std::vector<uint8_t> alphabet;
    
    Interval_pair left_extend(Interval_pair intervals, char c, int64_t cumul_rank_c) const;
    Interval_pair right_extend(Interval_pair intervals, char c, int64_t cumul_rank_c) const;
    std::vector<uint8_t> get_string_alphabet(const uint8_t* s) const;
    int64_t strlen(const uint8_t* str) const;
    int64_t compute_cumulative_char_rank_in_interval(const sdsl::wt_huff<t_bitvector>& wt, char c, Interval I) const;
    std::vector<uint8_t> get_interval_symbols(const sdsl::wt_huff<t_bitvector>& wt, Interval I) const;
    void get_interval_symbols(const sdsl::wt_huff<t_bitvector>& wt, Interval I, sdsl::int_vector_size_type& nExtensions, 
                              std::vector<uint8_t>& symbols, std::vector<uint64_t>& ranks_i, std::vector<uint64_t>& ranks_j) const;
    void count_smaller_chars(const sdsl::wt_huff<t_bitvector>& bwt, std::vector<uint8_t>& alphabet, std::vector<int64_t>& counts, Interval I) const;

public:

    // The input string must not contain the END byte
    static const uint8_t END = 0x01; // End of string marker.
    BD_BWT_index(const uint8_t* input);
    
    int64_t size() const { return forward_bwt.size();}
    uint8_t forward_bwt_at(int64_t index) const { return forward_bwt[index]; }
    uint8_t backward_bwt_at(int64_t index) const { return reverse_bwt[index]; }
    const std::vector<int64_t>& get_cumulative_char_count() const { return cumulative_char_count; }
    const std::vector<uint8_t>& get_alphabet() const { return alphabet; }

    // Returns an interval of size zero if extension not possible or if the given interval has size 0
    Interval_pair left_extend(Interval_pair intervals, char c) const;

    // Returns an interval of size zero if extension not possible or if the given interval has size 0
    Interval_pair right_extend(Interval_pair intervals, char c) const;

    // Let s be a suffix of length k with lexicographic rank lex_rank among all suffixes of the input string.
    // Returns the lexicographic rank of the suffix with length k+1 if k is not equal to the length of the
    // input string counting the terminating symbol, or the lexicographic rank of the suffix with length 0 otherwise.
    int64_t backward_step(int64_t lex_rank) const;
    
    // Let s be a prefix of length k with colexicographic rank colex_rank among all prefixes of the input string.
    // Returns the colexicographical rank of the prefix with length k+1 if k is not equal to the length of the
    // input string counting the terminating symbol, or the colexicographic rank of the prefix with length 0 otherwise.
    int64_t forward_step(int64_t colex_rank) const;
    
    bool is_right_maximal(Interval_pair I) const;
    bool is_left_maximal(Interval_pair I) const;
    bool interval_is_supermaximal(Interval_pair I) const;

};

template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::compute_cumulative_char_rank_in_interval(const sdsl::wt_huff<t_bitvector>& wt, char c, Interval I) const{
    int64_t ans = 0;
    if(I.size() == 0) return 0;
    
    // Sum of ranks of all characters that are lexicographically smaller than c
    for(uint8_t d : alphabet){
        if(d == c) break;
        ans += wt.rank(I.right + 1,d) - wt.rank(I.left,d);
    }
    return ans;
}


template<class t_bitvector>
std::vector<uint8_t> BD_BWT_index<t_bitvector>::get_interval_symbols(const sdsl::wt_huff<t_bitvector>& wt, Interval I) const{
    if(I.size() == 0){
        std::vector<uint8_t> empty;
        return empty;
    }
    
    sdsl::int_vector_size_type nExtensions;
    std::vector<uint8_t> symbols(wt.sigma);
    std::vector<uint64_t> ranks_i(wt.sigma);
    std::vector<uint64_t> ranks_j(wt.sigma);
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
void BD_BWT_index<t_bitvector>::get_interval_symbols(const sdsl::wt_huff<t_bitvector>& wt, Interval I, sdsl::int_vector_size_type& nExtensions, std::vector<uint8_t>& symbols,
 std::vector<uint64_t>& ranks_i, std::vector<uint64_t>& ranks_j) const{
    if(I.size() == 0){
        nExtensions = 0;
        return;
    }
    
    wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
}


// Takes a backward step in the forward bwt
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::backward_step(int64_t lex_rank) const{
    uint8_t c = forward_bwt[lex_rank];
    return cumulative_char_count[c] + forward_bwt.rank(lex_rank, c);
}

// Takes a backward step in the reverse bwt
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::forward_step(int64_t colex_rank) const{
    uint8_t c = reverse_bwt[colex_rank];
    return cumulative_char_count[c] + reverse_bwt.rank(colex_rank, c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::left_extend(Interval_pair intervals, char c) const{
    int64_t cumul_rank_c = compute_cumulative_char_rank_in_interval(forward_bwt, c, intervals.forward);
    return left_extend(intervals,c,cumul_rank_c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::right_extend(Interval_pair intervals, char c) const{
    int64_t cumul_rank_c = compute_cumulative_char_rank_in_interval(reverse_bwt, c, intervals.reverse);
    return right_extend(intervals,c,cumul_rank_c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::left_extend(Interval_pair intervals, char c, int64_t cumul_rank_c) const{
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
Interval_pair BD_BWT_index<t_bitvector>::right_extend(Interval_pair intervals, char c, int64_t cumul_rank_c) const{
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
void BD_BWT_index<t_bitvector>::count_smaller_chars(const sdsl::wt_huff<t_bitvector>& bwt, std::vector<uint8_t>& alphabet, std::vector<int64_t>& counts, Interval I) const{
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
bool BD_BWT_index<t_bitvector>::is_right_maximal(Interval_pair I) const{
    
    // An interval is right-maximal iff it has more than one possible right extension
    std::vector<uint8_t> symbols = get_interval_symbols(reverse_bwt, I.reverse);
    return (symbols.size() >= 2);
}

template<class t_bitvector>
bool BD_BWT_index<t_bitvector>::is_left_maximal(Interval_pair I) const{
    
    // An interval is left-maximal iff it has more than one possible left extension
    std::vector<uint8_t> symbols = get_interval_symbols(forward_bwt, I.forward);
    return (symbols.size() >= 2);
}

// Returns the alphabet in sorted order
template<class t_bitvector>
std::vector<uint8_t> BD_BWT_index<t_bitvector>::get_string_alphabet(const uint8_t* s) const{
    
    std::vector<bool> found(256,false);
    while(*s != 0){
        found[*s] = true;
        s++;
    }

    std::vector<uint8_t> alphabet;
    for(int i = 0; i < 256; i++){
        if(found[i]) alphabet.push_back((uint8_t)i);
    }
    
    return alphabet;
}

// strlen(const uint8_t*) is not in the standard library
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::strlen(const uint8_t* str) const{
    const uint8_t* start = str;
    while(*str != 0) str++;
    return str - start;
}

template<class t_bitvector>
BD_BWT_index<t_bitvector>::BD_BWT_index(const uint8_t* input){
    if(*input == 0) throw std::runtime_error("Tried to construct BD_BWT_index for an empty string");
    int64_t n = strlen(input);
    
    if(std::find(input, input+n, END) != input + n){
        std::stringstream error;
        error << "Input string contains forbidden byte " << std::hex << END;
        throw std::runtime_error(error.str());
    }
    
    // Build the two bwts
    uint8_t* forward = (uint8_t*) malloc(sizeof(uint8_t) * (n + 1));
    uint8_t* backward = (uint8_t*) malloc(sizeof(uint8_t) * (n + 1));
    for(int64_t i = 0; i < n; i++){
        forward[i] = input[i];
        backward[n-1-i] = input[i];
    }
    forward[n] = END;
    backward[n] = END;

    uint8_t* forward_transform = bwt_dbwt(forward,n,END);
    free(forward);
        
    uint8_t* backward_transform = bwt_dbwt(backward,n,END);
    free(backward);
    
    // Build wavelet trees
    construct_im(this->forward_bwt, (const char*)forward_transform, 1); // Must cast to signed char* or else breaks. File a bug report to sdsl?
    construct_im(this->reverse_bwt, (const char*)backward_transform, 1); // Must cast to signed char* or else breaks. File a bug report to sdsl?
    
    this->alphabet = get_string_alphabet(forward_transform);
    
    free(forward_transform);
    free(backward_transform); 
    
    // Compute cumulative character counts
    this->cumulative_char_count.resize(256);
    count_smaller_chars(forward_bwt,alphabet,cumulative_char_count,Interval(0,forward_bwt.size()-1));
    
}

#endif
