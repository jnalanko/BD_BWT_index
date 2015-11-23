#ifndef SDSL_ITERATE_HH
#define SDSL_ITERATE_HH

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <vector>
#include <utility>
#include <deque>
#include <string>
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

#include "BD_BWT_index_impl.hh"

#endif
