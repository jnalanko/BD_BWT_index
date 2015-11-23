#ifndef SDSL_ITERATE_HH
#define SDSL_ITERATE_HH

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <vector>
#include <utility>
#include <deque>
#include <string>
#include "Interval.hh"

// Bidirectional BWT index
template<class t_bitvector>
class BD_BWT_index{
private:
    sdsl::wt_huff<t_bitvector> forward_bwt;
    sdsl::wt_huff<t_bitvector> reverse_bwt;
    
    // cumulative_char_count[i] = total number of characters in the string with ascii value strictly less than i
    std::vector<int64_t> cumulative_char_count;
    std::vector<uint8_t> alphabet;
    
    Interval_pair left_extend(Interval_pair intervals, char c, int64_t cumul_rank_c);
    Interval_pair right_extend(Interval_pair intervals, char c, int64_t cumul_rank_c);
public:
    
    BD_BWT_index(std::string& input);
    
    int64_t size() { return forward_bwt.size(); }
    char forward_bwt_at(int64_t index) { return forward_bwt[index]; }
    char backward_bwt_at(int64_t index) { return reverse_bwt[index]; }
    std::vector<int64_t> get_cumulative_char_count() { return cumulative_char_count; }
    const std::vector<uint8_t>& get_alphabet() { return alphabet; }

    // Returns an interval of size zero if extension not possible or if the given interval has size 0
    Interval_pair left_extend(Interval_pair intervals, char c);

    // Returns an interval of size zero if extension not possible or if the given interval has size 0
    Interval_pair right_extend(Interval_pair intervals, char c);

    int64_t backward_step(int64_t pos);

    bool is_right_maximal(Interval_pair I);
    bool is_left_maximal(Interval_pair I);
    bool interval_is_supermaximal(Interval_pair I);

};

#include "BD_BWT_index_impl.hh"

#endif
