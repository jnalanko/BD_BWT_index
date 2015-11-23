#ifndef ITERATORS_HH
#define ITERATORS_HH

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include "BD_BWT_index.hh"
#include <algorithm>

/**
 * Class BD_BWT_index_iterator
 * 
 * Iterates the suffix link tree of the given index.
 * 
 */
template<class t_bitvector>
class BD_BWT_index_iterator{
    
public:
    
    class Stack_frame{
    public:
        Interval_pair intervals;  // forward interval, reverse interval
        int64_t depth; // depth in the suffix link tree
        char extension; // the label on the arc between this node and its parent
        Stack_frame(Interval_pair intervals, int64_t depth, char extension) : intervals(intervals), depth(depth), extension(extension) {}
        Stack_frame(){}
    };
    
    BD_BWT_index<t_bitvector>* index;
    
    // Iteration state
    std::deque<Stack_frame> iteration_stack;
    Stack_frame current; // A stack frame containing information about the current node
    std::string label; // The string on the path from the root to the current node
    
    // Reused space between iterations
    std::vector<int64_t> char_counts;
    
    BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index) : index(index), char_counts(256) {
        Interval empty_string(0,index->size()-1);
        iteration_stack.push_back(Stack_frame(Interval_pair(empty_string,empty_string), 0, 0));
        current = iteration_stack.back();
        label = "";
    }
    
    /**
     * @brief Go to the next suffix link tree node
     * 
     * Updates the member variables iteration_stack, current and label.
     * 
     * @return False if all nodes have been iterated, else true
     */
    bool next(); // Go to the next node. Return false if none found
    
    /**
     * @brief Go to the next depth k suffix link tree node
     * 
     * Updates the member variables iteration_stack, current and label.
     * 
     * @return False if all nodes have been iterated, else true
     */
    bool next(int64_t k);
    
private:
    void push_right_maximal_children(Stack_frame f);
    void update_label(Stack_frame f);
};


#include "Iterators_impl.hh"

#endif