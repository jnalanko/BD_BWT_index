#include <vector>
#include "BD_BWT_index.hh"
#include "Iterators.hh"
#include <iostream>

template<class t_bitvector>
void BD_BWT_index_iterator<t_bitvector>::push_right_maximal_children(Stack_frame f){
    for(uint8_t c : index->get_alphabet()){
        if(c == BD_BWT_index<t_bitvector>::END) continue;
        Interval_pair child = index->left_extend(f.intervals,c);
        if(child.forward.size() == 0) continue; // Extension not possible
        if(index->is_right_maximal(child)){
            // Add child to stack
            iteration_stack.push_back(Stack_frame(child,f.depth+1,c));
        }
    }    
}

template<class t_bitvector>
void BD_BWT_index_iterator<t_bitvector>::update_label(Stack_frame f){
    while(label.size() > 0 && label.size() >= f.depth) // Unwind stack
        label.pop_back();
    if(f.extension != 0)
        label.push_back(current.extension);    
}

template<class t_bitvector>
bool BD_BWT_index_iterator<t_bitvector>::next(int64_t k){
    
    while(true){
        if(iteration_stack.empty()) return false;
        
        current = iteration_stack.back(); iteration_stack.pop_back();
        update_label(current);
        
        if(current.depth == k) return true; // Stop recursing to children and give control back

        push_right_maximal_children(current);

    }
}

template<class t_bitvector>
bool BD_BWT_index_iterator<t_bitvector>::next(){
    if(iteration_stack.empty()) return false;
    
    current = iteration_stack.back(); iteration_stack.pop_back();
    update_label(current);
    push_right_maximal_children(current);
        
    return true;
}

