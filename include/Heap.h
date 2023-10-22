/**
 * Partially adapted from https://github.com/facebookresearch/faiss
 */

/*
 * C++ support for heaps. The set of functions is tailored for efficient
 * similarity search.
 *
 * There is no specific object for a heap, and the functions that operate on a
 * single heap are inlined, because heaps are often small. More complex
 * functions are implemented in Heaps.cpp
 *
 * All heap functions rely on a C template class that define the type of the
 * keys and values and their ordering (increasing with CMax and decreasing with
 * Cmin). The C types are defined in ordered_key_value.h
 */

#ifndef Heap_h
#define Heap_h

#include <climits>
#include <cmath>
#include <cstring>

#include <stdint.h>
#include <cassert>
#include <cstdio>

#include <limits>

#include <utility>

template<typename T_, typename TI_>
struct CMax;

template<typename T_, typename TI_>
struct CMax {
    typedef T_ T;
    typedef TI_ TI;
    inline static bool cmp(T a, T b) {
        return a > b;
    }
};

/*******************************************************************
 * Basic heap ops: push and pop
 *******************************************************************/

/** Pops the top element from the heap defined by bh_val[0..k-1] and
 * bh_ids[0..k-1].  on output the element at k-1 is undefined.
 */
template <class C>
inline void heap_pop(size_t k, typename C::T* bh_val, typename C::TI* bh_ids) {
    bh_val--; /* Use 1-based indexing for easier node->child translation */
    bh_ids--;
    typename C::T val = bh_val[k];
    typename C::TI id = bh_ids[k];
    size_t i = 1, i1, i2;
    while (1) {
        i1 = i << 1;
        i2 = i1 + 1;
        if (i1 > k)
            break;
        if ((i2 == k + 1) ||
            C::cmp(bh_val[i1], bh_val[i2])) {
            if (C::cmp(val, bh_val[i1])) {
                break;
            }
            bh_val[i] = bh_val[i1];
            bh_ids[i] = bh_ids[i1];
            i = i1;
        } else {
            if (C::cmp(val, bh_val[i2])) {
                break;
            }
            bh_val[i] = bh_val[i2];
            bh_ids[i] = bh_ids[i2];
            i = i2;
        }
    }
    bh_val[i] = bh_val[k];
    bh_ids[i] = bh_ids[k];
}

/** Pushes the element (val, ids) into the heap bh_val[0..k-2] and
 * bh_ids[0..k-2].  on output the element at k-1 is defined.
 */
template <class C>
inline void heap_push(
        size_t k,
        typename C::T* bh_val,
        typename C::TI* bh_ids,
        typename C::T val,
        typename C::TI id) {
    bh_val--; /* Use 1-based indexing for easier node->child translation */
    bh_ids--;
    size_t i = k, i_father;
    while (i > 1) {
        i_father = i >> 1;
        if (!C::cmp(val, bh_val[i_father])) {
            /* the heap structure is ok */
            break;
        }
        bh_val[i] = bh_val[i_father];
        bh_ids[i] = bh_ids[i_father];
        i = i_father;
    }
    bh_val[i] = val;
    bh_ids[i] = id;
}


/*******************************************************************
 * Operations on heap arrays
 *******************************************************************/

/** a template structure for a set of [min|max]-heaps it is tailored
 * so that the actual data of the heaps can just live in compact
 * arrays.
 */
template <typename C>
struct HeapArray {
    typedef typename C::TI TI;
    typedef typename C::T T;

    size_t k;  ///< allocated size per heap
    TI* ids;   ///< identifiers (size nh * k)
    T* val;    ///< values (distances or similarities), size nh * k
};


#endif /* Heap_h */
