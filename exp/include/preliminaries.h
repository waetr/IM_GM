//
// Created by asuka on 2022/7/15.
//

#ifndef EXP_MODELS_H
#define EXP_MODELS_H

#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <sstream>
#include <set>
#include <queue>
#include <stack>
#include <bitset>

#define MAX_NODE_SIZE 50000000
#define graph_type int8_t
#define DIRECTED_G 0
#define UNDIRECTED_G 1

#define model_type int8_t
#define NONE 0
#define IC 1
#define LT 2
#define IC_M 3

#define IM_solver int8_t
#define ENUMERATION 0
#define DEGREE 1
#define PAGERANK 2
#define CELF_NORMAL 3
#define DEGREE_ADVANCED 4
#define PAGERANK_ADVANCED 5
#define CELF_ADVANCED 6
#define IMM_NORMAL 7
#define IMM_ADVANCED 8
#define OPIM_NORMAL 9
#define OPIM_ADVANCED 10
#define CELF_THRESHOLD 11
#define CELF_THRESHOLD1 12
#define CELF_THRESHOLD2 13

typedef int64_t int64;
typedef std::pair<int64, int64 > bi_node;

const std::string solver_name[] = {"ENUMERATION", "DEGREE", "PAGERANK", "CELF", "DEGREE_ADVANCED", "PAGERANK_ADVANCED",
                                   "CELF_ADVANCED", "IMM_NORMAL", "IMM_ADVANCED", "OPIM_NORMAL", "OPIM_ADVANCED",
                                   "CELF_THRESHOLD", "CELF_THRESHOLD1", "CELF_THRESHOLD2"};

std::random_device rd__;
std::minstd_rand random_engine(rd__());
std::mt19937 mt19937engine(rd__());
std::uniform_real_distribution<double> distrib(0.0, 1.0);

std::ofstream stdFileOut;
int8_t verbose_flag, local_mg;
int64_t MC_iteration_rounds = 10000;

double logcnk(int n, int k) {
    k = std::min(n, k);
    double ans = 0;
    for (int i = n - k + 1; i <= n; i++) {
        ans += log(i);
    }
    for (int i = 1; i <= k; i++) {
        ans -= log(i);
    }
    return std::max(ans, 0.0);
}

template<class T>
T sqr(T x) {
    return x * x;
}

/*!
 * @brief Random number generator that generates real number between [0,1)
 * @return A random real number between [0,1)
 */
inline double random_real() {
    return distrib(random_engine);
}

/*!
 * @brief calculate the interval from start time.
 * @param start : start timestamp
 * @return the length of the interval
 */
double time_by(double start) {
    return (clock() - start) / CLOCKS_PER_SEC;
}

/*!
 * @brief Print all elements in a vector.
 * @param S : the set
 */
void print_set(std::vector<int64> &S, const std::string &Prefix = "") {
    std::vector<int64> S_ordered = S;
    std::sort(S_ordered.begin(), S_ordered.end());
    std::cout << Prefix;
    std::cout << "{";
    for (int64 i = 0; i < S_ordered.size(); i++) {
        std::cout << S_ordered[i];
        if (i != S_ordered.size() - 1) std::cout << ",";
    }
    std::cout << "}[" << S_ordered.size() << "]";
}

/*!
 * @brief Print the set if the size of the set is not exceeding 100.
 * @param S : the set
 */
void print_small_set(std::vector<int64> &S, const std::string &Prefix = "") {
    if (S.size() <= 100) print_set(S, Prefix);
    else {
        std::cout << Prefix;
        std::cout << "{...}[" << S.size() << "]";
    }
}

/*!
 * @brief Print all elements to the file in a vector.
 * @param S : the set
 */
void print_set_f(std::vector<int64> &S, const std::string &Prefix = "") {
    if (S.size() > 100) return;
    std::vector<int64> S_ordered = S;
    std::sort(S_ordered.begin(), S_ordered.end());
    stdFileOut << Prefix;
    stdFileOut << "{";
    for (int64 i = 0; i < S_ordered.size(); i++) {
        stdFileOut << S_ordered[i];
        if (i != S_ordered.size() - 1) stdFileOut << ",";
    }
    stdFileOut << "}[" << S_ordered.size() << "]";
}

inline int64 fast_read(std::ifstream &inFile) {
    int64 x = 0;
    auto ch = inFile.get();
    if (ch == EOF) return -1;
    while (ch < '0' || ch > '9') {
        if (ch == EOF) return -1;
        ch = inFile.get();
    }
    while (ch >= '0' && ch <= '9') {
        x = x * 10 + ch - '0';
        ch = inFile.get();
    }
    return x;
}

void convert_binode_into_node(std::vector<bi_node> &source, std::vector<int64> &target) {
    target.clear();
    for(auto e : source) target.emplace_back(e.first);
}

void convert_node_into_binode(std::vector<int64> &source, std::vector<bi_node> &target) {
    target.clear();
    for(auto e : source) target.emplace_back(e, -1);
}

#endif //EXP_MODELS_H