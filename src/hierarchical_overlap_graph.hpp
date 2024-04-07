#ifndef HIERARCHICAL_OVERLAP_GRAPH_HPP
#define HIERARCHICAL_OVERLAP_GRAPH_HPP

#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <functional>
#include <vector>
#include <queue>
#include <unordered_map>
#include <stack>
#include <set>

using std::queue;
using std::set;
using std::stack;
using std::string;
using std::unordered_map;
using std::vector;

#include "util.hpp"

class HOG
{
private:
    uint32_t _size{0};

    vector<uint32_t> parent{};
    vector<uint32_t> sibling{};
    vector<uint32_t> failure_link{};

    vector<string> reads{};
    vector<uint32_t> read_id{};
    vector<uint32_t> read_length{};

    vector<bool> marked{}; // bHog in the paper
    vector<bool> leaves{};

    vector<vector<uint32_t>> R_l{};

    uint32_t get_first_child(uint32_t node);
    uint32_t get_last_child(uint32_t node);
    uint32_t get_child(uint32_t node, char c);

    void add_reads();
    void add_failure_links();
    void mark_ehog();

    void contract();

    void mark_Rl();
    void mark_hog();

    string get_label(uint32_t node);
    void print_node(uint32_t);

public:
    HOG(vector<string> &in_reads);
    uint32_t size() { return _size; }
    void print();

}; // end class HOG

#endif // HIERARCHICAL_OVERLAP_GRAPH_HPP