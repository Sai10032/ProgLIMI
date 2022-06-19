#ifndef _BOOST_GRAPH_H
#define _BOOST_GRAPH_H
#include <algorithm>
#include <bitset>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/utility.hpp>
#include <iostream>
#include <iterator>
#include <queue>
#include <string>
#include <utility>
#include <vector>
// using namespace ogdf;
using namespace std;
using namespace boost;

#define INF 1000000
struct node_id_t {
    typedef vertex_property_tag kind;
};

struct real_id_t {
    typedef vertex_property_tag kind;
};
struct parent_t {
    typedef vertex_property_tag kind;
};
struct is_monitor_t {
    typedef vertex_property_tag kind;
};
struct is_virtual_t {
    typedef vertex_property_tag kind;
};

struct is_SDN_t {
    typedef vertex_property_tag kind;
};

struct order_t {
    typedef vertex_property_tag kind;
};

struct weight_t {
    typedef edge_property_tag kind;
};

struct mark_t {
    typedef edge_property_tag kind;
};
struct link_index_t {
    typedef edge_property_tag kind;
};

struct edge_component_t {
    enum { num = 555 };
    typedef edge_property_tag kind;
};

typedef property<
    node_id_t, int,
    property<
        is_monitor_t, bool,
        property<is_SDN_t, bool,
                 property<is_virtual_t, bool,
                          property<real_id_t, unsigned, property<order_t, double, property<parent_t, unsigned>>>>>>>
    Vproperty;

typedef property<edge_weight_t, int,
                 property<edge_component_t, std::size_t, property<link_index_t, unsigned, property<mark_t, bool>>>>
    Eproperty;
typedef adjacency_list<vecS, vecS, undirectedS, Vproperty, Eproperty> BGraph;
typedef graph_traits<BGraph>::vertex_descriptor Vertex_d;
typedef graph_traits<BGraph>::vertex_iterator Vertex_i;
typedef graph_traits<BGraph>::edge_descriptor Edge_d;
typedef graph_traits<BGraph>::edge_iterator Edge_i;
typedef graph_traits<BGraph>::out_edge_iterator Gout_Edge_i;
typedef graph_traits<BGraph>::adjacency_iterator Adj_Vertex_i;
typedef int Weight_T;
typedef list<int> PathType;
typedef list<int>::iterator PathType_iter;
void Compute_Shortest_Path_From_OneNode(BGraph& bgraph, Vertex_d& s, Vertex_d& d, PathType& path);
void Set_Weight_to_One(BGraph& bg);
#endif