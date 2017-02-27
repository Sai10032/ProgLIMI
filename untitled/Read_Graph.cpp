#include "Read_Graph.h"

BGraph ReadTopology(string Filename) {

    BGraph bgraph;
    ifstream fin;
    ofstream fout;
    fin.open(Filename.c_str());
    if (!fin) {
        cerr << "文件打开失败" << endl;
        exit(1);
    }
    int numSDNnodes;
    int SDNnode;
    vector<int> SDNnodes;
    fin >> numSDNnodes;
    for (int i = 0; i < numSDNnodes; i++) {

        fin >> SDNnode;
        SDNnodes.push_back(SDNnode);
    }

    unsigned int source;
    unsigned int target;
    unsigned int numRedundantEdge;
    numRedundantEdge = 0;
    while (fin >> source >> target) {
        bool suc;
        Edge_d ed;
        if ((source < num_vertices(bgraph)) && (target < num_vertices(bgraph))) {
            boost::tie(ed, suc) = edge(source, target, bgraph);
            if (!suc) {
                add_edge(source, target, bgraph);
            } else
                numRedundantEdge++;
        } else {
            add_edge(source, target, bgraph);
        }
    }
    cout << "冗余边数目：" << numRedundantEdge << endl;

    property_map<BGraph, is_SDN_t>::type Is_SDN = get(is_SDN_t(), bgraph);
    Vertex_i vb, ve;
    boost::tie(vb, ve) = vertices(bgraph);
    for (; vb != ve; vb++) {
        Is_SDN[*vb] = false;
    }
    for (int i = 0; i < numSDNnodes; i++) {
        Is_SDN[SDNnodes[i]] = true;
    }
    fin.close();

    return bgraph;
}