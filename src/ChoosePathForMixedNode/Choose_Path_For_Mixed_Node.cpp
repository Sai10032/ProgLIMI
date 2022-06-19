#include <boost/graph/biconnected_components.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include "Boost_Graph.h"
#include "Path_Selection.h"

int main(int /*argc*/, char* /*argv*/[]) {
    string filename;
    string fileresultname;
    string matrixname;

    vector<string> vecIndex;

    string prefixTopo;
    string prefixResult;

    prefixTopo = "GEANT23_SDN_Hybrid/";
    prefixResult = "GEANT23_SDN_Hybrid_result/";

    for (unsigned index = 0; index <= 10; index++) {
        char strIndex[4];
        string tmpstr;
        sprintf(strIndex, "%d", index * 10);
        tmpstr = strIndex;

        filename = prefixTopo + tmpstr + "_percent.txt";
        fileresultname = prefixResult + tmpstr + "_percent.txt";

        cout << "Processing file: " << filename << endl;
        cout << "Target file: " << fileresultname << endl;
        Main_Process_For_Path_Selection(filename, fileresultname);
        cout << "*************************" << endl;
    }

    return 0;
}
