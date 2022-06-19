#include "Boost_Graph.h"
#include "Read_Graph.h"
#include "Add_ECMP_Paths.h"
//#include <ogdf/basic/Graph.h>
#include <boost/numeric/ublas/matrix.hpp>
//using namespace ogdf;

void Set_SDNnodes_As_Monitors(BGraph& bgraph);
void Get_Monitors(BGraph bgraph,list<Vertex_d>& Moniters);
void Get_SDNnodes(BGraph bgraph,list<Vertex_d>& SDNnodes);
void Set_Link_Cost(BGraph& bgraph);
void Get_Vertex_Index(BGraph& bgraph, vector<Vertex_d>& nodeIndex);
void Compute_Shortest_Path_From_OneNode(BGraph& bgraph, Vertex_d& s, Vertex_d& d, PathType& path);

void Set_Graph_Index(BGraph& bgraph);

void Construct_Shortest_Paths_Between_SDNnodes(BGraph& bgraph, multimap<double, PathType>& Indepedent_Paths);
void Put_Edge_Paths_To_Set(BGraph& bgraph,multimap<double,PathType>& Indepedent_Paths);

void Change_Path_to_Matrix(BGraph& bgraph, PathType& path, boost::numeric::ublas::matrix<long double>& matrix1);
long double Compute_Squre( boost::numeric::ublas::matrix<long double>& v);
void Compute_The_QMatrix();
void Compute_R12Matrix(list<unsigned>& record);
void Compute_Prod_Of_RR_and_R12(boost::numeric::ublas::matrix<long double>& temp_matrix);
void Get_Paths_From_Nodes(BGraph& bgraph, list<Vertex_d>& monitors, Vertex_d v, std::queue<PathType>& paths);
int Find_Linear_Independent_Paths_For_One_Node(BGraph& bgraph, std::queue<PathType>& paths, multimap<double, PathType>& Output_Paths);

void Choose_Shortest_Independent_paths(BGraph& bgraph, multimap<double, PathType>& Indepedent_Paths, multimap<double, PathType>& Output_Paths);
void Output_Paths_In_File(BGraph& bgraph, multimap<double, PathType>& Output_Paths, string fileresultname);
void Check_Linear_Independent(BGraph& bgraph, multimap<double, PathType>& Output_Paths, string matrixname);

void Put_One_Step_Circle_Paths_To_Set(BGraph& bgraph,multimap<double,PathType>& Indepedent_Paths);

void Main_Process_For_Path_Selection(string filename, string fileresultname);