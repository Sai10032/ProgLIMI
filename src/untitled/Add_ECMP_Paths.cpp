#include "Add_ECMP_Paths.h"

void Compute_ECMP_Paths_From_Between_Two_Nodes(BGraph& bgraph,unsigned s, unsigned d,list<PathType>& paths)
{
	if(s == d)
	{
		cout<<"计算ECMP路径时源节点和目的节点相同！！！！"<<endl;
		int i;
		cin>>i;
		return;
	}
	std::vector<double> distance(num_vertices(bgraph));
	property_map < BGraph, edge_weight_t>::type
        edgeweight = get(edge_weight, bgraph);
	std::vector<Vertex_d> parent(boost::num_vertices(bgraph));
	boost::dijkstra_shortest_paths(bgraph, d, predecessor_map(boost::make_iterator_property_map(parent.begin(), get(boost::vertex_index, bgraph))).
                          distance_map(boost::make_iterator_property_map(distance.begin(), get(boost::vertex_index, bgraph))));
	
	map<unsigned,list<PathType>> pathrecorder;
	list<PathType> emptypath;
	PathType pp;
	pp.push_back(s);
	emptypath.push_back(pp);
	pathrecorder.insert(make_pair(s,emptypath));
	std::priority_queue<pair<double,unsigned>> vistnodes;
	vector<bool> visited(num_vertices(bgraph),false);
	vistnodes.push(make_pair(distance[s],s));
	visited[s] = true;
	while(!vistnodes.empty())
	{
		pair<double,unsigned> node = vistnodes.top();
		vistnodes.pop();
		Adj_Vertex_i adj_b, adj_e;
		boost::tie(adj_b,adj_e) = boost::adjacent_vertices(node.second,bgraph);
		for(;adj_b != adj_e; adj_b++)
		{
			Edge_d ed; 
			bool suc;
			boost::tie(ed,suc) = edge(node.second,*adj_b,bgraph);
			double erro = edgeweight[ed] + distance[*adj_b] - distance[node.second];
			if(erro < 0.0000001)
			{
				list<PathType>::iterator pathiter = pathrecorder[node.second].begin();
				if(!visited[*adj_b])
				{
					vistnodes.push(make_pair(distance[*adj_b],*adj_b));
					visited[*adj_b] = true;
				}
				for(;pathiter != pathrecorder[node.second].end();pathiter++)
				{
					PathType pathtemp = *pathiter;
					pathtemp.push_back(*adj_b);
					pathrecorder[*adj_b].push_back(pathtemp);
				}
			}
		}
	}
	paths = pathrecorder[d];

	return;
}
