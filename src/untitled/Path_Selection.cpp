#include "Path_Selection.h"

boost::numeric::ublas::matrix<long double> R, RR, R12,G,Q;
long double R22;
vector<list<unsigned>> G1;

void Set_SDNnodes_As_Monitors(BGraph& bgraph)
{
	property_map < BGraph, is_SDN_t>::type Is_SDN = get(is_SDN_t(), bgraph);
	property_map < BGraph, is_monitor_t>::type Is_monitor = get(is_monitor_t(), bgraph);

	Vertex_i vb,ve;
	boost::tie(vb,ve) = vertices(bgraph);
	for(;vb!=ve;vb++){
		if(Is_SDN[*vb] == true){
			Is_monitor[*vb] = true;
		}else{
			Is_monitor[*vb] = false;
		}
	}
}

void Get_Monitors(BGraph bgraph, list<Vertex_d>& Moniters)
{
	property_map < BGraph, is_monitor_t>::type Is_monitor = get(is_monitor_t(), bgraph);
	Vertex_i vb,ve;
	boost::tie(vb,ve) = vertices(bgraph);
	for(;vb!=ve;vb++){
		if(Is_monitor[*vb] == true){
			Moniters.push_back(*vb);
		}
	}
}

void Get_SDNnodes(BGraph bgraph, list<Vertex_d>& SDNnodes)
{
	property_map < BGraph, is_SDN_t>::type Is_SDN = get(is_SDN_t(), bgraph);
	Vertex_i vb,ve;
	boost::tie(vb,ve) = vertices(bgraph);
	for(;vb!=ve;vb++){
		if(Is_SDN[*vb] == true){
			SDNnodes.push_back(*vb);
		}
	}
}
void Set_Link_Cost(BGraph& bgraph)
{
	Edge_i eb,ee;
	property_map < BGraph, edge_weight_t>::type Edge_Weight = get(edge_weight_t(), bgraph);
	boost::tie(eb,ee) = edges(bgraph);
	int c = 0;
	for(;eb != ee; eb++){
		Edge_Weight[*eb] = 1000 + c;
		c++;
	}
}

void Get_Vertex_Index(BGraph& bgraph, vector<Vertex_d>& nodeIndex)
{
	Vertex_i vb,ve;
	boost::tie(vb,ve) = vertices(bgraph);
	for(;vb!=ve;vb++){
		nodeIndex.push_back(*vb);
	}
}
bool Is_SDNnodes_In_Path(list<Vertex_d>& SDNnodes, PathType path)
{
	PathType_iter pb, pe;
	list<Vertex_d>::iterator ib, ie;
	path.pop_back();
	if(path.empty()){
		return false;
	}
	pe = path.end();
	ib = SDNnodes.begin();
	ie = SDNnodes.end();
	for(pb = path.begin(); pb!=pe; pb++){
		if(find(ib, ie, *pb)!= ie){
			return true;
		}else{
			continue;
		}
	}
	return false;
}

void Compute_Shortest_Path_From_OneNode(BGraph& bgraph, Vertex_d& s, Vertex_d& d, PathType& path)
{
	if(s == d)
		return;
	std::vector<Vertex_d> parent(boost::num_vertices(bgraph));
	std::vector<int> distance(num_vertices(bgraph));
	Vertex_d next_parent,last_parent;  
	
	boost::dijkstra_shortest_paths(bgraph, s, predecessor_map(boost::make_iterator_property_map(parent.begin(), get(boost::vertex_index, bgraph))).
                          distance_map(boost::make_iterator_property_map(distance.begin(), get(boost::vertex_index, bgraph))));
	vector<Vertex_d>::iterator v_iter;
	
	last_parent = d;
	path.push_back(last_parent);
	do{
		next_parent=parent[last_parent];
	    last_parent = next_parent;
		path.push_front(next_parent);
	}
	while(next_parent != s);
}

void Construct_Shortest_Paths_Between_SDNnodes(BGraph& bgraph, multimap<double, PathType>& Indepedent_Paths)
{
	list<Vertex_d> SDNnodes;
	list<Vertex_d>::iterator viter1, viter2;
	vector<Vertex_d> nodeIndex;
	PathType path;
	
	Set_Link_Cost(bgraph);	
	//��ȡ��������
	Get_Vertex_Index(bgraph, nodeIndex);
	//get SDN nodes
	Get_SDNnodes(bgraph, SDNnodes);
	SDNnodes.sort();

	Adj_Vertex_i adj_i, adj_ie;
	for(viter1 = SDNnodes.begin(); viter1 != SDNnodes.end(); viter1++){
		for(viter2 = SDNnodes.begin(); viter2 != SDNnodes.end(); viter2++){
			if (*viter1 != *viter2){
				boost::tie(adj_i, adj_ie) = adjacent_vertices(*viter1, bgraph);
				for(; adj_i != adj_ie; adj_i++){
					if(nodeIndex[*adj_i] == *viter2){
						path.push_back(*viter1);
						path.push_back(*viter2);
						Indepedent_Paths.insert(make_pair(1, path));
						path.clear();
					}else{
						int k = nodeIndex[*adj_i];
						Compute_Shortest_Path_From_OneNode(bgraph, nodeIndex[*adj_i], *viter2, path);
						//�ж�path�ϳ���Ŀ�Ľ���Ƿ���SDN���
						if( Is_SDNnodes_In_Path(SDNnodes, path) ){
							path.clear();
							continue;
						}else{
							path.push_front(*viter1);
							Indepedent_Paths.insert(make_pair(1, path));
							path.clear();
						}
					}
				}
			}
		}
	}
	//�Ƿ���Ҫȥ��

}

void Set_Graph_Index(BGraph& bgraph)
{
	Edge_i eb, ee;
    unsigned i = 0;
	property_map < BGraph, link_index_t>::type linkindex = get(link_index_t(), bgraph);
	for(boost::tie(eb,ee) = edges(bgraph); eb != ee; eb++){
		linkindex[*eb] = i;
		i++;
	}
}

void Put_Edge_Paths_To_Set(BGraph& bgraph,multimap<double,PathType>& Indepedent_Paths)
{
	Edge_i eb, ee;
	boost::tie(eb,ee) = edges(bgraph);
	for(; eb != ee; eb++){
		Vertex_d s,d;
		s = source(*eb, bgraph);
		d = target(*eb, bgraph);
		PathType path;
		path.push_back(s);
		path.push_back(d);
		Indepedent_Paths.insert(make_pair(3, path));
		path.clear();
	}
}

void Change_Path_to_Matrix(BGraph& bgraph,PathType& path, boost::numeric::ublas::matrix<long double>& matrix1)
{
	property_map < BGraph, link_index_t>::type linkindex = get(link_index_t(), bgraph);
	unsigned index;

	PathType_iter p_iter1, p_iter2;
	Edge_d ed;
	bool Is_Suc = false;
	for(unsigned col = 0; col < matrix1.size2(); col++){
		matrix1(0, col) = 0;
	}
	for(p_iter1 = path.begin();p_iter1 != path.end(); p_iter1++){
		if (p_iter1 == path.begin()){
			p_iter2 = p_iter1;
			continue;
		}
		boost::tie(ed, Is_Suc) = boost::edge(*p_iter2, *p_iter1, bgraph);
		p_iter2 = p_iter1;
		if(Is_Suc == true){
			index = linkindex[ed];
//			matrix1(0, index) = matrix1(0,index) + 1; // �޸Ĺ���ʹ������֧���ظ��ı�
			matrix1(0, index) = 1;
		}else{
			matrix1(0,index)  = 0;
		}
	}
	//output matrix to console
//	for(int i=0; i< matrix1.size1();i++){
//		for(int j=0; j<matrix1.size2(); j++){
//			cout<< matrix1(i, j)<< "  ";
//		}
//		cout<< endl;
//	}
}

long double Compute_Squre( boost::numeric::ublas::matrix<long double>& v)
{
	unsigned row = v.size1();
	unsigned col = v.size2();
	long double  sum = 0;
	for(unsigned i = 0; i < row; i++){
		for(unsigned j = 0; j < col; j++){
			sum = sum + v(i,j)*v(i,j);
			
		}
	}
	//cout<<endl;
	return sum;	
}

void Compute_The_QMatrix()
{
	unsigned rrcol = RR.size2() - 1;
	unsigned gcol;
	Q.resize(RR.size2(), G1.size());
	for(gcol = 0; gcol < G1.size(); gcol++){
		list<unsigned>::iterator Gcoliter;
		long double value = 0;
		for(Gcoliter = G1[gcol].begin(); Gcoliter != G1[gcol].end(); Gcoliter++){
			unsigned index = *Gcoliter;
			value += RR(index, rrcol);
		}
		Q(rrcol, gcol) = value;
	} 
}

void Compute_R12Matrix(list<unsigned>& record)
{
	list<unsigned>::iterator riter;
	/*
	for(unsigned num = 0; num < R12.size1(); num++){
		R12(num,0) = 0;
	}
	*/
	for(unsigned index = 0;index < Q.size1();index++){
		long double value = 0;
		for(riter = record.begin(); riter != record.end(); riter++){
			value += Q(index,*riter);
		}
		R12(index,0) = value;
	}
}

void Compute_Prod_Of_RR_and_R12(boost::numeric::ublas::matrix<long double>& temp_matrix)
{
	list<unsigned> record;
	for(unsigned raw = 0; raw < R12.size1(); raw++){
		double bb = R12(raw,0);
		if(bb < 0){
			bb = -bb;
		}
		if(bb > 0){
			record.push_back(raw); //�´��룬��һ���ڵ�0����һ������Ԫ��
		}
	}	
	list<unsigned>::iterator riter;
	for(unsigned index = 0; index < RR.size1(); index++){
		long double value = 0;
		for(riter = record.begin(); riter != record.end(); riter++){
			value += (RR(index,*riter)*R12(*riter,0));
		}
		temp_matrix(index,0) = value;
	}
}
/*
void Get_Paths_From_Nodes(BGraph& bgraph, list<Vertex_d>& monitors, PathType& path, std::queue<PathType>& paths1, std::queue<PathType>& paths2)
{
	Vertex_d v,w;
	v = *path.begin();
    w = *(++path.begin());
	PathType pathtemp;
	list<Vertex_d>::iterator viter;
	property_map < BGraph, is_monitor_t>::type Is_Monitor = get(is_monitor_t(), bgraph);
//	if(Is_Monitor[w]){
//		Vertex_d u = v;
//		v = w;
//		w = u;
//	}
	list<Vertex_d> SDNnodes;
	Get_SDNnodes(bgraph, SDNnodes);
	list<Vertex_d>::iterator ib, ie;
	ib = SDNnodes.begin();
	ie = SDNnodes.end();

	//��ȡ��������
	vector<Vertex_d> nodeIndex;
	Get_Vertex_Index(bgraph, nodeIndex);

	Adj_Vertex_i adj_i, adj_ie;
	if(!Is_Monitor[v]){
		list<Vertex_d>::iterator viter;
		for(viter= monitors.begin(); viter != monitors.end(); viter++){
			if(find(ib, ie, *viter) != ie){
				boost::tie(adj_i, adj_ie) = adjacent_vertices(*viter, bgraph);
				for(; adj_i != adj_ie; adj_i++){
					if(nodeIndex[*adj_i] == v){
						pathtemp.push_back(*viter);
						pathtemp.push_back(v);
						paths1.push(pathtemp);
						pathtemp.clear();
					}else{
						Compute_Shortest_Path_From_OneNode(bgraph, nodeIndex[*adj_i], v, pathtemp);
						//�ж�path�ϳ���Ŀ�Ľ���Ƿ���SDN���
						if( Is_SDNnodes_In_Path(SDNnodes, pathtemp) ){
							pathtemp.clear();
							continue;
						}else{
							pathtemp.push_front(*viter);
							paths1.push(pathtemp);
							pathtemp.clear();
						}
					}
				}
			}else{
				Compute_Shortest_Path_From_OneNode(bgraph, *viter, v, pathtemp);
				paths1.push(pathtemp);
				pathtemp.clear();
			}

		}
	}

	if(!Is_Monitor[w]){
		list<Vertex_d>::iterator viter;
		for(viter= monitors.begin(); viter != monitors.end(); viter++){
			if(find(ib, ie, *viter) != ie){
				boost::tie(adj_i, adj_ie) = adjacent_vertices(*viter, bgraph);
				for(; adj_i != adj_ie; adj_i++){
					if(nodeIndex[*adj_i] == w){
						pathtemp.push_back(*viter);
						pathtemp.push_back(w);
						paths2.push(pathtemp);
					}else{
						Compute_Shortest_Path_From_OneNode(bgraph, nodeIndex[*adj_i], w, pathtemp);
						//�ж�path�ϳ���Ŀ�Ľ���Ƿ���SDN���
						if( Is_SDNnodes_In_Path(SDNnodes, pathtemp) ){
							pathtemp.clear();
							continue;
						}else{
							pathtemp.push_front(*viter);
							paths2.push(pathtemp);
							pathtemp.clear();
						}
					}
				}
			}else{
				Compute_Shortest_Path_From_OneNode(bgraph, *viter, w, pathtemp);
				paths2.push(pathtemp);
				pathtemp.clear();
			}
		}
	}

	pathtemp.clear();
	pathtemp.push_back(v);
	pathtemp.push_back(w);
	paths2.push(pathtemp);
}
*/
void Get_Paths_From_Nodes(BGraph& bgraph, list<Vertex_d>& monitors, Vertex_d v, std::queue<PathType>& paths)
{
	PathType pathtemp;
	list<Vertex_d>::iterator viter;
	property_map < BGraph, is_monitor_t>::type Is_Monitor = get(is_monitor_t(), bgraph);
//	if(Is_Monitor[w]){
//		Vertex_d u = v;
//		v = w;
//		w = u;
//	}
	list<Vertex_d> SDNnodes;
	Get_SDNnodes(bgraph, SDNnodes);
	list<Vertex_d>::iterator ib, ie;
	ib = SDNnodes.begin();
	ie = SDNnodes.end();

	//��ȡ��������
	vector<Vertex_d> nodeIndex;
	Get_Vertex_Index(bgraph, nodeIndex);

	Adj_Vertex_i adj_i, adj_ie;
	if(!Is_Monitor[v]){
		list<Vertex_d>::iterator viter;
		for(viter= monitors.begin(); viter != monitors.end(); viter++){
			if(find(ib, ie, *viter) != ie){
				boost::tie(adj_i, adj_ie) = adjacent_vertices(*viter, bgraph);
				for(; adj_i != adj_ie; adj_i++){
					if(nodeIndex[*adj_i] == v){
						pathtemp.push_back(*viter);
						pathtemp.push_back(v);
						paths.push(pathtemp);
						pathtemp.clear();
					}else{
						Compute_Shortest_Path_From_OneNode(bgraph, nodeIndex[*adj_i], v, pathtemp);
						//�ж�path�ϳ���Ŀ�Ľ���Ƿ���SDN���
						if( Is_SDNnodes_In_Path(SDNnodes, pathtemp) ){
							pathtemp.clear();
							continue;
						}else{
							pathtemp.push_front(*viter);
							paths.push(pathtemp);
							pathtemp.clear();
						}
					}
				}
			}else{
				Compute_Shortest_Path_From_OneNode(bgraph, *viter, v, pathtemp);
				paths.push(pathtemp);
				pathtemp.clear();
			}

		}
	}
}


int Find_Linear_Independent_Paths_For_One_Node(BGraph& bgraph, std::queue<PathType>& paths, multimap<double, PathType>& Output_Paths)
{
	int counter;
	counter = 0;
	while(!paths.empty()){
		// cout<<paths.size()<<endl;
		PathType path = paths.front();
		paths.pop();
		boost::numeric::ublas::matrix<long double> v(1, num_edges(bgraph));
		Change_Path_to_Matrix(bgraph,path,v);
		//boost::numeric::ublas::matrix<long double> m_temp =	 boost::numeric::ublas::prod(boost::numeric::ublas::trans(RR),G);
		//R12 = boost::numeric::ublas::prod(m_temp, boost::numeric::ublas::trans(v));

		Compute_The_QMatrix();//�´���
		list<unsigned> record; //�´���
		//�������ѭ��ʽ�´���
		  
		for(unsigned colum = 0; colum < v.size2(); colum++){
			if(v(0,colum) > 0){
				record.push_back(colum); //�´��룬��һ���ڵ�0����һ������Ԫ��
			}
		}
		R12.resize(R.size1(),1);  //�´���
		Compute_R12Matrix(record); //�´���

		R22 = Compute_Squre(v) - Compute_Squre(R12);
		if(R22 < 0){
			R22 = -R22;
		}
		if(R22>0.000000001){
			unsigned k = R.size1()+1;
			
			R22 = sqrt(R22);
			Output_Paths.insert(make_pair(1,path));
			counter++;
			
			if(Output_Paths.size() >= num_edges(bgraph))
				return -1;
			
			G.resize(k,v.size2());
			for(unsigned colum = 0; colum < v.size2(); colum++){
				G(k-1,colum) = v(0, colum);
				if(v(0,colum) > 0){  //����ӵĴ���
					G1[colum].push_back(k-1); //�´��룬��һ���ڵ�K-1����һ������Ԫ��
				}
			}

			R.resize(k,k);
			R(k-1,k-1) = R22; 
			for(unsigned colum = 0; colum < k-1; colum++){
				R(k-1,colum) = 0;
				R(colum,k-1) = R12(colum,0);
			}
			 // boost::numeric::ublas::matrix<long double> temp_matrix =  boost::numeric::ublas::prod(RR,R12);
			boost::numeric::ublas::matrix<long double> temp_matrix(RR.size1(),R12.size2());
			Compute_Prod_Of_RR_and_R12(temp_matrix);
			temp_matrix*=(1/R22);

			RR.resize(k,k);
			RR(k-1,k-1) = 1/R22;
			for(unsigned colum = 0; colum < k-1; colum++){
				RR(k-1,colum) = 0;
				RR(colum,k-1) = -temp_matrix(colum,0);
			}  
		}
	}
	return counter;
}

void Choose_Shortest_Independent_paths(BGraph& bgraph, multimap<double, PathType>& Indepedent_Paths, multimap<double, PathType>& Output_Paths)
{
	boost::numeric::ublas::matrix<long double> v(1, num_edges(bgraph));
    list<Vertex_d> monitors;
	Get_Monitors(bgraph, monitors);
	multimap<double,PathType>::iterator path_iter;
	unsigned k = 1;
	int i = 0;
	int m = 0;
	int sign;
	long double ComputeSqure_R12;
	long double ComputeSqure_v;
	path_iter = Indepedent_Paths.begin();
	int num = 0;
	list<unsigned> templ;   
	G1.assign(num_edges(bgraph), templ);  
	for(path_iter = Indepedent_Paths.begin(); path_iter != Indepedent_Paths.end(); path_iter++){
		if(Output_Paths.size() >= num_edges(bgraph)){
			return;
		}
		num++;
		if(num%50 == 0){   
			cout<<Indepedent_Paths.size()<<endl;
			cout<<num<<endl;
		}
		Change_Path_to_Matrix(bgraph, path_iter->second, v);
		//Out_Put_Matrix(v);
		if(path_iter == Indepedent_Paths.begin()){
			Output_Paths.insert(make_pair(path_iter->first, path_iter->second));
			//��Ӵ��룬��ʱ����paths��dstû�б�ѡΪmonitor
			/********/
			property_map< BGraph, is_monitor_t>::type Is_Monitor = get(is_monitor_t(), bgraph);
			int firstPathTail;
			int firstPathLength;
			firstPathTail = (path_iter->second).back();
			if(Is_Monitor[firstPathTail] == false){
				Is_Monitor[firstPathTail] = true;
				monitors.push_back(firstPathTail);
			}
			/*******/
			G.resize(k, num_edges(bgraph));
			for(unsigned colum = 0; colum < v.size2(); colum++){
				G(0,colum) = v(0, colum);
				//����ӵĴ���
				if(v(0,colum) > 0){
					G1[colum].push_back(0); //�´��룬��һ���ڵ�0����һ������Ԫ��
				}
			}
			long double q = Compute_Squre(v);
			q = sqrt(q);
			RR.resize(k, k);
			RR(0,0) = 1/q;
			R.resize(1, 1);
			R(0,0) = q;
		}else{
			//boost::numeric::ublas::matrix<long double> m_temp =	 boost::numeric::ublas::prod(boost::numeric::ublas::trans(RR),G);
			Compute_The_QMatrix();//�´���
		
			list<unsigned> record; //�´���
			//�������ѭ��ʽ�´���
			for(unsigned colum = 0; colum < v.size2(); colum++){
				if(v(0,colum) > 0){
					record.push_back(colum); //�´��룬��һ���ڵ�0����һ������Ԫ��
				}
			}
			R12.resize(R.size1(), 1); 
			Compute_R12Matrix(record);
			
			//R12 = boost::numeric::ublas::prod(m_temp, boost::numeric::ublas::trans(v));	
			ComputeSqure_v = Compute_Squre(v);
			ComputeSqure_R12 = Compute_Squre(R12);
			R22 = ComputeSqure_v - ComputeSqure_R12;
			if(R22 < 0){
				R22 = -R22;
			}
			if(R22>0.000000001){
				if( path_iter->first < 3){
					R22 = sqrt(R22);
					Output_Paths.insert(make_pair(path_iter->first, path_iter->second));
					k = G.size1()+1;
					G.resize(k, v.size2());
					for(unsigned colum = 0; colum < v.size2(); colum++){
						G(k-1,colum) = v(0, colum);
						//����ӵĴ���
						if(v(0,colum) > 0){
							G1[colum].push_back(k-1); //�´��룬��һ���ڵ�K-1����һ������Ԫ��
						}
					}
					R.resize(k, k);
					R(k-1, k-1) = R22; 
					for(unsigned colum = 0; colum < k-1; colum++){
						R(k-1,colum) = 0;
						R(colum, k-1) = R12(colum, 0);
					}
				
				  // boost::numeric::ublas::matrix<long double> temp_matrix1 =  boost::numeric::ublas::prod(RR,R12);
				  
					boost::numeric::ublas::matrix<long double> temp_matrix(RR.size1(), R12.size2());
					Compute_Prod_Of_RR_and_R12(temp_matrix);
				  // cout<<endl;
				 
					temp_matrix*=(1/R22);
					RR.resize(k, k);
					RR(k-1, k-1) = 1/R22;
					for(unsigned colum = 0; colum < k-1; colum++){
						RR(k-1, colum) = 0;
						RR(colum, k-1) = -temp_matrix(colum, 0);
					}  
				}
				if(path_iter->first >= 3){
					if((path_iter->second).size()!=2){
                        cout<<"严重的错误"<<endl;
						exit(0);
					}
					PathType_iter ppiter = path_iter->second.begin();
					Vertex_d v, w;
					v = *ppiter;
					w = *(++ppiter);
					property_map< BGraph, is_monitor_t>::type Is_Monitor = get(is_monitor_t(), bgraph);
//					if(Is_Monitor[w]){
//						Vertex_d u = v;
//						v = w;
//						w = u;
//					}
					PathType path = path_iter->second;
					std::queue<PathType> paths1, paths2;
					Get_Paths_From_Nodes(bgraph, monitors, v, paths1);
					sign = Find_Linear_Independent_Paths_For_One_Node(bgraph, paths1, Output_Paths);
					if( sign <= -1){
						if(Is_Monitor[v] == false){
							Is_Monitor[v] = true;
							monitors.push_back(v);
						}
						return;
					}else if( sign > 0){
						if(Is_Monitor[v] == false){
							Is_Monitor[v] = true;
							monitors.push_back(v);
						}
					}else{
						cout<<"not moniter"<<endl;
					}
					Get_Paths_From_Nodes(bgraph,monitors, w, paths2);
					sign = Find_Linear_Independent_Paths_For_One_Node(bgraph, paths2, Output_Paths);
					if( sign <= -1){
						if(Is_Monitor[w] == false){
							Is_Monitor[w] = true;
							monitors.push_back(w);
						}
						return;
					}else if(sign >0){
						if(Is_Monitor[w] == false){
							Is_Monitor[w] = true;
							monitors.push_back(w);
						}
					}else{
						cout<<"not moniter"<<endl;
					}
				}
			}
		}
	}
}
void Output_Paths_In_File(BGraph& bgraph, multimap<double, PathType>& Output_Paths, string fileresultname)
{
	list<Vertex_d> monitors;
	list<Vertex_d> SDNnodes;
	list<Vertex_d>::iterator vb, ve;
	multimap<double,PathType>::iterator path_iter;
	ofstream fout;
	PathType pathtmp;
	PathType_iter pb, pe;

	string filename = fileresultname.c_str();
	fout.open(fileresultname.c_str(), ofstream::out);
	Get_Monitors(bgraph, monitors);

	Get_SDNnodes(bgraph, SDNnodes);
	if(SDNnodes.size() == 0){
		monitors.push_front(0);
	}
	//output monitors
	fout<< "monitors:"<< endl;
	ve = monitors.end();
	for(vb = monitors.begin(); vb != ve; vb++){
		fout<< *vb<<" " ;
	}
	fout<< endl;
//	fout<< "monitors数目:"<< endl;
//	fout<< monitors.size()<< endl;
	//output SDNnodes
	fout<< "SDN node:"<< endl;
	ve = SDNnodes.end();
	for(vb = SDNnodes.begin(); vb != ve; vb++){
		fout<< *vb<<" " ;
	}
	fout<< endl;

	//outpath links
	Edge_i eb, ee;
	fout<< "links:" << endl;
	for(boost::tie(eb,ee) = edges(bgraph); eb != ee; eb++){
		fout<< source(*eb, bgraph)<< "   ";
		fout<< target(*eb, bgraph)<< endl;
	}
	//output path
	fout<< "paths:" << endl;
	for(path_iter = Output_Paths.begin(); path_iter != Output_Paths.end(); path_iter++){
		pb = (*path_iter).second.begin();
		pe = (*path_iter).second.end();
		for(; pb != pe; pb++){
			fout<< *pb<< " ";
		}
		fout<< endl;
	}

	fout.close();
}
void Set_Edge_ID(BGraph& bg, vector<vector<int>>& edgeID){
	Edge_i eb, ee;
	boost::tie(eb, ee) = boost::edges(bg);
	int i=0;
	int src, dst;
	for(; eb != ee; eb++, i++){
		src = source(*eb, bg);
		dst = target(*eb,bg);
		edgeID[src][dst] = i;
		edgeID[dst][src] = i;
	}
}
void Check_Linear_Independent(BGraph& bgraph, multimap<double, PathType>& Output_Paths, string matrixname)
{
	//construct edgeID
	int vnum = num_vertices(bgraph);
	int edgeNum = num_edges(bgraph);
	vector<vector<int>> edgeID;
	for(int i=0; i<vnum; i++){
		vector<int> tmpvec(vnum, -1);
		edgeID.push_back(tmpvec);
	}
	Set_Edge_ID(bgraph, edgeID);
	//initialize pathFs
	vector<vector<int>> path_0_1;
	for(int i=0; i<Output_Paths.size(); i++){
		vector<int> tmpvec(edgeNum, 0);
		path_0_1.push_back(tmpvec);
	}
	multimap<double,PathType>::iterator path_iter;
	int k=0; 
	int src, dst;
	int t;
	//construct pathFs
	for(path_iter = Output_Paths.begin(); path_iter != Output_Paths.end(); path_iter++){
		PathType tmpPath = path_iter->second;
		list<int>::iterator vb, ve, vt;
		ve = tmpPath.end();
		vt = tmpPath.begin();
		vt++;
		for(vb = tmpPath.begin(); vt!=ve ;vb++){
			src = *vb;
			dst = *vt;
			vt++;
			t = edgeID[src][dst];
			path_0_1[k][t] = 1;
		}
		k++;
	}
	//output result to file
	ofstream fout;
	fout.open(matrixname.c_str(), ofstream::out);
	fout<< "a=[";
	for(int i=0; i<edgeNum; i++){
		for(int j=0; j<edgeNum; j++){
			fout<< " "<< path_0_1[i][j] ;
		}
		fout<< ";";
		fout<< endl;
	}
	fout<< "]";
	fout.close();

}
void Put_One_Step_Circle_Paths_To_Set(BGraph& bgraph,multimap<double,PathType>& Indepedent_Paths)
{
	list<Vertex_d> SDNnodes;
	list<Vertex_d>::iterator viter1;
	PathType path;
	property_map < BGraph, is_SDN_t>::type Is_SDN = get(is_SDN_t(), bgraph);
	//get SDN nodes
	Get_SDNnodes(bgraph, SDNnodes);
	Adj_Vertex_i adj_i, adj_ie;
	for(viter1 = SDNnodes.begin(); viter1 != SDNnodes.end(); viter1++){
		boost::tie(adj_i, adj_ie) = adjacent_vertices(*viter1, bgraph);
		for(; adj_i != adj_ie; adj_i++){
			unsigned tmpnode = *adj_i;
			if(!Is_SDN[*adj_i]){
				path.push_back(*viter1);
				path.push_back(*adj_i);
				path.push_back(*viter1);
				Indepedent_Paths.insert(make_pair(1, path));
				path.clear();
			}
		}
	}
}
void Main_Process_For_Path_Selection(string filename, string fileresultname)
{
	BGraph bgraph;
	bgraph = ReadTopology(filename);
	Set_SDNnodes_As_Monitors(bgraph);

	property_map < BGraph, is_monitor_t>::type Is_monitor = get(is_monitor_t(), bgraph);

	list<Vertex_d> monitors;
	Get_Monitors(bgraph, monitors);
	multimap<double, PathType> Indepedent_Paths;
	multimap<double, PathType> Output_Paths;
	multimap<double, PathType> Indepedent_Paths_No_Ecmp;

	Set_Graph_Index(bgraph);
	Set_Link_Cost(bgraph);

	Construct_Shortest_Paths_Between_SDNnodes(bgraph, Indepedent_Paths);
	Put_Edge_Paths_To_Set(bgraph, Indepedent_Paths);
	//������ 0 1 0��֮���path����Indepedent_Paths
	Put_One_Step_Circle_Paths_To_Set(bgraph, Indepedent_Paths);
	Choose_Shortest_Independent_paths(bgraph, Indepedent_Paths, Output_Paths);	
//	Check_Linear_Independent(bgraph, Output_Paths, matrixname);
	Output_Paths_In_File(bgraph, Output_Paths, fileresultname);
}

