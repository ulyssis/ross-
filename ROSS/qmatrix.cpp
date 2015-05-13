/*
 * Aug 15, in the output potential clusters, there are not duplicates
 */
#include "clusters.h"

namespace Clustering
{


	/*
	 * 2014-5-16
	 * output channel availability to file channel_availability.txt

	std::ofstream myfile;
	myfile.open ("channel_availability.txt");
	for (unsigned int i=0; i< v_all.size(); i++){
		for (unsigned int j=0; j<ChDim; j++){
				if(j!=ChDim-1){
						myfile << v_all[i][j] << ",";
						}
						else{
							myfile << v_all[i][j] << "\n";
								}

				}
			}
	myfile.close();
	 * */






// output one single big file: all possible clusters
void Clusters::legitimateClusters(Points & ps, unsigned int numCR, unsigned int seed, Channel_allnode v_all){

	std::cout << "working..." << "\n";

	// output single node cluster
	std::ofstream myfile3;
	std::string nameEnd1 = SSTR(numCR);
	std::string nameEnd2 = SSTR(seed);
	std::string filename = "potential_clusters_" + nameEnd1 + "_" + nameEnd2 + ".txt";
	myfile3.open(filename.c_str());

for(unsigned int clusterSize=1; clusterSize<4; clusterSize++){

	std::ofstream myfile; // num_channel_each_node_
	std::string nameEnd1 = SSTR(numCR);
	std::string nameEnd2 = SSTR(seed);
	std::string filename = "num_channel_each_node_" + nameEnd1 + "_" + nameEnd2 + ".txt";
	myfile.open(filename.c_str());
// generate number of channels on each node
	for (unsigned int i=0; i< v_all.size(); i++){
		if(i!=v_all.size()-1){
			myfile << std::accumulate(v_all[i].begin(), v_all[i].end(), 0) << " ";
			}
		else{
			myfile << "\n";
			}
				}
	myfile.close();



	//--------------1-------------
	if (clusterSize ==1){
		std::cout << "working...1" << "\n";
		for (unsigned int i = 0; i < numCR; i++){

			for (unsigned int j=0; j<numCR; j++){
					if(j==i){
							myfile3 <<  std::accumulate(v_all[i].begin(), v_all[i].end(), 0);
							}
							else{
								myfile3 << 0;
								}
					if(j!= numCR-1){
					 	myfile3 << ",";
					 	}
					 	else{
						myfile3 << "\n";
					 	}
					}
			}
		}

//---------------2-------------
	if (clusterSize ==2) {
	    std::cout << "working...2..." << "\n";
	    PoolofNeighborhoods recorded; /*store elegite clusters,
					    then the new cluster will be examined based on it.
					    if the cluster has been in 'record', drop.*/
	    for (unsigned int i1 = 0; i1 < numCR; i1++){
	    for (unsigned int i2 = 0; i2 < numCR; i2++){

	    //int x=legitimate(ps[i1],ps[i2]);
	    // judge whether nodes can be in the same cluster!
	    if(i1!=i2){
		    int x=_dis(i1,i2)<10?1:0;

		    if(x){				//if the two nodes can reach each other
			    // calculate the number of ICC
			    Cluster group;

			    group.push_back(i1);
			    group.push_back(i2);
			    int n = CaculateNumCCC(group);
			    // record the legitimate clusters
		    //	myfile << i1 << " " << i2 << ",";
			    // record the utility matrix
			    if(n){
				    std::sort(group.begin(), group.end());
				    int m = checkChosen(recorded, group);
				    if(!m){								// i1, i2 are not already chosen!
					    recorded.push_back(group);
					    for (unsigned int j = 0; j < numCR; j++){
						    if(j==i1 || j==i2) {
							myfile3 << n;
							}
						    else{
							    myfile3 << 0;
							    }
						    if(j!=numCR-1) {
							myfile3 << ",";
							}
						    else {
							myfile3 << "\n";
							}
						}
						}
					}
				}
			}
		    }
		    }
	}


//-----------3--------------------
	if (clusterSize ==3) {
	std::cout << "working...3..." << "\n";
	PoolofNeighborhoods recorded; //vector of vectors
	int flag =0;
		//if there is one node which can reach all the others, then this group is legitimate
	for (unsigned int i1 = 0; i1 < numCR; i1++){
		for (unsigned int i2 = 0; i2 < numCR; i2++){
		    for (unsigned int i3 = 0; i3 < numCR; i3++){
			Cluster group;
			group.push_back(i1);
			group.push_back(i2);
			group.push_back(i3);
			std::sort(group.begin(), group.end());
			unsigned int unique = checkUniqueness(group);
			if(unique) {
			    flag = inOneCluster(group);
			    if (flag) {
				int n = CaculateNumCCC(group);
				if(n){
				    int m = checkChosen(recorded, group);
				    if(!m) {								// i1, i2 are not already chosen!
					recorded.push_back(group);
					for (unsigned int j = 0; j < numCR; j++){
					    if(j==i1 || j==i2 || j==i3) {
						myfile3 << n;
						}
						else {
						myfile3 << 0;
						}
					    if(j!=numCR-1){
						myfile3 << ",";
						}
					    else {
						myfile3 << "\n";
						}
					    }
					}//unique
				    }// there is common channels within the cluster
			    }// in one cluster
			}
		    }
		}
	    }
	}





	//-----------4--------------------
	if (clusterSize ==4) {
		std::cout << "working...4..." << "\n";
		PoolofNeighborhoods recorded; //vector of vectors
		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					for (unsigned int i4 = 0; i4 < numCR; i4++){
						Cluster group;
						PoolofNeighborhoods recorded; //vector of vectors
						group.push_back(i1);
						group.push_back(i2);
						group.push_back(i3);
						group.push_back(i4);
						std::sort(group.begin(), group.end());
						unsigned int unique = checkUniqueness(group);
						if(unique){
							flag = inOneCluster(group);
							if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
							std::cout << "doing 4...." << "\n";
									int n = CaculateNumCCC(group);
									if (n){									// there is at least one ICC
										std::sort(group.begin(), group.end());
										int m = checkChosen(recorded, group);
										if(!m){								// i1, i2 are not already chosen!

				//							myfile2 << "there is common channel" <<"\n";
											// record the utility matrix
											for (unsigned int j = 0; j < numCR; j++){
												if(j==i1 || j==i2 || j==i3 || j==i4){
															myfile3 << n;
															}
															else{
																myfile3 << 0;
																}
												if(j!=numCR-1){
															myfile3 << ",";
															}
													else{
														myfile3 << "\n";
														}
												}
											}
									}// there is common channels within the cluster
								}//geographically close
						}//unique
						}
					}
					}
					}
	}




	//-----------5--------------------
		if (clusterSize ==5){
		std::cout << "working...5..." << "\n";
		PoolofNeighborhoods recorded; //vector of vectors
		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					for (unsigned int i4 = 0; i4 < numCR; i4++){
						for (unsigned int i5 = 0; i5 < numCR; i5++){
							Cluster group;
							PoolofNeighborhoods recorded; //vector of vectors
							group.push_back(i1);
							group.push_back(i2);
							group.push_back(i3);
							group.push_back(i4);
							group.push_back(i5);
							std::sort(group.begin(), group.end());
							unsigned int unique = checkUniqueness(group);
							if(unique){

							flag = inOneCluster(group);

							if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
													std::cout << "doing 5...." << "\n";
		//						myfile2 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ";
									int n = CaculateNumCCC(group);
									if (n){									// there is at least one ICC
										std::sort(group.begin(), group.end());
										int m = checkChosen(recorded, group);
										if(!m){								// i1, i2 are not already chosen!
											// record the utility matrix
											for (unsigned int j = 0; j < numCR; j++){
												if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5){
															myfile3 << n;
															}
															else{
																myfile3 << 0;
																}
												if(j!=numCR-1){
															myfile3 << ",";
															}
													else{
														myfile3 << "\n";
														}
													}
												}
										}// there is common channels within the cluster
									}//geographically close
							}//unique
						}
						}
						}
						}
						}
		}



//	//-----------6--------------------
//		if (clusterSize ==6){
//										std::cout << "working...6..." << "\n";
//		int flag =0;
//			//if there is one node which can reach all the others, then this group is legitimate
//		for (unsigned int i1 = 0; i1 < numCR; i1++){
//			for (unsigned int i2 = 0; i2 < numCR; i2++){
//				for (unsigned int i3 = 0; i3 < numCR; i3++){
//					for (unsigned int i4 = 0; i4 < numCR; i4++){
//						for (unsigned int i5 = 0; i5 < numCR; i5++){
//							for (unsigned int i6 = 0; i6 < numCR; i6++){
//							Cluster group;
//							group.push_back(i1);
//							group.push_back(i2);
//							group.push_back(i3);
//							group.push_back(i4);
//							group.push_back(i5);
//							group.push_back(i6);
//							unsigned int unique = checkUniqueness(group);
//							if(unique){
//
//							flag = inOneCluster(group);
//
//							if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
//													std::cout << "doing 6..." << "\n";
//		//						myfile2 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ";
//									int n = CaculateNumCCC(group);
//									if (n){									// there is at least one ICC
//									// record the utility matrix
//									for (unsigned int j = 0; j < numCR; j++){
//										if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5 || j==i6){
//													myfile3 << n;
//													}
//													else{
//														myfile3 << 0;
//														}
//										if(j!=numCR-1){
//													myfile3 << ",";
//													}
//											else{
//												myfile3 << "\n";
//												}
//											}
//										}// there is common channels within the cluster
//									}//geographically close
//							}//unique
//						}
//						}
//						}
//						}
//						}
//						}
//		}
//
//
//
//
//
//	//-----------7--------------------
//	if (clusterSize ==7){
//		std::cout << "working...7..." << "\n";
//
//		int flag =0;
//			//if there is one node which can reach all the others, then this group is legitimate
//		for (unsigned int i1 = 0; i1 < numCR; i1++){
//			for (unsigned int i2 = 0; i2 < numCR; i2++){
//				for (unsigned int i3 = 0; i3 < numCR; i3++){
//					for (unsigned int i4 = 0; i4 < numCR; i4++){
//						for (unsigned int i5 = 0; i5 < numCR; i5++){
//							for (unsigned int i6 = 0; i6 < numCR; i6++){
//								for (unsigned int i7 = 0; i7 < numCR; i7++){
//								Cluster group;
//								group.push_back(i1);
//								group.push_back(i2);
//								group.push_back(i3);
//								group.push_back(i4);
//								group.push_back(i5);
//								group.push_back(i6);
//								group.push_back(i7);
//								unsigned int unique = checkUniqueness(group);
//								if(unique){
//									flag = inOneCluster(group);
//									if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
//															std::cout << "doing 7..." << "\n";
//											int n = CaculateNumCCC(group);
//											if (n){									// there is at least one ICC
//											// record the utility matrix
//											for (unsigned int j = 0; j < numCR; j++){
//												if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5 || j==i6 || j ==i7){
//															myfile3 << n;
//															}
//															else{
//																myfile3 << 0;
//																}
//												if(j!=numCR-1){
//															myfile3 << ",";
//															}
//													else{
//														myfile3 << "\n";
//														}
//												}
//											}// there is common channels within the cluster
//										}//geographically close
//								}//unique
//							}
//							}
//							}
//							}
//							}
//							}
//							}
//		}

	}//traversal all cluster sizes

	myfile3.close();
}





/*
 * output the possible clusters (size =5) and their number of common channels!
 * legitimate cluster means there is at least one cluster member which can reach the others.
 */



//check whether group has been chosen before!
inline unsigned int Clusters::checkChosen(PoolofNeighborhoods recorded, Cluster group){
	if(recorded.size() >0){
		for (unsigned int i =0; i< recorded.size(); i++){
			if (std::equal(group.begin(), group.end(), recorded[i].begin()) ){
				return 1;
				}

//			unsigned int k =0, flag =1;
//			while(k!=group.size()+1){
//					if(recorded[i][k] != group[k])
//					{
//						return false;
//					}
//					k++;
//				}
			}
		return 0;
		}
}

//decide all the nodes in group are different!
inline unsigned int Clusters::checkUniqueness(Cluster group){
		unsigned int unique;
		unsigned int initialSize = group.size();
		std::vector<unsigned int>::iterator it = std::unique(group.begin(), group.end());
		group.resize(std::distance(group.begin(),it));

		if (initialSize != group.size()) {
			unique = 0;
			}
			else{
			    unique = 1;
			    }
		return unique;
	}

//decide whether this group of nodes are geographically close, i.e. whether exist one node which can reach all the rest nodes
inline unsigned int Clusters::inOneCluster(Cluster group)
{
	for (unsigned int i = 0; i < group.size()-1; i++)
		{
			unsigned int reachAnyOne=1;
			for (unsigned int j = 0; j < group.size()-1; j++)
				{
				// judge whether nodes can be in the same cluster!
				unsigned int x=_dis(group[i],group[j])<10?1:0;
				reachAnyOne = reachAnyOne*x;
				}
				if(reachAnyOne) // node i reaches all the nodes in the group
				{
					return 1;
					}
		}

}


//输出操作符 << 重载
	// single point output
	std::ostream& operator<<(std::ostream& o, const Point& p)
	{
		o << "{ ";
		BOOST_FOREACH(Point::value_type x, p)
		{
			o <<" "<< x;
		}
	        o << " }";
	        o<< std::endl;
		//o << " }<<"\n";
		return o;
	}

	// single cluster output
	std::ostream& operator<<(std::ostream& o, const Cluster& c)
	{
		o << "[ ";
		BOOST_FOREACH(PointId pid, c)
		{
			o << " " << pid;
		}
		o << " ]";

		return o;
	}

	// clusters output
	std::ostream& operator<<(std::ostream& o, const Clusters& cs)
	{
		ClusterId cid = 1;
//----------di
//		std::cout<<"the size of clusters is: "<<sizeof cs<<std::endl;
//----------
		// interate for ||cs._clusters|| times.
		BOOST_FOREACH(Cluster c, cs._clusters)
		{
			o << "c(" << cid++ << ")=";

			BOOST_FOREACH(PointId pid, c)
			{
				o << cs._ps[pid];
			}
			o << std::endl;
		}
		return o;
	}

};
