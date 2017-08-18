#include "qmatrix.h"
#include <boost/foreach.hpp>
#include <functional>
#include <numeric>
#include <algorithm>
#include <cstdlib>

#include <stdlib.h>
#include <time.h>

namespace Clustering
{

	Channel_onenode CaculateInit;	//a mediate voctor for caculating bi

	PoolofNeighborhoods Pool_Neighbors;		//所走node的group的集合
	PoolofNeighborhoods Pool_Groups;		//basin nodes的group的集合, 和c中的顺序一致。
	PoolofNeighborhoods Pool_Groups_copy;
	PoolofNeighborhoods Pool_MotherGroups(100,Neighbors());	//所有nodes的mother group的集合, 包括那些只有一个mother group的node
	PoolofNeighborhoods Pool_CompetitiveMotherGroups;	//overlapping node的mother group的集合
	PoolofNeighborhoods backuppool_original_group(100,Neighbors());//store original group in the phaseII

	robustness_vector ai_vector;
	robustness_vector bi_vector;

	Cluster c;   // a new cluster for storing nodes serving as clusterheads
	Cluster NodesAfterSearch(100,0);
	Cluster overlapping_nodes;
	Cluster real_overlapping_nodes;

	std::vector<int> Nodesvisited(100,0);

	Channel_allnode v_all;
	Channel_allnode v_all_copy;
	Channel_allnode v_all_PR;
	Channel_allnode v_all_PR_copy;

	//phaseI 的遍历遍数
	int Num=0;

	Points PRGroup;


//decide the location of nodes, the number is num_points.
void randomInit	(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR)
	{
//	srand(time(NULL));	//as long as not been excuted all together, this part of program will have new initialization!

    		for (unsigned int j = 0; j < NumCR; j++)		//each node node within the possible scope
		{
		  Point p(GeoDim);
if (NumCR < 100)
			  {
			  for (unsigned int i = 0; i < GeoDim; i++)
				{
					//points distributed in the plane of 100m X 100m
				p[i]=50*(rand() * (1.0) / RAND_MAX);
	//				std::cout <<p(i) << ' ';
				}
			  }
		  else
		  {
		//-------uniformly distribution!--------
				  unsigned int j_0=j % 100;
				  unsigned int j_1=floor(j_0/10);
				  unsigned int j_2=j_0 % 10;


				  //points distributed in the plane of 100m X 100m
/*				  p[0]=j_1*10 + 4 + 2*(rand() * (1.0) / RAND_MAX);
				  p[1]=j_2*10 + 4 + 2*(rand() * (1.0) / RAND_MAX);*/
				  p[0]=j_1*10 + 10*(rand() * (1.0) / RAND_MAX);
				  p[1]=j_2*10 + 10*(rand() * (1.0) / RAND_MAX);
		  }

//		  std::cout<<")"<<"\t";
		ps.push_back(p);
	//		std::cout << std::endl;
		}

// initialize PR nodes' position!
		for (unsigned int j = 0; j < NumPR; j++)		//each node node within the possible scope
		{
		  Point p(GeoDim);
	    //		  std::cout<<"point "<<j<<"=(";
		for (unsigned int i = 0; i < GeoDim; i++)
		    {
				//points distributed in the plane of 100m X 100m
		    p[i]=50*(rand() * (1.0) / RAND_MAX);
	    //				std::cout <<p(i) << ' ';
		    }
	    //		  std::cout<<")"<<"\t";
		PRGroup.push_back(p);
	    //		std::cout << std::endl;
		}
	}










//----------------------------------------------------------------------
//	1. Initialize channels on each node;
//	2. caculate ai and bi on each node;
//	3. Channels are assigned to nodes totlaly randomly, which needs modifying.
//	4. separate this function into channelRandnomInit and caculate_ai_bi, if have time.
//-----------------------------------------------------------------------

//--Aug.10
// channel Initialization on CR and PR.
// parameters:
// - cardinality of PR and CR (ChDim)    1 X parameter
// - radius of PR and CR (RadiusPR, RadiusCR)	2 X parameter
// output v_all 以及每个node的 a_i 和 b_i.

void Clusters::channelRandomInit(unsigned int ChDim, double RadiusPR, double RadiusCR)
{
    //对所有的PR node进行channel初始化
//    srand(time(NULL));
    for(unsigned int k=0; k<PRGroup.size(); k++)
	{
	Channel_onenode v_one(ChDim,1);
//	srand(time(NULL));
	//( value % 100 )  is in the range 0 to 99
	unsigned int x= rand() % 10;	//0-9
//	    double y= (rand() * 10.0) / RAND_MAX;
//	    unsigned int x= (unsigned int)floor(y);
	v_one[x]=0;
	v_all_PR.push_back(v_one);
	}

    //对所有的CR node进行channel初始化
    for(unsigned int l=0; l<NumOfNodes; l++)
	{
	Channel_onenode v_one(ChDim,1);
	for(unsigned int k=0; k<PRGroup.size(); k++)
	    {
	    double dis= sqrt( (_ps[l][0]-PRGroup[k][0]) * (_ps[l][0]-PRGroup[k][0]) + (_ps[l][1]-PRGroup[k][1]) * (_ps[l][1]-PRGroup[k][1]) );
	    if(dis <= RadiusPR)
		{
		for(unsigned int m=0; m<v_one.size(); m++)
		    {
		    v_one[m]=v_one[m]*v_all_PR[k][m];
		    }
		}
	    }
	v_all.push_back(v_one);             //v_all is the whole set of channels on every node!
	v_all_copy=v_all;
	}
	
	
	
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
}

void Clusters::findNeighbors(double radius)
{
	for (unsigned int i=0; i < _ps.size(); i++){
		for (unsigned int j=0; j < _ps.size(); j++)
			{
			if 	((j != i ) && (_dis(i, j) < radius))
				{
				Pool_Neighbors[i].push_back(j);
				}
			}
		}

	
	unsigned int sum_num_neighbors =0;
	for (unsigned int i=0; i<Pool_Neighbors.size(); i++){
		 sum_num_neighbors += Pool_Neighbors[i].size();
		}
	double x = sum_num_neighbors/Pool_Neighbors.size()*1.0;
	std::ofstream myfile4;
	myfile4.open("average_num_neighbor.txt");
	myfile4 << x <<"\n";
	myfile4.close();

//	//output the coordinates of the CR nodes!
//	std::cout<<"\n"<<std::endl;
//	std::cout<<"Node "<<pid<<"'s coordinates are (";
//	for(unsigned int i=0; i<ps[pid])

//	unsigned int w=Pool_Neighbors.size();
//	std::cout<<"there are "<<w<<" neighbors...of course, this number is the same with _ps"<<"\n";

};


// output one single big file: all possible clusters
void Clusters::legitimateClusters(Points & ps, unsigned int numCR, unsigned int seed){

						std::cout << "working..." << "\n";

	// output single node cluster
	std::ofstream myfile3;
	std::string nameEnd1 = SSTR(numCR);
	std::string nameEnd2 = SSTR(seed);
	std::string filename = "potential_clusters_" + nameEnd1 + "_" + nameEnd2 + ".txt";
	myfile3.open(filename.c_str());
	
for(unsigned int clusterSize=1; clusterSize<5; clusterSize++){

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
		if (clusterSize ==2){
								std::cout << "working...2..." << "\n";
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
					for (unsigned int j = 0; j < numCR; j++){
						if(j==i1 || j==i2){
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
					}
				}
			}
			}
		}


	//-----------3--------------------
		if (clusterSize ==3){
								std::cout << "working...3..." << "\n";
		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					Cluster group;
					group.push_back(i1);
					group.push_back(i2);
					group.push_back(i3);
					unsigned int unique = checkUniqueness(group);
					if(unique){
				
					flag = inOneCluster(group);

					if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
						std::cout << "doing 3" << "\n";
							int n = CaculateNumCCC(group);
							if (n){									// there is at least one ICC
								// record the utility matrix
								for (unsigned int j = 0; j < numCR; j++){
									if(j==i1 || j==i2 || j==i3 ){
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
								}// there is common channels within the cluster
							}//geographically close
					}//unique
						}
						}
						}
		}



	//-----------4--------------------
	if (clusterSize ==4){
		std::cout << "working...4..." << "\n";

		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					for (unsigned int i4 = 0; i4 < numCR; i4++){
						Cluster group;
						group.push_back(i1);
						group.push_back(i2);
						group.push_back(i3);
						group.push_back(i4);
						unsigned int unique = checkUniqueness(group);
						if(unique){
							flag = inOneCluster(group);
							if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
							std::cout << "doing 4...." << "\n";
									int n = CaculateNumCCC(group);
									if (n){									// there is at least one ICC
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
		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					for (unsigned int i4 = 0; i4 < numCR; i4++){
						for (unsigned int i5 = 0; i5 < numCR; i5++){
							Cluster group;
							group.push_back(i1);
							group.push_back(i2);
							group.push_back(i3);
							group.push_back(i4);
							group.push_back(i5);
							unsigned int unique = checkUniqueness(group);
							if(unique){
						
							flag = inOneCluster(group);

							if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
													std::cout << "doing 5...." << "\n";
		//						myfile2 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ";
									int n = CaculateNumCCC(group);
									if (n){									// there is at least one ICC
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
										}// there is common channels within the cluster
									}//geographically close
							}//unique
						}
						}
						}
						}
						}
		}




	//-----------6--------------------
		if (clusterSize ==6){
										std::cout << "working...6..." << "\n";
		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					for (unsigned int i4 = 0; i4 < numCR; i4++){
						for (unsigned int i5 = 0; i5 < numCR; i5++){
							for (unsigned int i6 = 0; i6 < numCR; i6++){
							Cluster group;
							group.push_back(i1);
							group.push_back(i2);
							group.push_back(i3);
							group.push_back(i4);
							group.push_back(i5);
							group.push_back(i6);
							unsigned int unique = checkUniqueness(group);
							if(unique){
						
							flag = inOneCluster(group);

							if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
													std::cout << "doing 6..." << "\n";
		//						myfile2 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ";
									int n = CaculateNumCCC(group);
									if (n){									// there is at least one ICC
									// record the utility matrix
									for (unsigned int j = 0; j < numCR; j++){
										if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5 || j==i6){
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
										}// there is common channels within the cluster
									}//geographically close
							}//unique
						}
						}
						}
						}
						}
						}
		}





	//-----------7--------------------
	if (clusterSize ==7){
		std::cout << "working...7..." << "\n";

		int flag =0;
			//if there is one node which can reach all the others, then this group is legitimate
		for (unsigned int i1 = 0; i1 < numCR; i1++){
			for (unsigned int i2 = 0; i2 < numCR; i2++){
				for (unsigned int i3 = 0; i3 < numCR; i3++){
					for (unsigned int i4 = 0; i4 < numCR; i4++){
						for (unsigned int i5 = 0; i5 < numCR; i5++){
							for (unsigned int i6 = 0; i6 < numCR; i6++){
								for (unsigned int i7 = 0; i7 < numCR; i7++){
								Cluster group;
								group.push_back(i1);
								group.push_back(i2);
								group.push_back(i3);
								group.push_back(i4);
								group.push_back(i5);
								group.push_back(i6);
								group.push_back(i7);
								unsigned int unique = checkUniqueness(group);
								if(unique){
									flag = inOneCluster(group);
									if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
															std::cout << "doing 7..." << "\n";
											int n = CaculateNumCCC(group);
											if (n){									// there is at least one ICC
											// record the utility matrix
											for (unsigned int j = 0; j < numCR; j++){
												if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5 || j==i6 || j ==i7){
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
											}// there is common channels within the cluster
										}//geographically close
								}//unique
							}
							}
							}
							}
							}
							}
							}
		}

	}//traversal all cluster sizes
	
	myfile3.close();
}





/*
 * output the possible clusters (size =5) and their number of common channels!
 * legitimate cluster means there is at least one cluster member which can reach the others.
 */
void Clusters::legitimateClusters(Points & ps, unsigned int clusterSize)
{
//	std::ofstream myfile;
//	myfile.open ("legitimateClusters.txt");
	unsigned numCR =0;
	
	if (clusterSize ==2){
	std::ofstream myfile2;
	myfile2.open("numberCC_100_2.txt");
	for (unsigned int i1 = 0; i1 < numCR; i1++){
	for (unsigned int i2 = i1+1; i2 < numCR; i2++){

	//int x=legitimate(ps[i1],ps[i2]);
	// judge whether nodes can be in the same cluster!
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
	for (unsigned int j = 0; j < numCR; j++){
		if(j==i1 || j==i2){
					myfile2 << n;
					}
					else{
						myfile2 << 0;
						}
		if(j!=numCR-1){
					myfile2 << ",";
					}
			else{
				myfile2 << "\n";
				}
		}
	}
}
}
	myfile2.close();
}

//traversal is complished!

	if (clusterSize ==4){
	std::ofstream myfile2;
	myfile2.open("numberCC_200_4.txt");
	int flag =0;
		//if there is one node which can reach all the others, then this group is legitimate
	for (unsigned int i1 = 0; i1 < numCR; i1++){
		for (unsigned int i2 = i1+1; i2 < numCR; i2++){
			for (unsigned int i3 = i2+1; i3 < numCR; i3++){
				for (unsigned int i4 = i3+1; i4 < numCR; i4++){
					Cluster group;
					group.push_back(i1);
					group.push_back(i2);
					group.push_back(i3);
					group.push_back(i4);
					flag = inOneCluster(group);
					if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
//							myfile2 << "close nodes are found!" <<"\n";
							int n = CaculateNumCCC(group);
							if (n){									// there is at least one ICC
//							myfile2 << "there is common channel" <<"\n";
							// record the utility matrix
							for (unsigned int j = 0; j < numCR; j++){
								if(j==i1 || j==i2 || j==i3 || j==i4){
											myfile2 << n;
											}
											else{
												myfile2 << 0;
												}
								if(j!=numCR-1){
											myfile2 << ",";
											}
									else{
										myfile2 << "\n";
										}
								}
							}
							}
					}
				}
				}
				}
	myfile2.close();
	}



	if (clusterSize ==6){
	std::ofstream myfile2;
	myfile2.open("numberCC_300_6.txt");
	int flag =0;
		//if there is one node which can reach all the others, then this group is legitimate
	for (unsigned int i1 = 0; i1 < numCR; i1++){
		for (unsigned int i2 = 0; i2 < numCR; i2++){
			for (unsigned int i3 = 0; i3 < numCR; i3++){
				for (unsigned int i4 = 0; i4 < numCR; i4++){
					for (unsigned int i5 = 0; i5 < numCR; i5++){
						for (unsigned int i6 = 0; i6 < numCR; i6++){
					Cluster group;
					group.push_back(i1);
					group.push_back(i2);
					group.push_back(i3);
					group.push_back(i4);
					group.push_back(i5);
					group.push_back(i6);
					unsigned int unique = checkUniqueness(group);
					if(unique){
						
					flag = inOneCluster(group);
//myfile2 << flag << " ";
					if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
//						myfile2 << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ";
							int n = CaculateNumCCC(group);
							if (n){									// there is at least one ICC
							// record the utility matrix
							for (unsigned int j = 0; j < numCR; j++){
								if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5 || j==i6){
											myfile2 << n;
											}
											else{
												myfile2 << 0;
												}
								if(j!=numCR-1){
											myfile2 << ",";
											}
									else{
										myfile2 << "\n";
										}
								}
							}// there is common channels within the cluster
							}//geographically close
					}//unique
					}
							
//myfile2 << "\n";							
							
					}
				}
				}
				}
				}
		myfile2.close();
	}
	
	
	if (clusterSize ==7){
	std::ofstream myfile2;
	myfile2.open("numberCC_400_7.txt");
	int flag =0;
		//if there is one node which can reach all the others, then this group is legitimate
	for (unsigned int i1 = 0; i1 < numCR; i1++){
		for (unsigned int i2 = i1+1; i2 < numCR; i2++){
			for (unsigned int i3 = i2+1; i3 < numCR; i3++){
				for (unsigned int i4 = i3+1; i4 < numCR; i4++){
					for (unsigned int i5 = i4+1; i5 < numCR; i5++){
						for (unsigned int i6 = i5+1; i6 < numCR; i6++){
							for (unsigned int i7 = i6+1; i7 < numCR; i7++){
					Cluster group;
					group.push_back(i1);
					group.push_back(i2);
					group.push_back(i3);
					group.push_back(i4);
					group.push_back(i5);
					group.push_back(i6);
					group.push_back(i7);
					flag = inOneCluster(group);
					myfile2 << flag << " ";
					if (flag){										// at least one node in the 4 nodes can reach all the other 3 nodes.
							int n = CaculateNumCCC(group);
							if (n){									// there is at least one ICC
							// record the utility matrix
							for (unsigned int j = 0; j < numCR; j++){
								if(j==i1 || j==i2 || j==i3 || j==i4 || j==i5 || j==i6 || j ==i7){
											myfile2 << n;
											}
											else{
												myfile2 << 0;
												}
								if(j!=numCR-1){
											myfile2 << ",";
											}
									else{
										myfile2 << "\n";
										}
								}
							}
							}
							
							myfile2 << "\n";
							
					}
					}
					}
					}
					}
					}
					}
	myfile2.close();
	}
}





//decide all the nodes in group are different!
inline unsigned int Clusters::checkUniqueness(Cluster group){
		unsigned int unique;
		unsigned int initialSize = group.size();
		std::vector<unsigned int>::iterator it = std::unique(group.begin(), group.end());
		group.resize(std::distance(group.begin(),it));

		if (initialSize != group.size())
			{
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

//------------------------------------
//------------------------------------
//这个function是为了计算一个group中common channel的数量
//variable is a group of nodes, the last element in the group is head node!
inline unsigned int Clusters::CaculateNumCCC(Cluster group)
    {
    unsigned int NumCCC;
    if(group.size()>0)
		{
		Channel_onenode ChannelOnHead=v_all_copy[group.back()];//对于需要使用此函数的group，把head node加在最后！
		for (unsigned int i = 0; i < group.size()-1; i++)
			{
			Channel_onenode ChannelOnNeighbor=v_all_copy[group[i]];
			for (unsigned int j=0;j<10;j++)
			  {
			  ChannelOnHead[j]=ChannelOnHead[j]*ChannelOnNeighbor[j];
			  }
			}
		NumCCC=std::accumulate(ChannelOnHead.begin(),ChannelOnHead.end(),0);
		}
    else
		{
		NumCCC=0;
		}
    return NumCCC;
    }
//-----------------------------------

//	=========	Caculate CV of a vector of numbers	=========
    inline double Clusters::CaculateCV(Cluster group){
	//caculate the average value
	double average=(std::accumulate( group.begin(), group.end(), 0 ))/group.size();
	double Standard_deviation;
	double inter1=0;
	double cv;

	for(unsigned int i=0; i<group.size(); i++)
	    {
	    inter1=inter1+pow((group[i]-average),2.0);
	    }
	Standard_deviation=sqrt(inter1/group.size());
	cv=Standard_deviation/average;

	return cv;
	}

//------------------------------------------
//---	For every node, compare its a_i with surrouding nodes,
//---	if	< every of them ---> basin node
//---		<=them, compare b_i of the nodes having identical a_i.
		//this is not sbsolutly right, but is also not wrong
//---		> even one of them --->not basin node

void Clusters::ClusteringPhaseI(){
 
    }



void Clusters::OutputClustersAfterPhaseI()
    {
}







//--------------------------------Clusters::ClusteringPhaseII()---------------------------------------------------------
//---	step1: operation is implemented upon the sets of groups after phaseI. To obtain the sets of monthergroups of
//---	       overlapping nodes, we start to search from the 1st element in 1st group.
//---

void Clusters::ClusteringPhaseII(){

}//phaseII end


void Clusters::clearall()
    {
    Channel_onenode CaculateInit;	//a mediate voctor for caculating bi

    Pool_Neighbors.clear();		//所走node的group的集合
    Pool_Groups.clear();		//basin nodes的group的集合
    Pool_Groups_copy.clear();

    for(unsigned int i=0; i<Pool_MotherGroups.size(); i++)
	{
	Pool_MotherGroups[i].clear();
	}


    Pool_CompetitiveMotherGroups.clear();	//overlapping node的mother group的集合

    for(unsigned int i=0; i<backuppool_original_group.size(); i++)
	{
	backuppool_original_group[i].clear();
	}


    ai_vector.clear();
    bi_vector.clear();

    c.clear();   // a new cluster for storing nodes serving as clusterheads

    for(unsigned int i=0; i<NodesAfterSearch.size(); i++)
	{
	NodesAfterSearch[i]=0;
	}


    overlapping_nodes.clear();
    real_overlapping_nodes.clear();

    for(unsigned int i=0; i<Nodesvisited.size(); i++)
	{
	Nodesvisited[i]=0;
	}


    v_all.clear();
    v_all_copy.clear();
    v_all_PR.clear();
    v_all_PR_copy.clear();

    //phaseI 的遍历遍数
    int Num=0;

    PRGroup.clear();

    for(unsigned int i=0; i<size_distr.size(); i++)
	{
	size_distr[i]=0;
	}


}



/*
//output
         std::cout<<"following nodes have included by more than one cluster:"<<"\n"<<std::endl;
         for(unsigned int i=0;i<Pool_MotherGroups.size();i++)
    	 {
    	 std::cout<<"node "<<overlapping_nodes[i]<<"----> (";
    	 for(unsigned int j=0;j<Pool_MotherGroups[i].size();j++)
    	 {std::cout<<pool_mothergroups[i][j]<<"\t";}
    	 std::cout<<"\n";
    	 }
*/

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