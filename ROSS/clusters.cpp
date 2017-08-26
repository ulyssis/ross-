/*
 * where is the function legitimateClusters? 21.07.2014
 */
#include "clusters.h"
#include <boost/foreach.hpp>
#include <functional>
#include <numeric>
#include <algorithm>
#include <cstdlib>

#include <stdlib.h>
#include <time.h>

namespace Clustering
{
//number of CR nodes should be set here!!
//int numCR = 100;

	Channel_onenode CaculateInit;	//a mediate vector for caculating bi

	PoolofNeighborhoods Pool_Neighbors;		//所走node的group的集合
	PoolofNeighborhoods Pool_Groups;		//basin nodes的group的集合, 和c中的顺序一致。
	PoolofNeighborhoods Pool_Groups_copy;
	PoolofNeighborhoods Pool_MotherGroups(2000,Neighbors());	//所有nodes的mother group的集合, 包括那些只有一个mother group的node
	PoolofNeighborhoods Pool_CompetitiveMotherGroups;	//overlapping node的mother group的集合
	PoolofNeighborhoods backuppool_original_group(2000,Neighbors());//store original group in the phaseII
	PoolofNeighborhoods Pool_Groups_trueindex(500,Neighbors());

	robustness_vector ai_vector;
	robustness_vector bi_vector;

	Cluster c;   // a new cluster for storing nodes serving as clusterheads
	Cluster NodesAfterSearch(2000,0);
	Cluster overlapping_nodes;
	Cluster real_overlapping_nodes;

	std::vector<int> Nodesvisited(2000,0);

	Channel_allnode v_all;
	Channel_allnode v_all_copy;
	Channel_allnode v_all_PR;
	Channel_allnode v_all_PR_copy;

	//phaseI 的遍历遍数
	int Num=0;

	Points PRGroup;
	unsigned int network_geo_edge;

//decide the location of nodes, the number is num_points.
void randomInit(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR, unsigned int network_geo_size)
	{
//	srand(time(NULL));	//as long as not been excuted all together, this part of program will have new initialization!

    network_geo_edge = network_geo_size;
    		for (unsigned int j = 0; j < NumCR; j++)		//each node node within the possible scope
		{
		  Point p(GeoDim);
//		  //std::cout<<"point "<<j<<"=(";

		  if (NumCR < 100) {
		      /*
		       * the 4 rows and 4 columns
		       */
		      unsigned int m=4;

		      unsigned int j_0=j % 16;
		      unsigned int j_1= floor (j_0/m);
		      unsigned int j_2=j_0 % m;

		      float side = network_geo_size/m;

		      p[0]=j_1* side + side * (rand() * (1.0) / RAND_MAX); // x coordinate
		      p[1]=j_2* side + side * (rand() * (1.0) / RAND_MAX);	// y coordinate

		      for (unsigned int i = 0; i < GeoDim; i++)
			    {
			    p[i]=network_geo_size *(rand() * (1.0) / RAND_MAX);
			    }
		      }

		  else {
		/*
		 * -------quasi uniformly distribution!--------
		 * divide the plane into 10x10 blocks.
		 *
		 */

		  unsigned int j_0=j % 100;
		  unsigned int j_1=floor(j_0/10);
		  unsigned int j_2=j_0 % 10;


		  //points distributed in the plane of 100m X 100m
		  // NOTE: only when numCR is 100*n, (n is integer), the distribution is even!

/*				  p[0]=j_1*10 + 4 + 2*(rand() * (1.0) / RAND_MAX);
		  p[1]=j_2*10 + 4 + 2*(rand() * (1.0) / RAND_MAX);*/
		  float side = network_geo_size/10;

		  p[0]=j_1* side + side * (rand() * (1.0) / RAND_MAX); // x coordinate
		  p[1]=j_2* side + side * (rand() * (1.0) / RAND_MAX);	// y coordinate
		  }


		//--------------------


//		  //std::cout<<")"<<"\t";
		ps.push_back(p);
	//		//std::cout << std::endl;
		}

// initialize PR nodes' position!
		for (unsigned int j = 0; j < NumPR; j++)		//each node node within the possible scope
		{
		  Point p(GeoDim);
	    //		  //std::cout<<"point "<<j<<"=(";
		for (unsigned int i = 0; i < GeoDim; i++)
		    {
				//points distributed in the plane of 100m X 100m
		    p[i]=network_geo_size *(rand() * (1.0) / RAND_MAX);
	    //				//std::cout <<p(i) << ' ';
		    }
	    //		  //std::cout<<")"<<"\t";
		PRGroup.push_back(p);
	    //		//std::cout << std::endl;
		}
	}


// Neighbors is a vector of int
// using this structure, I can formulate into a clustering scheme based on distance!

inline Neighbors Clusters::findNeighbors(PointId pid, double radius)
{
	Neighbors ne;
	//std::cout<<"Node "<<pid<<"'s neighbors are (";
	for (unsigned int j=0; j < _ps.size(); j++)
	{
		if 	((pid != j ) && (_dis(pid, j) < radius))
		{
			ne.push_back(j);
//			//std::cout << "dis(" << pid  << "," << j << ")" << _dis(pid, j) << "<" << threshold << std::endl;
			//std::cout<<j<<"\t";

		}
	}

	//output ne:
	//std::cout <<")"<<std::endl;

	Pool_Neighbors.push_back(ne);

//	//output the coordinates of the CR nodes!
//	//std::cout<<"\n"<<std::endl;
//	//std::cout<<"Node "<<pid<<"'s coordinates are (";
//	for(unsigned int i=0; i<ps[pid])

//	unsigned int w=Pool_Neighbors.size();
//	//std::cout<<"there are "<<w<<" neighbors...of course, this number is the same with _ps"<<"\n";
	return ne;
};

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

Channel_allnode Clusters::channelRandomInit(unsigned int ChDim, double RadiusPR, double RadiusCR)
{
    //对所有的PR node进行channel初始化, each PR works on one channel
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



#ifdef COMPARE_SPECTRUM_WITH_GROUND_TRUTH_METHOD2
    for(unsigned int i=0; i<NumOfNodes; i++ ){
	for(unsigned int k=0; k<ChDim; k++ ){
	    if(!v_all[i][k]){
		float random=(double)rand()/(RAND_MAX);
		v_all[i][k] = (random > FALSE_NEGATIVE_RATE_METHOD2) ? 0:1;
	    }
	}
    }
#endif

//caculate ai_vector and bi_vector after initialing available channels on nodes!
//ai_vector is a vector which stores the values of ai of each node.

    for(unsigned int k=0; k<NumOfNodes; k++)
	{
	int ai=0,bi=0;	  		//对每一个node k，初始化robustness value: ai bi

	Neighbors ne1;
	ne1=findNeighbors(k,RadiusCR);
//output channel infomation!

	for(unsigned int n_ch=0; n_ch<v_all[k].size(); n_ch++)
	    {
	    //std::cout<<v_all[k][n_ch]<<"\t";
	    }
	//std::cout<<"\n";

	CaculateInit=v_all[k];
	Channel_onenode ChannelVectorForbi=v_all[k];
	//begin to caculate ai and bi on node k.
	//calculate ai,
	// ai and bi are spectrum connectivity degree and local connectivity in the paper respectively.
	//---
	for (unsigned int i = 0; i < ne1.size(); i++)
	    {
	    int inner=inner_product( CaculateInit.begin(),CaculateInit.end(),v_all[ne1[i]].begin(), 0);
	    ai=ai+inner;

	    //Calculate bi
	    Channel_onenode ChannelOnNeighbor=v_all[ne1[i]];

	    for (unsigned int j=0; j<ChDim; j++)
		{
		ChannelVectorForbi[j]=ChannelVectorForbi[j]*ChannelOnNeighbor[j];
		}
	    }

	    for(unsigned int h=0; h<ChDim; h++)
		{
		bi=bi+ChannelVectorForbi[h];
		}
	    //std::cout<<"a_"<<k<<"="<<ai<<"\t";
	    //std::cout<<"b_"<<k<<"="<<bi<<"\t";
	    //std::cout<<"\n";
	    //put ai and bi into vector
	    ai_vector.push_back(ai);
	    bi_vector.push_back(bi);


	    /*
	     * calculate average_num_channel_pernode
	     */
	    unsigned int num_pernode;
	    Cluster set_num_pernode;
	    for(unsigned int i=0; i<v_all.size(); i++)
			{
	    	num_pernode = std::accumulate(v_all[i].begin(), v_all[i].end(),0);
	    	set_num_pernode.push_back(num_pernode);
			}
	    average_num_channel_pernode = std::accumulate(set_num_pernode.begin(), set_num_pernode.end(),0)/set_num_pernode.size();
	}
    return v_all;
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
	Channel_onenode ChannelOnHead=v_all[group.back()];//对于需要使用此函数的group，把head node加在最后！
	for (unsigned int i = 0; i < group.size()-1; i++)
	{

	Channel_onenode ChannelOnNeighbor=v_all[group[i]];
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


inline unsigned int Clusters::CaculateNumCCC_v_all_copy(Cluster group)
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

//	=========	Caculate CI of a vector of numbers	=========
    double Clusters::CaculateCV(Cluster group){
	//caculate the average value
    	if(!group.size()){
			double average=(std::accumulate( group.begin(), group.end(), 0 ))/group.size();
			double Standard_deviation;
			double inter1=0;
			double cv;

			for(unsigned int i=0; i<group.size(); i++)
				{
				inter1=inter1+pow((group[i]-average),2.0);
				}
			Standard_deviation=sqrt(inter1/group.size());
			cv=1.96*Standard_deviation/pow(group.size(),0.5);

			return cv;
    	}
    	else
    		return 0;
	}

//	=========	Caculate CI of a vector of numbers	=========

	double Clusters::CaculateCV(doubleVector group){
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
	cv=1.96*Standard_deviation/pow(group.size(),0.5);

	return cv;
	}



//------------------------------------------
//---	For every node, compare its a_i with surrouding nodes,
//---	if	> every of them ---> basin node
//---		== them, compare b_i of the nodes having identical a_i.
		//this is not sbsolutly right, but is also not wrong
//---		< even one of them --->not basin node

#ifdef CLUSTER_SIZE_CONTROL
	void Clusters::ClusteringPhaseI(unsigned int optimalSize)
#endif
#ifndef CLUSTER_SIZE_CONTROL
	void Clusters::ClusteringPhaseI()
#endif
	{
#ifndef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME

    unsigned int number_visited =0;
    //-----phaseI sp1
    Cluster tem_group;


    while(number_visited !=NumOfNodes)
		{
	std::cout<< "xxxxx" << "\n";
		for(unsigned int i=0; i<NumOfNodes; i++)	//traversal nodes for clusterhead!
			{
			//test
			unsigned int test272=Nodesvisited[i];
			if(!Nodesvisited[i])		//if node i is not included, then...
				{
				//judge whether the channel of this node is all 0 !!
				unsigned int BeaconOfAllZero = std::accumulate(v_all[i].begin(),v_all[i].end(),0);
				if(BeaconOfAllZero>0)
					{
						//x, y are just flag denoting certain "qualifications"
					unsigned int x=0,y=0,z=0;
					Cluster neighboring_identical_ai;
					Cluster Neighbor=Pool_Neighbors[i];	//node i's neighborhood, members are represented by index No.

					for(unsigned int j=0;j<Neighbor.size();j++)	//get into neighborhood of node i, to check whether it is clusterhead.
											//judge formation of clustering via comparison of ai an bi.
						{
						unsigned int _j=Neighbor[j];
						//判断相等用==
						x=x+((ai_vector[i]>ai_vector[_j])+(ai_vector[i]==Neighbor[_j]));
						y=y+(ai_vector[i]>ai_vector[_j]);
						}

					/*
					 * First situation
					 * a_i is strictly smaller than all its neighbors, node i is cluster head!
					 */
					if(!x)
						{
						//std::cout<<"node "<<i<<" is basin node and clusterhead"<<"\n";

						//===	add	===
						//	phaseI sp1:	shrink the new cluster which has no ccc, to make every cluster has
						//	at least one CCC;
						//
						Cluster changingmem = Neighbor;	//node i's neighbors
						changingmem.push_back(i);	//this is a complete neighborhood with head node!, head in at the last place!
						unsigned int beacon = CaculateNumCCC(changingmem);

						/*
						 * shrink
						 */
	#ifdef CLUSTER_SIZE_CONTROL
						while(!beacon || changingmem.size() > optimalSize+1) // changingmem needs to satisfy two conditions: 1. there is common channel, 2. clusters size smaller or equal to optimalSize+1
	#endif
	#ifndef CLUSTER_SIZE_CONTROL
						while(!beacon)
	#endif
							{
							Cluster set_cc;
							//get the set of inner_product between node i and neighbors
							for(unsigned int j=0; j< changingmem.size()-1; j++)
								{
								unsigned int cc=inner_product(v_all_copy[i].begin(),v_all_copy[i].end(),v_all_copy[changingmem[j]].begin(), 0);
								set_cc.push_back(cc);
								}

							unsigned int min_cc = *min_element(set_cc.begin(),set_cc.end()); //min_element returns the smallest element

							Cluster set_cc_min; // to contain the smallest elements which may be more than one.

							//find the set of smallest inner_product
							for(unsigned int k=0; k < set_cc.size(); k++)
								{
								if(min_cc==set_cc[k])
									{
									set_cc_min.push_back(k);	//store the position index of element in set_cc/changingmem
									}
								}

							//if there are more than one neighbors having the smallest inner_product,
							//consider the change of #CCC, kick off the one which can bring the biggest increase of #CCC
							if(set_cc_min.size()>1)
								{
								Cluster set_cc_min_delta;
								for(unsigned int l=0; l < set_cc_min.size(); l++)
									{
									tem_group=changingmem;
									std::vector<unsigned int>::iterator where = std::find(tem_group.begin(), tem_group.end(), changingmem[set_cc_min[l]]);
									tem_group.erase(where);
									unsigned int delta=CaculateNumCCC(tem_group);
									set_cc_min_delta.push_back(delta);
									}
								unsigned int maxdelta=*max_element(set_cc_min_delta.begin(),set_cc_min_delta.end());

								unsigned int m=0;
								while(set_cc_min_delta[m]!=maxdelta)
									{
									m++;
									}			//m 也是 是set_cc_min_delta中值最大的元素的pisition index， 同时也是set_cc_min
								std::vector<unsigned int>::iterator where = std::find(changingmem.begin(), changingmem.end(), changingmem[set_cc_min[m]]);
								changingmem.erase(where);
								}

							else
								{
								std::vector<unsigned int>::iterator where = std::find(changingmem.begin(), changingmem.end(), changingmem[set_cc_min[0]]);
								changingmem.erase(where);
								}

							beacon=CaculateNumCCC(changingmem);
							}

						//对于Pool_Groups中的每一个成员，不能包含已知的cluster heads
						Cluster changingmem_copy=changingmem;
						for(unsigned int n=0; n<changingmem_copy.size(); n++)
							{
							for(unsigned int p=0; p<c.size();p++)
								{
								if(changingmem_copy[n]==c[p])
									{
									std::vector<unsigned int>::iterator where =std::find(changingmem.begin(), changingmem.end(), changingmem_copy[n]);
									changingmem.erase(where);
									}
								}
							}

						// change node's ai, which is in changingmem.
	//					setAi(changingmem);


						Pool_Groups.push_back(changingmem);
						//


						for(unsigned int j=0;j<changingmem.size();j++) //1. improve ai compulsorily; 2. label nodes of i's neighbood tobe 'cluster head impossible'
							{
							Nodesvisited[changingmem[j]]=1;
							}
						c.push_back(i);			//put this node(head) into set 'c'.
						}

					/*
					 * the second situation to judge a basin node
					 * a_i equal to certain neighboring nodes.
					 */
					else if(!y) // x=1, y =0, there are equalities
						{
						for(unsigned int j=0;j<Neighbor.size();j++)
							{
							unsigned int t=Neighbor[j];	//t is still the index No.
							if (ai_vector[i]==ai_vector[t])
								{
								neighboring_identical_ai.push_back(t);	//find the set of nodes with identical ai,把有相同ai的点找出来！
								}
							}
					   //找出了node i周围有相同ai的点！ 但是还得看，这些点是不是具备了成为head的可能，也就是“周围点的a_i没有比它的a_i大的”！！

						Cluster basinpossible_neighboring_identical_ai;
						unsigned int yy=0;
						for(unsigned int k=0; k<neighboring_identical_ai.size(); k++)
							{
							unsigned int a=neighboring_identical_ai[k];
							for(unsigned int l=0; l < Pool_Neighbors[a].size(); l++)
								{
								unsigned int b=Pool_Neighbors[a][l];
								//xx is impossible to equal to 0!
								//xx=xx+((ai_vector[b]>ai_vector[a])+(ai_vector[b]==Neighbor[a]));
								yy=yy+(ai_vector[b] < ai_vector[a]);
								}
							if(!yy)	//找到了那些临近的，具有相同ai并且有可能成为head的neighbors
								{
								basinpossible_neighboring_identical_ai.push_back(a);
								}
							}

						/*
						 * if basinpossible_neighboring_identical_ai is not empty, then check z, if empty, then z=0.
						 */
						for(unsigned int m=0; m<basinpossible_neighboring_identical_ai.size(); m++)
							{
							unsigned int t=basinpossible_neighboring_identical_ai[m];
							z=z+(bi_vector[i]<bi_vector[t]);
							}

						if(!z)	//z=0 means node bi_i>=bi_neighbors
							{
							//std::cout<<"node "<<i<<" is cluster head although certain neighboring nodes having identical ai"<<"\n";

							//===	add	===
							//	phaseI sp1:	decrease the new cluster which has no ccc.
							Cluster changingmem = Neighbor;
							changingmem.push_back(i);	//this is a complete neighorhood with head node!
							unsigned int beacon = CaculateNumCCC(changingmem);

							/*
							 * shrink
							 */
	#ifdef CLUSTER_SIZE_CONTROL
						while(!beacon || changingmem.size() > optimalSize+1) // changingmem needs to satisfy two conditions: 1. there is common channel, 2. clusters size smaller or equal to optimalSize+1
	#endif
	#ifndef CLUSTER_SIZE_CONTROL
						while(!beacon)
	#endif
								{
								Cluster set_cc;
								for(unsigned int j=0; j< changingmem.size()-1; j++)
									{
									unsigned int cc=inner_product(v_all_copy[i].begin(),v_all_copy[i].end(),v_all_copy[changingmem[j]].begin(), 0);
									set_cc.push_back(cc);
									}
				/*				for(unsigned int k=0; k<set_cc.size();k++)
									{
									//std::cout<<set_cc[k]<<"\t";
									}
								 //std::cout<<std::endl;
				*/
								unsigned int min_cc=*min_element(set_cc.begin(),set_cc.end());

								Cluster set_cc_min;
								for(unsigned int k=0; k < set_cc.size(); k++)
									{
									if (min_cc == set_cc[k])
										{
										set_cc_min.push_back(k);	//store the position index of element in set_cc/changingmem
										}
									}

								if(set_cc_min.size()>1)
									{
									Cluster set_cc_min_delta;

									for(unsigned int l=0; l < set_cc_min.size(); l++)
									{
									tem_group=changingmem;
									std::vector<unsigned int>::iterator where = std::find(tem_group.begin(), tem_group.end(), changingmem[set_cc_min[l]]);
									tem_group.erase(where);
									unsigned int delta=CaculateNumCCC(tem_group);
									set_cc_min_delta.push_back(delta);
									}
									unsigned int maxdelta=*max_element(set_cc_min_delta.begin(),set_cc_min_delta.end());

										unsigned int m=0;
										while(set_cc_min_delta[m]!=maxdelta)
										{
										m++;
										}
										std::vector<unsigned int>::iterator where = std::find(changingmem.begin(), changingmem.end(), changingmem[set_cc_min[m]]);
										changingmem.erase(where);
									}
								else
									{
									std::vector<unsigned int>::iterator where = std::find(changingmem.begin(), changingmem.end(), changingmem[set_cc_min[0]]);
									changingmem.erase(where);
									}
								beacon=CaculateNumCCC(changingmem);
								}

							//对于Pool_Groups中的每一个成员，不能包含已知的cluster heads
							Cluster changingmem_copy=changingmem;
							for(unsigned int n=0; n<changingmem_copy.size(); n++)
								{
								for(unsigned int p=0; p<c.size();p++)
									{
									if(changingmem_copy[n]==c[p])
										{
										std::vector<unsigned int>::iterator where =std::find(changingmem.begin(), changingmem.end(), changingmem_copy[n]);
										changingmem.erase(where);
										}
									}
								}

	//						setAi(changingmem);

							Pool_Groups.push_back(changingmem);

							for(unsigned int j=0;j<changingmem.size();j++)
								{
								unsigned s=changingmem[j];
								Nodesvisited[s]=1;
								}
							c.push_back(i);
							}

				//		else	//i的nrighboring node j， 其a和点i的相同，b也比点i的大， 但是因为不一定能够满足head的条件，所以不认定其为cluster head，
				//		    {	//对于j的判断，留在下一个循环中进行。而对i的判断也将停止，并认定点i也暂时不是cluster head。

						}//Identical ai exist! finding out cluster head！



					}//avoid all channels are 0!
				else if (!BeaconOfAllZero)
					{
					Nodesvisited[i]=1;
					}
				}//avoid repeated access!
			}// for each node(or each node with neighborhood of it)
		number_visited = CheckVisitComplete();
		//std::cout<<std::endl;
		}
    //std::cout<<"The number of traversal is "<< Num <<" times!"<<std::endl;
//---record the real groups after phase I
    Pool_Groups_copy=Pool_Groups;
#endif
}

/*
 * return the number of visited nodes!
 * and
 * set ai as big number;
 */
unsigned int Clusters::CheckVisitComplete()
    {
    //phaseI 的遍历遍数
	Num++;
	int n=0;
	for(unsigned int i=0;i< Nodesvisited.size();i++)
	    {
	    if(Nodesvisited[i]==1)
		{
		ai_vector[i]=100000;
		}
	    n=n+Nodesvisited[i];
	    }
	return n;
    }

void Clusters::setAi(robustness_vector robustness_vector)
	{
	for(unsigned int i=0;i< robustness_vector.size();i++)
		{
		robustness_vector[i]=0;
		}
	}


void Clusters::OutputClustersAfterPhaseI()
    {
    int NumOfAllNodes=0;
	//std::cout<<"there are "<<Pool_Groups.size()<<" group(s)"<<"\n";
//	//std::cout<<"there are "<<c.size()<<" group(s)"<<"\n";
	for(unsigned int i=0;i<Pool_Groups.size();i++)
	    {
	    //std::cout<<"Group("<<c[i]<<")= "<<"\t";
	    for(unsigned int j=0;j<Pool_Groups[i].size();j++)
		{
		//std::cout<< Pool_Groups[i][j] <<"\t";
		NumOfAllNodes++;
		}
	    //std::cout<<"\n";
	    }
//output the total number of nodes in all groups, including ovelapping nodes!
	std::cout<<"the number of all nodes involved is "<<NumOfAllNodes<<std::endl;

//	//output all element of ai_vector after phaseI
//	for(unsigned int l=0;l<100;l++)
//	    {
//	    std::cout<<"ai_vector["<<l<<"]="<<ai_vector[l]<<"\t"<<std::endl;
//	    }



//	//output nodes within groups, make sure there is no node which is not included into a cluster1
//	std::vector<unsigned int> NodesInTheGroups;
//	for(unsigned int i=0;i<Pool_Groups.size();i++)
//	    {
//	    for(unsigned int j=0;j<Pool_Groups[i].size();j++)
//		for(unsigned int k=0;k<NodesInTheGroups.size();k++)
//		    {
//		    int test=(Pool_Groups[i][j]==NodesInTheGroups[k]);
//		    }
//		    NodesInTheGroups.push_back(j);
//	    }

	//output sum of ai_vector, if sum_ai should be 100000, which means all nodes are visited!
	int sum_ai=0;
	for(unsigned int i=0;i<ai_vector.size();i++)
	    {
	    sum_ai=sum_ai+ai_vector[i];
	    }
	std::cout<<"the sum of ai on all nodes after phaseI is "<<sum_ai<<std::endl;

	// output the formed cluster after phase I
	    std::cout<<"formed clusters after phase I is: "<<"\n";

	    for(unsigned int i=0;i<Pool_Groups.size();i++)
	    	    {
	    	    for(unsigned int j=0; j<Pool_Groups[i].size(); j++)
	    		{
	    		std::cout<< Pool_Groups[i][j] <<"\t";
	    		}
	    	    std::cout<<"\n";
	    	    }
}







//--------------------------------Clusters::ClusteringPhaseII()---------------------------------------------------------
//---	step1: operation is implemented upon the sets of groups after phaseI. To obtain the sets of monthergroups of
//---	       overlapping nodes, we start to search from the 1st element in 1st group.
//---
#ifndef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME
void Clusters::ClusteringPhaseII()
#endif

#ifdef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME
    void Clusters::ClusteringPhaseII(unsigned int seed)
#endif
    {
    //std::cout<<"following node is included by more than one cluster:"<<"\n"<<std::endl;
#ifndef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME
    for(unsigned int h=0;h<Pool_Groups.size();h++)
	{
	for(unsigned int i=0;i<Pool_Groups[h].size();i++)	//对于每个group里面的node：Pool_Groups[h][i]
	    {
	    unsigned int index=Pool_Groups[h][i];	//index is index of node.
//	    int beacon=0;
//	    int a=NodesAfterSearch[index];
	    if(!NodesAfterSearch[index])
		{
	    //---------
	    for(unsigned int j=0;j<Pool_Groups.size();j++)	//遍历每个group里面的node：Pool_Groups[j][k]
		{
		for(unsigned int k=0;k<Pool_Groups[j].size();k++)
		    {
			if(Pool_Groups[j][k]==index)
			{
			//this set contains the nodes which belong to more than one group
			overlapping_nodes.push_back(index);
			//put group index 'j' into node 'Pool_Groups[h][i]''s mothergroup, if node i find
			//can be simplified by using back()

//			int RealIndexOfGroup=Pool_Groups[j].size()-1;
			Pool_MotherGroups[index].push_back(Pool_Groups[j].back());	//collect the index of groups containing the node, notice that, element of Pool_MotherGroups is corresponding to index of node!
			//该node在查找完之后被labeled！
			NodesAfterSearch[index]=1;
			}
		    }
		}

	    if(Pool_MotherGroups[index].size()>1)
		{
		Pool_CompetitiveMotherGroups.push_back(Pool_MotherGroups[index]);	//后面没有用到
		real_overlapping_nodes.push_back(index);				//step2 有用到

	    //output the belongings of overlapping nodes
		//std::cout<<"node "<<index<<" --belong to--> (";
		for(unsigned int j=0;j<Pool_MotherGroups[index].size();j++)
		    {
		    //std::cout<<Pool_MotherGroups[index][j]<<"\t";
		    }
		//std::cout<<")"<<"\t"<<"belongs to "<<Pool_MotherGroups[index].size()<<" groups"<<"\n";;
		}

		}
	    }
	 }

    //number of overlapping nodes
    numOverlappingNode = real_overlapping_nodes.size();


    //Pool_MotherGroups[i] contains the index of cluster heads which claim i
    unsigned int totalNumClaimingCluster=0;
    for(unsigned int i=0; i<Pool_MotherGroups.size(); i++){
    	if (Pool_MotherGroups[i].size() >1){
    		totalNumClaimingCluster +=Pool_MotherGroups[i].size();
    	}
    }
    //average number of claiming clusters.
    average_ClaimingCluster = totalNumClaimingCluster/(numOverlappingNode*1.0);


/*
    //std::cout<<"the number of nodes within all the groups(including repeated nodes) are "<<overlapping_nodes.size()<<"\n";
    //std::cout<<"the number of overlapping nodes are "<<real_overlapping_nodes.size()<<"\n";
    //std::cout<<"the number of possible mother groups are "<<Pool_CompetitiveMotherGroups.size()<<"\n";
    //std::cout<<"the number of primitive possible mother groups are "<<Pool_MotherGroups.size()<<"\n";
*/

#ifndef GREEDY
//-----------------------------------------------------------
//---	step2 of clarification begin from here---------------
//---	for each overlapping node, when it comes to the time to decide which cluster to go, relevant clusters will shrink with the
//---	departure of the overlapping node! then the delta_bi can be calculated by:
//---							delta_bi=bi_new - bi_vector[index_of_competitive_mothergroup]
//---	if delta_bi equal, compare the number of common channels between overlapping nodes and clustered! the bigger will be chosen!
//---	if still equal, size of possible mothers cluster will be considered!
//---
//---
//---
    for(unsigned int h=0;h<real_overlapping_nodes.size();h++)//遍历real_overlapping_nodes中每个元素，即每个overlapping的node
	{

	unsigned int index_of_overlapping_node = real_overlapping_nodes[h];	//get authentic index for overlapping node!
	//	begin to exame and compare every possible mother groups
	//	define set of delta_bi
	Cluster set_delta_bi;
	//	define set of num_hn
	Cluster num_hn;

	int delta_bi;
	unsigned int bi_new;
	unsigned int bi_old;
	Cluster set_tem;
	//------
	//start to compare possible groups!
//output to test!
	//std::cout<<"As to the overlapping node of "<<index_of_overlapping_node<<"\n";
	//std::cout<<"the size of its mother group is "<<"\n";
	for(unsigned int i=0;i<Pool_MotherGroups[index_of_overlapping_node].size();i++)
	    {
	    unsigned int index_of_competitive_mothergroup = Pool_MotherGroups[index_of_overlapping_node][i];	//get index for possible mother group!
	    //backup把原来的cluster放入backuppool_original_group中第index_of_competitive_mothergroup个元素！
	    //此后，Pool_Groups不在是phaseI结束后的情况！！！

	    //查找 index_of_competitive_mothergroup 在 c 中的位置！这个位置也就是index_of_competitive_mothergroup 在Pool_Groups中的位置！
	    unsigned int counter=0;
	    unsigned int ii=0;
	    while(c[ii]!=index_of_competitive_mothergroup)
		{
		counter++;
		ii++;
		}
	    backuppool_original_group[index_of_competitive_mothergroup]=Pool_Groups[counter];
	    bi_old = CaculateNumCCC(Pool_Groups[counter]);
//output to test!
	    //std::cout<<Pool_Groups[counter].size()<<"\t";
	    //delete overlapping nodes
	    std::vector<unsigned int>::iterator where = std::find(Pool_Groups[counter].begin(), Pool_Groups[counter].end(), index_of_overlapping_node);
	    Pool_Groups[counter].erase(where);	//new group, after deleting an overlapping node

//	    //delete clusterhead(should have not been deleted, but in actual operation, clusterhead will be inner_product for 2 times), to save caculation time!
//	    std::vector<unsigned int>::iterator where1 = std::find(Pool_Groups[counter].begin(), Pool_Groups[counter].end(), index_of_competitive_mothergroup);
//	    Pool_Groups[counter].erase(where1);

	    //caculate some metrics of relavant clusters(after shrinking):
	    //two metrics:
	    //		1. delta_bi: num of common channels
	    //		2. num_hn: num of common channels between  index_of_overlapping_node and index_of_competitive_mothergroup;

	    bi_new = CaculateNumCCC(Pool_Groups[counter]);

	    delta_bi = bi_new - bi_old;
	    set_delta_bi.push_back(delta_bi);	//set_delta_bi中元素的个数和顺序同Pool_MotherGroups[index_of_overlapping_node]中的是一致的，只是具体的index不一样，可资利用！
	    //if there is smallest delta_bi exists, keep it inside the group, and recover group membership, then end!
	    //if not, compare num_hn, and find the biggest,

	    //cluster recover!
	    Pool_Groups[counter]=backuppool_original_group[index_of_competitive_mothergroup];
	    }



//output to test!
	//std::cout<<std::endl;
	//std::cout<<"the size of set_delta_bi of node "<<"index_of_overlapping_node is "<<set_delta_bi.size()<<","<<"\n";
	//std::cout<<"thay are: "<<"\n";
//	for (unsigned int aa=0; aa < set_delta_bi.size(); aa++)
//	    {
//	    //std::cout<<set_delta_bi[aa]<<"\t";
////	    unsigned int test=set_delta_bi[aa];
//	    }
	//std::cout<<std::endl;


	//here:	begin to compare metrics
	//find the smallest delta_bi
	Cluster couter;
	Cluster couter2;
	unsigned int stone=set_delta_bi[0];
	for(unsigned int m=0;m<set_delta_bi.size();m++)
	    {
	    if(stone>set_delta_bi[m])
		{
		stone=set_delta_bi[m];
		}
	    }	//stone is the smallest element in set_delta_bi.

	for(unsigned int m=0; m < set_delta_bi.size(); m++)
	    {
	    if (stone == set_delta_bi[m])
		{
		couter.push_back(m);//record relative position index in set_delta_bi, which are about 'mother groups' which have equal delta.
		}
	    }

	 unsigned int index_realmother;


	 //if true, there are identical delta_bi, with same small value.
	if(couter.size()>1)
	    {
	    for(unsigned int p=0;p<couter.size();p++)
		{
		//set_tem is a group of index of possible mother clusters' heads.
		//couter[p]就是在set_delta_bi中的次序, for one overlapping node, this pool stores the index of mother groups, which are index of nodes in the same time!
		set_tem.push_back(Pool_MotherGroups[index_of_overlapping_node][couter[p]]);
		}


	    Channel_onenode CaculateInit=v_all_copy[index_of_overlapping_node];
	    Cluster set_inner_product;	//get delta_bi
	    //caculate number of common channel between index_of_overlapping_node and relavent cluster heads.find the biggest!
	    for(unsigned int r = 0; r < set_tem.size(); r++)
		{
		 int inner=inner_product(CaculateInit.begin(),CaculateInit.end(),v_all_copy[set_tem[r]].begin(), 0);
		 set_inner_product.push_back(inner);
		}

	    unsigned int stone2=set_inner_product[0];
	    for(unsigned int q=0; q < set_inner_product.size(); q++)
		{
		 if(stone2 < set_inner_product[q])	//stone2 是诸多#cc between overlapping node and possble mother clusters 中的最大值, biggest value of inner
		     {
		     stone2 = set_inner_product[q];
		     }
		}
//	    for(unsigned int q=0; q < set_inner_product.size(); q++)
//		{
//		//std::cout<<"the two number of #common channels are: "<<set_inner_product[q]<<"\t"<<"\n";
//		}
	    for(unsigned int r=0; r < set_inner_product.size(); r++)
		{
		if (stone2 == set_inner_product[r])
		    {
		    //***********
		    //couter2 stores index of mother cluster heads, which have the same hn and delta_bi
		    //***********
		    couter2.push_back(set_tem[r]);
		    }
		}

	    unsigned int smallestsize = Pool_Groups[couter2[0]].size();


	    //there are more than one mother cluster having the identical smallest delta_bi and biggest hn
	    if(couter2.size()>1)
		{
		//find the one with smallest size out of the possible mother clusters which have identical delta_bi and hn
		for(unsigned int s=0; s<couter2.size(); s++)
		    {
		    if(smallestsize < Pool_Groups[couter2[s]].size())
			{
			smallestsize = Pool_Groups[couter2[s]].size();
			}
		    }

		for(unsigned int t=0; t<couter2.size(); t++)
		    {
		    if(Pool_Groups[couter2[t]].size()==smallestsize)
			{
			index_realmother=couter2[t];
			//std::cout<<"Justice about size was done!"<<"\n";
			}
		    }
		}
	    else
		{
		index_realmother=couter2[0];
		}

	    //找到delta_bi最小的，就把overlapping node留在那里.因为目前所有的group都是原处的规模，所以要把其他group里的overlapping ndoe删除！
	    for(unsigned int j=0; j<Pool_MotherGroups[index_of_overlapping_node].size();j++)
		{
		unsigned int index_cc=Pool_MotherGroups[index_of_overlapping_node][j];
		if(index_cc!=index_realmother)
		    {
		    unsigned int ii_3=0;
		    while(c[ii_3]!=index_cc)
			{
			ii_3++;
			}
		    std::vector<unsigned int>::iterator where = std::find(Pool_Groups[ii_3].begin(), Pool_Groups[ii_3].end(), index_of_overlapping_node);
		    Pool_Groups[ii_3].erase(where);
		    }
		}
	    }
	else
	    {
	    // stone record the smallest delta_bi, and including the location of the delta_bi and the group!
	    unsigned int bb=0;
	    while(set_delta_bi[bb]!=stone)
		{
		bb++;
		}

	    //group(index_bb) remain! Pool_MotherGroups[index_of_overlapping_node]中其他的cluster删点！
	    index_realmother = Pool_MotherGroups[index_of_overlapping_node][bb];

	    for(unsigned int j=0; j<Pool_MotherGroups[index_of_overlapping_node].size();j++)
		{
		unsigned int index_cc=Pool_MotherGroups[index_of_overlapping_node][j];
		if(index_cc!=index_realmother)
		    {
		    unsigned int ii_3=0;
		    while(c[ii_3]!=index_cc)
			{
			ii_3++;
			}
		    std::vector<unsigned int>::iterator where = std::find(Pool_Groups[ii_3].begin(), Pool_Groups[ii_3].end(), index_of_overlapping_node);
		    Pool_Groups[ii_3].erase(where);
		    }
		}
	    }
	}//step2 of phaseII
#endif

/*
 * greedy scheme, converge in the end!
 */
#ifdef GREEDY
    //-----------------------------------------------------------
    //---	step2 of clarification begin from here---------------
    //---	LOCAL FAST ALGORITHM:
    //---	for each overlapping node, when it comes to the time to decide which cluster to go, relavant clusters will shrink with the
    //---	departure of the overlapping node! then the delta_bi can be caculated by:
    //---							delta_bi=bi_new - bi_vector[index_of_competitive_mothergroup]
    //---	if delta_bi equal, compare the number of common channels between overlapping nodes and clusterhead! the bigger will be chosen!
    //---	if still equal, size of possible mothers cluster will be considered!
    //---
    //---
    //---	LOCAL GREEDY ALGORITHM:
    //---	rewrite the code according to paper!

        //this vector contains the correspondings between debatable node and its temporal cluster
        Cluster tem_mum(1000);
        Cluster set_delta_bi;
        Cluster set_delta_icc;//decrements on relavant clusters in one iteration.
        Cluster set_same_delta_icc;//stores the index of C\subset S_i, which have the same increment of ICCs
        Cluster set_clustersize;


    	//this substep1 put debatable node into only one group.
    	for(unsigned int h=0;h<real_overlapping_nodes.size();h++)//遍历real_overlapping_nodes中每个元素，即每个overlapping的node
    	    {
    	    unsigned int index_of_overlapping_node = real_overlapping_nodes[h];	//get authentic index for overlapping node!
    	    //	begin to exame and compare every possible mother groups
    	    //	define set of delta_bi to store the increment of ICCs of the clusters\cup Cs

    	    int delta_bi;
    	    unsigned int bi_new;
    	    unsigned int bi_old;
    	    Cluster set_tem;
    	    //------
    	    //start to compare possible groups!

    		//output to test!
    		//	    std::cout<<"As to the overlapping node of "<<index_of_overlapping_node<<"\n";
    		//	    std::cout<<"the size of its mother group is "<<"\n";

    	    for(unsigned int i=0;i<Pool_MotherGroups[index_of_overlapping_node].size();i++)
    		{
    		unsigned int index_of_competitive_mothergroup = Pool_MotherGroups[index_of_overlapping_node][i];	//get index for possible mother group!

    		//查找 index_of_competitive_mothergroup 在 c 中的位置！这个位置也就是index_of_competitive_mothergroup 在Pool_Groups中的位置！
    		unsigned int counter=0;
    		unsigned int ii=0;
    		while(c[ii]!=index_of_competitive_mothergroup)
    		    {
    		    counter++;
    		    ii++;
    		    }
    		//backup把原来的cluster放入backuppool_original_group中第index_of_competitive_mothergroup个元素！
    		//此后，Pool_Groups不在是phaseI结束后的情况！！！
    		backuppool_original_group[index_of_competitive_mothergroup]=Pool_Groups[counter];
    		bi_old = CaculateNumCCC(Pool_Groups[counter]);

    		//output to test!
    		//		std::cout<<Pool_Groups[counter].size()<<"\t";

    		//delete overlapping nodes
    		std::vector<unsigned int>::iterator where = std::find(Pool_Groups[counter].begin(), Pool_Groups[counter].end(), index_of_overlapping_node);
    		Pool_Groups[counter].erase(where);	//new group, after deleting the overlapping node

    		    //	    //delete clusterhead(should have not been deleted, but in actual operation, clusterhead will be inner_product for 2 times), to save caculation time!
    		    //	    std::vector<unsigned int>::iterator where1 = std::find(Pool_Groups[counter].begin(), Pool_Groups[counter].end(), index_of_competitive_mothergroup);
    		    //	    Pool_Groups[counter].erase(where1);

    		//caculate some metrics of relavant clusters(after shrinking):
    		//two metrics:
    		//		1. delta_bi: num of common channels
    		//		2. clustersize;

    		bi_new = CaculateNumCCC(Pool_Groups[counter]);

    		delta_bi = bi_new - bi_old;
    		//+++++++++++	set_delta_bi obtained here	++++++++++
    		set_delta_bi.push_back(delta_bi);	//set_delta_bi中元素的个数和顺序同Pool_MotherGroups[index_of_overlapping_node]中的是一致的，只是具体的index不一样，可资利用！
    		//if there is smallest delta_bi exists, keep it inside the group, and recover group membership, then end!
    		//if not, compare num_hn, and find the biggest,

    		//cluster recover!
    		Pool_Groups[counter]=backuppool_original_group[index_of_competitive_mothergroup];
    		}



//    	    //output to test!
//    	    std::cout<<std::endl;
//    	    std::cout<<"the size of set_delta_bi of node "<<index_of_overlapping_node<<" is "<<set_delta_bi.size()<<","<<"\n";
//    	    std::cout<<"they are: "<<"\n";
//    	    for (unsigned int aa=0; aa < set_delta_bi.size(); aa++)
//    		{
//    		std::cout<<set_delta_bi[aa]<<"\t";
//    		}
//    	    std::cout<<std::endl;
    	    //there is one big problem, which we don't like:
    	    //debatable node is more likely leaves smaller clusters and stay in bigger clusters, because leaving smalelr cluster can more
    	    //likely increase the number of common channels!!!!!


    	    //here: begin to compare metrics
    	    //find the smallest delta_bi
    	    Cluster couter;
    	    //	Cluster couter2;
    	    unsigned int stone=set_delta_bi[0];
    	    for(unsigned int m=0;m<set_delta_bi.size();m++)
    		{
    		if(stone>set_delta_bi[m])
    		    {
    		    stone=set_delta_bi[m];
    		    }
    		}	//stone is the smallest element in set_delta_bi.

    	    for(unsigned int m=0; m < set_delta_bi.size(); m++)
    		{
    		if (stone == set_delta_bi[m])
    		    {
    		    couter.push_back(m);//record m (relative position index in set_delta_bi, corresponding to the 'mother group' in Pool_MotherGroups[index_of_overlapping_node].)
    		    }
    		}

    	     unsigned int index_realmother;
    	    //	 unsigned int m=couter[0];

    	     if(couter.size()>1)//means: there are more than one mother clusters with the identical small delta_bi.
    		{

    		 unsigned int smallestsize = Pool_Groups[Pool_MotherGroups[index_of_overlapping_node][couter[0]]].size();

    		 //find the one with smallest size out of the possible mother clusters which have identical delta_bi and hn
    		 for(unsigned int s=0; s<couter.size(); s++)
    		    {
    		     if(smallestsize > Pool_Groups[Pool_MotherGroups[index_of_overlapping_node][couter[s]]].size())
    			 {
    			 smallestsize = Pool_Groups[Pool_MotherGroups[index_of_overlapping_node][couter[s]]].size();
    			 }
    		    }

    		 for(unsigned int t=0; t<couter.size(); t++)
    		    {
    		     if(Pool_Groups[Pool_MotherGroups[index_of_overlapping_node][couter[t]]].size()==smallestsize)
    			 {
    			 index_realmother=Pool_MotherGroups[index_of_overlapping_node][couter[t]];
    			 tem_mum[index_of_overlapping_node]=index_realmother;
    //			 std::cout<<"Justice about size was done!"<<"\n";
//    			 std::cout<<"debatable node "<<index_of_overlapping_node<<" belongs to cluster "<<index_realmother<<"\n";
    			 //just ignore the repeated output that"debatable node xx belongs to cluster xxx", because anyway, each debatable
    			 //will join into only one cluster finally!

    			 //get rid of this node from the chosen cluster!

    //				unsigned int counter2=0;
    ////				unsigned int ii=0;
    //				while(c[counter2]!=index_of_competitive_mothergroup)
    //				    {
    //				    counter2++;
    ////				    ii++;
    //				    }
    //
    //			std::vector<unsigned int>::iterator where = std::find(Pool_Groups[counter2].begin(), Pool_Groups[counter2].end(), index_of_overlapping_node);
    //			Pool_Groups[counter2].erase(where);
    			 }
    		    }
    		}
    	    else
    		{
    		index_realmother=Pool_MotherGroups[index_of_overlapping_node][couter[0]];
    		tem_mum[index_of_overlapping_node]=index_realmother;
    		std::cout<<"debatable node "<<index_of_overlapping_node<<" belongs to cluster "<<index_realmother<<"\n";

    		}

    	    //找到delta_bi最小的当中clustersize最小的那个，就把overlapping node留在那里.因为目前所有的group都是原处的规模，所以要把其他mother groups里的overlapping ndoe删除！
    	    for(unsigned int j=0; j<Pool_MotherGroups[index_of_overlapping_node].size();j++)
    		{
    		unsigned int index_cc=Pool_MotherGroups[index_of_overlapping_node][j];
    		if(index_cc!=index_realmother)
    		    {
    		    unsigned int ii_3=0;
    		    while(c[ii_3]!=index_cc)
    			{
    			ii_3++;
    			}
    		    std::vector<unsigned int>::iterator where = std::find(Pool_Groups[ii_3].begin(), Pool_Groups[ii_3].end(), index_of_overlapping_node);
    		    Pool_Groups[ii_3].erase(where);
    		    }
    		}
    	    set_delta_bi.clear();
    	    }	//The initial membership clarification is finished.
    		//All the debatable nodes have clarified their membership after local fast algorithm.


    	unsigned int NewNumOfAllNodes = 0;
    	std::cout<<"there are "<<Pool_Groups.size()<<"group(s)"<<"\n";
    	//	std::cout<<"there are "<<c.size()<<" group(s)"<<"\n";
    		for(unsigned int i=0;i<Pool_Groups.size();i++)
    		    {
    		    std::cout<<"Group("<<c[i]<<")= "<<"\t";
    		    for(unsigned int j=0;j<Pool_Groups[i].size();j++)
    			{
    			std::cout<< Pool_Groups[i][j] <<"\t";
    			NewNumOfAllNodes++;
    			}
    		    std::cout<<"\n";
    		    }
    	std::cout<<"There are "<<NewNumOfAllNodes<<" nodes after initial clarification!"<<"\n";

    	//After this change, a cluster of node i is Pool_Groups_trueindex[i]
    	//Pool_Groups	--->	Pool_Groups_trueindex

    //	unsigned int index=0;

    	for(unsigned int u=0;u<Pool_Groups_trueindex.size();u++)
    	    {
    	    Pool_Groups_trueindex[u].clear();
    	    }

/*
 * Pool_Groups_trueindex is vector of clusters,
 * the index of element in Pool_Groups_trueindex is the cluster head,
 * which means in Pool_Groups_trueindex, there are empty element (also is vector).
 */

    	for(unsigned int i=0;i<Pool_Groups.size();i++)
    	    {
    	    for(unsigned int j=0; j<Pool_Groups[i].size(); j++)
    		{
    		for (unsigned int k=0; k<c.size(); k++)
    		    {
    		    if(c[k]==Pool_Groups[i][j])
    			{
    			Pool_Groups_trueindex[c[k]]=Pool_Groups[i];
    			}
    		    }
    		}
    	    }

    	////////////////////////////////////////////////////////////////////////////////
    	//update membership
    	////////////////////////////////////////////////////////////////////////////////
    	//Now, it is time to update their membership to improve sum_

    	unsigned int flag_update=1;
    	unsigned int index_realmother_update;


while(flag_update)
	{
	flag_update=0;
	UpdateRound++;
	for(unsigned int i=0; i<real_overlapping_nodes.size(); i++)
	    {
	    unsigned int debatable_node = real_overlapping_nodes[i];	//get authentic index for overlapping node!
	    //	begin to exame and compare every possible mother groups
	    //	define set of delta_bi
	    Cluster new_mems;	//set of index of clusters which can possiblily be debatable node's new clusters.
	    //delta_a is the increment of ICCs in C, i\in C.
	    //delta_a_0 is the K_C
	    //delta_a_1 is the K_{C\i}
//test
//	    unsigned int clusterhead=tem_mum[debatable_node];

	    unsigned int indexofmothercluster=tem_mum[debatable_node];
	    unsigned int delta_a_0=CaculateNumCCC(Pool_Groups_trueindex[tem_mum[debatable_node]]);
	    for(unsigned i=0; i<Pool_Groups_trueindex[tem_mum[debatable_node]].size();i++)
		{
		unsigned int x=Pool_Groups_trueindex[tem_mum[debatable_node]][i];
		}
	    //delete the debatable node from the initial cluster to see the increment of ICCs
	    std::vector<unsigned int>::iterator w1 = std::find(Pool_Groups_trueindex[tem_mum[debatable_node]].begin(), Pool_Groups_trueindex[tem_mum[debatable_node]].end(), debatable_node);
	    Pool_Groups_trueindex[tem_mum[debatable_node]].erase(w1);
	    unsigned int delta_a_1=CaculateNumCCC(Pool_Groups_trueindex[tem_mum[debatable_node]]);
	    unsigned int delta_a=delta_a_1-delta_a_0;
	    //recover Pool_Groups[tem_mum[debatable_node]]:
	    Pool_Groups_trueindex[tem_mum[debatable_node]].push_back(debatable_node);

	    //get rid of cluster C from the S, C\subset S, i\in C.
	    std::vector<unsigned int>::iterator w2 = std::find(Pool_MotherGroups[debatable_node].begin(), Pool_MotherGroups[debatable_node].end(), tem_mum[debatable_node]);
	    Pool_MotherGroups[debatable_node].erase(w2);


	    for(unsigned int j=0; j<Pool_MotherGroups[debatable_node].size();j++)
		{
		unsigned int delta_b_0=CaculateNumCCC(Pool_Groups_trueindex[Pool_MotherGroups[debatable_node][j]]);
		Pool_Groups_trueindex[Pool_MotherGroups[debatable_node][j]].push_back(debatable_node);
		unsigned int delta_b_1=CaculateNumCCC(Pool_Groups_trueindex[Pool_MotherGroups[debatable_node][j]]);
		unsigned int delta_b=delta_b_0-delta_b_1;
		set_delta_icc.push_back(delta_b);

		Pool_Groups_trueindex[Pool_MotherGroups[debatable_node][j]].pop_back();
		}
	    //recover Pool_MotherGroups[debatable_node]!
	    Pool_MotherGroups[debatable_node].push_back(tem_mum[debatable_node]);
	    //

	    if(*min_element(set_delta_icc.begin(),set_delta_icc.end())<delta_a)
		//this means update brings profits, the minial decrement of ICCs is smaller than the increment.
		//notice: in the following, 'set_same_delta_icc' stores the index of C\subset S_i, which have the same increment of ICCs in the process of update
		{
		//debatable node quits from the previous cluster:
		flag_update=1;
		std::vector<unsigned int>::iterator w3 = std::find(Pool_Groups_trueindex[tem_mum[debatable_node]].begin(), Pool_Groups_trueindex[tem_mum[debatable_node]].end(), debatable_node);
		Pool_Groups_trueindex[tem_mum[debatable_node]].erase(w3);

		for(unsigned int k=0; k<set_delta_icc.size(); k++)
		    {
		    if(*min_element(set_delta_icc.begin(),set_delta_icc.end()) == set_delta_icc[k])
			{
			set_same_delta_icc.push_back(Pool_MotherGroups[debatable_node][k]);
			}
		    }

		//find the cluster C\subset set_same_delta_icc according to cluster size.
		for(unsigned int l=0; l<set_same_delta_icc.size(); l++)
		    {
		    set_clustersize.push_back(Pool_Groups_trueindex[set_same_delta_icc[l]].size());
		    }
		for(unsigned int m=0; m<set_clustersize.size();m++)
		    {
		    if(*min_element(set_clustersize.begin(),set_clustersize.end())==set_clustersize[m])
			{
			index_realmother_update=set_same_delta_icc[m];
			}
		    }
		/////debatable node joins into cluster of index_realmother_update, and quit from cluster
		tem_mum[debatable_node]=index_realmother_update;
		Pool_Groups_trueindex[tem_mum[debatable_node]].push_back(debatable_node);
		}

	    //update can happen if debatable node can transfer from a larger cluster to a smaller one!
	    else if(*min_element(set_delta_icc.begin(),set_delta_icc.end())==delta_a)
		{
		for(unsigned int n=0; n<set_delta_icc.size(); n++)
		    {
		    if(*min_element(set_delta_icc.begin(),set_delta_icc.end()) == set_delta_icc[n])
			{
			set_same_delta_icc.push_back(Pool_MotherGroups[debatable_node][n]);
			}
		    }
		for(unsigned int l=0; l<set_same_delta_icc.size(); l++)
		    {
//		    unsigned int index7=set_same_delta_icc[l];
//		    unsigned int size1=Pool_Groups_trueindex[set_same_delta_icc[l]].size();
		    set_clustersize.push_back(Pool_Groups_trueindex[set_same_delta_icc[l]].size());
		    }
		unsigned int smallersize=*min_element(set_clustersize.begin(),set_clustersize.end());
		//	\exist C\subset S_i, with C.size()=C_{tem_mum}.size(), Delta C=DeltaC_{tem_mum}.
		if(smallersize<Pool_Groups_trueindex[tem_mum[debatable_node]].size()-1)
		    {
		    flag_update=1;

		    //debatable node quits from the previous cluster:
		    std::vector<unsigned int>::iterator w4 = std::find(Pool_Groups_trueindex[tem_mum[debatable_node]].begin(), Pool_Groups_trueindex[tem_mum[debatable_node]].end(), debatable_node);
		    Pool_Groups_trueindex[tem_mum[debatable_node]].erase(w4);

		    for(unsigned int m=0; m<set_clustersize.size(); m++)
			{
			if(smallersize==set_clustersize[m])
			    {
			    index_realmother_update=set_same_delta_icc[m];
			    }
			}
		    tem_mum[debatable_node]=index_realmother_update;
		    Pool_Groups_trueindex[tem_mum[debatable_node]].push_back(debatable_node);
		    }

		}//end elseif

		set_delta_icc.clear();
		set_clustersize.clear();
		set_same_delta_icc.clear();

	    }//finish all debatable nodes!

	/*
	 * convert Pool_Groups_trueindex --> Pool_Groups
	 */
	for(unsigned int i=0; i < Pool_Groups.size(); i++) {
	    Pool_Groups[i].clear();
	}

	/*
	 * i for Pool_Groups_trueindex whose size is a const
	 * j for Pool_Groups
	 * k for Pool_Groups[j]
	 */
	std::cout<<"Pool_Groups_trueindex.size() =" << Pool_Groups_trueindex.size() <<"\n";
	unsigned int j=0;
	for(unsigned int i=0; i < Pool_Groups_trueindex.size(); i++) {
	    if(!Pool_Groups_trueindex[i].empty()) {
		unsigned int k=0;
//		std::cout<<"Pool_Groups_trueindex[" << i << "].size() =" << Pool_Groups_trueindex[i].size() <<"\n";
		while (k != Pool_Groups_trueindex[i].size()) {
//		    if(i==14 && k == 0)
//			{
//			int checkvalue=Pool_Groups_trueindex[i][k];
//			}
//		    std::cout<<"Pool_Groups_trueindex[" << i << "][" << k << "] =" << Pool_Groups_trueindex[i][k] <<"\n";
		    Pool_Groups[j].push_back(Pool_Groups_trueindex[i][k]);
		    k++;

		    }
		    j++;
		}
	}
	std::cout<<"Pool_Groups.size() =" << Pool_Groups.size() <<"\n";
	}//update complete!
#endif
//    std::cout<<"\n";
    std::cout<<"formed clusters after phase II is: "<<"\n";
//    std::cout<<"there are "<<Pool_Groups.size()<<" groups"<<"\n";
    int NumOfAllNodes=0;
    for(unsigned int i=0;i<Pool_Groups.size();i++)
    	    {
    	    //std::cout<<"Group("<<c[i]<<")= "<<"\t";
//	    std::cout<<"Pool_Groups.size() is " << Pool_Groups.size() << "\n";
    	    for(unsigned int j=0; j<Pool_Groups[i].size(); j++)
    		{
//    		std::cout<<"Pool_Groups[" << i << "].size() is " << Pool_Groups[i].size() << "\n";
    		std::cout<< Pool_Groups[i][j] <<"\t";
    		NumOfAllNodes++;
    		}
    	    std::cout<<"\n";
    	    }
    std::cout<<"Now, after clarification, the number of all nodes involved is "<<NumOfAllNodes<<"\n";
#endif



#ifdef SURVIVAL_CLUSTERS
#ifdef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME
    /*
     * read the clustering result generated from MATLAB.
     * modify Pool_Groups!
     */
    survial_ratio_from_centralized(seed);
#endif

    /*
     * newPU is the number of added PUs every time
     */
    unsigned int oldPU =30;
    unsigned int newPU =10;// number of PRs added each time, the total number is decided as nuwPU*(a number)
    unsigned int ChDim =10;
    unsigned int RadiusPR = 20;

    std::ofstream myfile;
    myfile.open("survival.cls", std::ios::app);

    /*
     * the total number of PUs is controlled here. currently, add PUs for 10 times.
     */
    unsigned int num_unclustered=0;
    for(unsigned int i=0; i< Pool_Groups.size(); i++) {
//		myfile << "Pool_Groups[" << i << "].size()= " << Pool_Groups[i].size() << "\n";

	if(Pool_Groups[i].size()==1)
	    {
	    num_unclustered++;
	    }
    }
    myfile << num_unclustered << "\t";

    while(PRGroup.size() < oldPU + newPU*30 ){

	//add new PUs in the network, localization
	    for (unsigned int j = 0; j < newPU; j++)		//each node node within the possible scope
    		    {
		    Point p(2);
		    for (unsigned int i = 0; i < 2; i++)
			{
			//points distributed in the plane of 100m X 100m
			p[i]=network_geo_edge *(rand() * (1.0) / RAND_MAX);
			}
		    PRGroup.push_back(p);

	    /*
	     * initialize channel usage of the new PUs
	     */
	    for(unsigned int k=PRGroup.size(); k<PRGroup.size()+newPU; k++)
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
    		}


	    /*
	     * remaning clusters under initial PUs and new PUs
	     * check each cluster in Pool_Groups
	     */


	    num_unclustered=0;
//	    myfile << "Pool_Groups.size() = " << Pool_Groups.size() << "\n";

	    int checksize = Pool_Groups.size();
	    for(unsigned int i=0; i< Pool_Groups.size(); i++) {
//		myfile << "Pool_Groups[" << i << "].size()= " << Pool_Groups[i].size() << "\n";

		if(Pool_Groups[i].size()==1)
		    {
		    num_unclustered++;
		    }
		else
		    {
		    /*
		     * check whether this cluster remains
		     */
		    unsigned int survival = survival_check(Pool_Groups[i], RadiusPR, ChDim);
		    if(!survival)
			{
			num_unclustered += Pool_Groups[i].size();
			}
		    }
		}

	    /*
	     * print out the number of unclustered nodes to file
	     */
	    myfile << num_unclustered << "\t";
    }
    myfile << "\n";

#endif





#ifdef COMPARE_SPECTRUM_WITH_GROUND_TRUTH_METHOD2
    for(unsigned int i=0; i< Pool_Groups.size(); i++) {

	if(Pool_Groups[i].size()>1)
	    {
	    /*
	     * check whether the non-singleton cluster remains
	     */
	    unsigned ccc= CaculateNumCCC(Pool_Groups[i]);
	    unsigned ccc_copy= CaculateNumCCC_v_all_copy(Pool_Groups[i]);
	    if(!ccc_copy && ccc){
		for(unsigned int k =0; k <Pool_Groups[i].size(); k++){
		    Cluster singleton_cluster;
		    singleton_cluster.push_back(Pool_Groups[i][k]);
		    Pool_Groups.push_back(singleton_cluster);
		}
		    Pool_Groups.erase(Pool_Groups.begin()+i);
	    }

	    }
	}
#endif

    /*
     * Check the clusters in Pool_Groups, whose sizes are larger than 1.
     *
     * check whether the cluster sustains, if not, renew Pool_Groups, if yes, the CCC in the cluster.
     */
#ifdef _DOUBLE_CHECK_UNACCURATE_SPECTRUM_SENSING

    for(unsigned int i=0; i< Pool_Groups.size(); i++) {

	if(Pool_Groups[i].size()>1)
	    {
	    /*
	     * check whether the non-singleton cluster remains
	     */
	    float false_negative;
	    false_negative = FALSE_NEGATIVE_SPECTRUM_SENSING;
	    cluster_with_accurate_spectrum(Pool_Groups[i], i, false_negative, 10);
	    }
	}
#endif


#ifdef COMPARE_SPECTRUM_WITH_GROUND_TRUTH_METHOD2
v_all = v_all_copy;
#endif

//====
    //std::cout<<"\n";
    //std::cout<<"number of ccc and outwards cc are as follows: "<<"\n";
    double SumNumCCC=0; // sum of CCC over clusters whose sizes are bigger than 1
    Cluster set_SumNumCCC; // sum of CCC over clusters whose sizes are bigger than 1

    double SumNumCCC_allsize=0;
    Cluster set_SumNumCCC_allsize;

    double SumNumOutCC=0;
    Cluster set_SumNumOutCC;

    unsigned int SumOCC_onecluster=0;
    unsigned int SumOCC_onenode=0;

//output #CCC and #OCC!
    for(unsigned int i=0;i<Pool_Groups.size();i++)
	{
	//std::cout<<"inner ccc and outwards cc of Group("<<c[i]<<")= "<<"\t";
	unsigned int NumCCC=CaculateNumCCC(Pool_Groups[i]);
	if(Pool_Groups[i].size()>1)
		{
		SumNumCCC += NumCCC;
		set_SumNumCCC.push_back(NumCCC);
		}
	SumNumCCC_allsize += NumCCC;
	set_SumNumCCC_allsize.push_back(NumCCC);

	prod_size_ccc += NumCCC * Pool_Groups[i].size();

//	//output #CCC
//	//std::cout<<NumCCC<<"\t";
//	unsigned int inner=0;
//
//	for(unsigned int j=0; j<Pool_Groups[i].size(); j++)	//Pool_Groups[i][j] is a member of group i, we want to find the number of ccc between it and other nodes!
//								//actually, Pool_Groups[i][j] can possibly be clusterhead, but that is OK!
//	    {
//	    Cluster tem_neighbor= Pool_Neighbors[Pool_Groups[i][j]];	//node Pool_Groups[i][j]'s neighborhood
//	    Cluster tem_group = Pool_Groups[i];				//the group where node Pool_Groups[i][j] is
//	    Cluster tem_tem_neighbor=tem_neighbor; //tem_tem_neighbor is unchanged!
//
//
//	    std::vector<unsigned int>::iterator where = std::find(tem_group.begin(), tem_group.end(), Pool_Groups[i][j]);
//	    tem_group.erase(where);
//
//	    //elimilate Pool_Groups[i][j]'s neighbors locating in the same group with Pool_Groups[i][j].
//	    for(unsigned int k=0; k<tem_group.size(); k++)
//		{
//		for(unsigned m=0; m<tem_neighbor.size(); m++)
//		    {
//		    if (tem_neighbor[m]==tem_group[k])
//			{
//			std::vector<unsigned int>::iterator where = std::find(tem_neighbor.begin(), tem_neighbor.end(), tem_group[k]);
//			tem_neighbor.erase(where);
//			}
//		    }
//		}
//	    tem_tem_neighbor=tem_neighbor;
//	    //this restriction is cancled, because cluster node can also connect with other node, not only border nodes!
//	    //elimilate Pool_Groups[i][j]'s neighbors which are clusterhead of other clusters
//
////	    for(unsigned int k=0; k<c.size(); k++)
////		{
////		for(unsigned m=0; m<tem_tem_neighbor.size(); m++)
////		    {
////		    if (tem_tem_neighbor[m]==c[k])
////			{
////			if(Pool_Groups[c[k]].size()>1)
////			    {
////			    std::vector<unsigned int>::iterator where = std::find(tem_neighbor.begin(), tem_neighbor.end(), tem_tem_neighbor[m]);
////			    tem_neighbor.erase(where);
////			    }
////			}
////
////		    }
////		}
//
//
////	    for(unsigned int l=0; l< tem_neighbor.size(); l++)
////		{
////		inner=inner+inner_product(v_all_copy[tem_neighbor[l]].begin(), v_all_copy[tem_neighbor[l]].end(),v_all_copy[Pool_Groups[i][j]].begin(), 0);
////		}
//	    //make sure the intermediate vector are new!
//	    Channel_onenode v_tem(10,0);
//	    Channel_onenode v_onebordernode(10,0);
//	    for(unsigned int l=0; l< tem_neighbor.size(); l++)
//		{
//		for(unsigned int m=0; m< 10; m++)
//		    {
//		    v_tem[m] = v_all_copy[Pool_Groups[i][j]][m]*v_all_copy[tem_neighbor[l]][m];
//		    }
//		for(unsigned int n=0; n< 10; n++)
//		    {
//		    v_onebordernode[n]=v_onebordernode[n]||v_tem[n];
//		    }
//		}
//	    //all the usable channels in "one node" on which conncections can be bulit!
//	    unsigned int SumOCC_onenode = std::accumulate( v_onebordernode.begin(), v_onebordernode.end(), 0 );
//	    //the sum of usable channels in cluster Pool_Groups[i]!
//	    SumOCC_onecluster = SumOCC_onecluster + SumOCC_onenode;
//	    }
//	set_SumNumOutCC.push_back(SumOCC_onecluster);
//
//	//output #OCC
//	//std::cout<<SumOCC_onecluster<<"\t";
//	//std::cout<<"\n";
//	SumOCC_onecluster=0;
//	SumOCC_onenode=0;
//	inner=0;
	}

//    protected member variats of class
//	float average_CCC;
//	float average_OCC;
//	float number_clusters;
//	float average_size;
//	float cv_size;
//	float cv_CCC;
//	float cv_OCC;
//    //
	for(unsigned int i=0; i<100;i++)
	    {
	    size_distr.push_back(0);
	    CCC_distr.push_back(0);
	    OCC_distr.push_back(0);
	    }


    //---here to record the distribution of #CCC
    //set_SumNumCCC
    for(unsigned int i=0; i < set_SumNumCCC.size(); i++)
	{
	CCC_distr[set_SumNumCCC[i]]=CCC_distr[set_SumNumCCC[i]]+1;
	}


//    /*
//     * Average CCC with consideration of non-single node cluster
//     */


    /*
     * Average CCC considering only clusters whose sizes are bigger than 1
     */
    average_CCC = SumNumCCC/set_SumNumCCC.size();
    average_CCC_allsize = SumNumCCC_allsize/set_SumNumCCC_allsize.size();

    //std::cout<<"Average number of inner CCC is "<<average_CCC<<"\n";

    unsigned int SumOCC = std::accumulate( set_SumNumOutCC.begin(), set_SumNumOutCC.end(), 0 );
    //---here to record the distribution of #OCC
    for(unsigned int i=0; i < set_SumNumOutCC.size(); i++)
	{
	OCC_distr[set_SumNumOutCC[i]]=OCC_distr[set_SumNumOutCC[i]]+1;
	}

    //the average usable channels in each cluster!
    average_OCC=SumOCC/(Pool_Groups.size()*1.0);
    //std::cout<<"Average number of outwards cc per node is "<< average_OCC <<"\n";

    number_clusters=Pool_Groups.size();
    //std::cout<<"number of clusters in this topology is: "<< number_clusters <<std::endl;
    //inline double Clusters::CaculateCV(Cluster group)

//    Cluster set_size;
    std::vector<unsigned int> set_size;
    for(unsigned int i=0; i<Pool_Groups.size(); i++)
		{
//    	std::cout<<"-------Pool_Groups[" << i << "].size()---------" << Pool_Groups[i].size() << "-------------------\n";
//    	std::cout<<"------set_size.size()----------" << set_size.size() << "-------------------\n";
//    	std::cout << typeid(set_size).name() << std::endl;
    	unsigned int x = Pool_Groups[i].size();
		set_size.push_back(x);
		}
    std::cout<<"-----------------1---------------" << "\n";
    //---here to record the distribution of cluster size
    for(unsigned int j=0; j<set_size.size(); j++)
		{
		size_distr[set_size[j]]=size_distr[set_size[j]]+1;
		}
    std::cout<<"--------------2-------------------" << "\n";
    average_size = std::accumulate(set_size.begin(),set_size.end(),0)/(set_size.size()*1.0);
    //std::cout<<"average size of clusters is: "<< average_size <<std::endl;

	average_size_allsize = SumNumCCC_allsize/set_SumNumCCC_allsize.size();
	    std::cout<<"--------------3-------------------" << "\n";


    //========output CV of cluster size
    cv_size=CaculateCV(set_size);
    //std::cout<<"CV of cluster size is: "<<cv_size<<std::endl;
    std::cout<<"--------------3-------------------" << "\n";

    //========output CV of CCC and OCC
    cv_CCC=CaculateCV(set_SumNumCCC);
    //std::cout<<"CV of CCC is: "<< cv_CCC <<std::endl;
//    cv_OCC=CaculateCV(set_SumNumOutCC);
//    //std::cout<<"CV of OCC is: "<< cv_OCC <<std::endl;
    std::cout<<"--------------3-------------------" << "\n";
    /*
     * 14.5
     * ======output average and CI of number of neighbors
     */
    Cluster group_neighborSize;
    for(unsigned int i=0; i < Pool_Neighbors.size(); i++)
    {
    	group_neighborSize.push_back(Pool_Neighbors[i].size());
    }
    average_neighborSize = std::accumulate( group_neighborSize.begin(), group_neighborSize.end(), 0 )/(group_neighborSize.size()*1.0);
    //std::cout<<"Average number of neighbors is "<< average_neighborSize <<"\n";
    cv_neighborSize = CaculateCV(group_neighborSize);
    //std::cout<<"CV of Average number of neighbors is: "<< cv_neighborSize <<std::endl;

    //std::cout<<"ROSS finished!"<<std::endl;
}//phaseII end

#ifdef _DOUBLE_CHECK_UNACCURATE_SPECTRUM_SENSING
/*
 * When the false negative ratio is 10%, it is 10% that the channel which is deemed as available (primary user is absent)
 * is actually not available.
 *
 * This function changes the Pool_Groups and v_all.
 */
void Clusters::cluster_with_accurate_spectrum(Cluster cluster_checked, unsigned int index, float false_negative, unsigned int ChDim){
    unsigned int flag =0;
    for(unsigned int m=0; m<ChDim; m++)
	{
	unsigned int availability_one_channel =1;
	for(unsigned int i =0; i<cluster_checked.size(); i++)
	    {
	    if(v_all[cluster_checked[i]][m]){
		v_all[cluster_checked[i]][m] = v_all[cluster_checked[i]][m] * (((double) rand() / (RAND_MAX))>false_negative ? 1:0);
	    }
	    availability_one_channel *= v_all[cluster_checked[i]][m];
	    }
	flag+=availability_one_channel;
	}

    if(!flag){
	// dissolve the cluster
	for(unsigned int i =0; i <cluster_checked.size(); i++){
	    Cluster singleton_cluster;
	    singleton_cluster.push_back(cluster_checked[i]);
	    Pool_Groups.push_back(singleton_cluster);
	}
	    Pool_Groups.erase(Pool_Groups.begin()+index);
    }
}
#endif

void Clusters::survial_ratio_from_centralized(unsigned int seed) {
	    /*
	     * clear Pool_Groups, note the memory is not released.
	     */
	for(unsigned int i=0; i< Pool_Groups.size(); i++) {
	    Pool_Groups[i].clear();
	}


	    /*
	     * read one file which correspondes to one topology
	     */

	    std::string fileName;
	    std::string seed_index = SSTR(seed);
	    fileName = "potential_clusters_20_" + seed_index + ".txt_3_0.5cluster.txt";

	    const char *str = fileName.c_str();
	    std::ifstream is(str);

	    int i =0;
	    std::string line;
	    while (std::getline(is, line)) {
		std::istringstream iss(line);
	    	int j=0;
	    	unsigned int x;
	    	Cluster tempVector;
	    	while (iss >> x) {
	    		tempVector.push_back(x-1); // decrease 1, because the id is from MATLAB where index starts from 1
	    		j++;
	    		}
	    	Pool_Groups.push_back(tempVector);
	    	i++;
	    }

	    /*-----------------
	     * add the singleton clusters in 20 node network
	     */

	    Cluster network;
	    unsigned int index =0;

	    /*
	     * construct a 20 node network
	     */
	    while (index < 20) {
		network.push_back(index);
		index++;
	    }

	    /*
	     * find out the singleton clusters
	     * traversal all 20 nodes, and delete the nodes which are already in non-singleton clustesrs
	     */
	    for (unsigned inti =0; i < network.size(); i++) {
		for (unsigned int j=0; j < Pool_Groups.size(); j++) {
		    for (unsigned int k=0; k < Pool_Groups[j].size(); k++) {
		    std::vector<unsigned int>::iterator where = std::find(network.begin(), network.end(), Pool_Groups[j][k]);
		    network.erase(where);
		    }
		}
	    }
	    /*
	     * add the singleton clusters into Pool_Groups.
	     */
	    for (unsigned inti =0; i < network.size(); i++) {
		Cluster singleton;
		singleton.push_back(network[i]);
		Pool_Groups.push_back(singleton);
	    }

	    std::cout << Pool_Groups.size() <<"\n";
	    is.close();                // close file

}




/*
 * v_all contains initial channel availability before adding new PUs
 * 21.07.2014
 */

unsigned int Clusters::survival_check(Cluster cluster_checked, unsigned int RadiusPR, unsigned int ChDim)
    {
    for(unsigned int i =0; i<cluster_checked.size(); i++)
	{
	unsigned int suID = cluster_checked[i];

	for(unsigned int k=0; k<PRGroup.size(); k++)
	    {
	    double dis= sqrt( (_ps[suID][0]-PRGroup[k][0]) * (_ps[suID][0]-PRGroup[k][0]) + (_ps[suID][1]-PRGroup[k][1]) * (_ps[suID][1]-PRGroup[k][1]) );
	    if(dis <= RadiusPR)
		    {
		    for(unsigned int m=0; m<ChDim; m++)
			    {
			    v_all[suID][m] = v_all[suID][m] * v_all_PR[k][m];
			    }
		    }
	    }
	}


    unsigned int flag =0;
    for(unsigned int m=0; m<ChDim; m++)
	{
	unsigned int availability_one_channel =1;
	for(unsigned int i =0; i<cluster_checked.size(); i++)
	    {
	    availability_one_channel *= v_all[cluster_checked[i]][m];
	    }
	flag+=availability_one_channel;
	}


    return flag;
    }


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
         //std::cout<<"following nodes have included by more than one cluster:"<<"\n"<<std::endl;
         for(unsigned int i=0;i<Pool_MotherGroups.size();i++)
    	 {
    	 //std::cout<<"node "<<overlapping_nodes[i]<<"----> (";
    	 for(unsigned int j=0;j<Pool_MotherGroups[i].size();j++)
    	 {//std::cout<<pool_mothergroups[i][j]<<"\t";}
    	 //std::cout<<"\n";
    	 }
*/

////输出操作符 << 重载
//	// single point output
//	std::ostream& operator<<(std::ostream& o, const Point& p)
//	{
//		o << "{ ";
//		BOOST_FOREACH(Point::value_type x, p)
//		{
//			o <<" "<< x;
//		}
//	        o << " }";
//	        o<< std::endl;
//		//o << " }<<"\n";
//		return o;
//	}
//
//	// single cluster output
//	std::ostream& operator<<(std::ostream& o, const Cluster& c)
//	{
//		o << "[ ";
//		BOOST_FOREACH(PointId pid, c)
//		{
//			o << " " << pid;
//		}
//		o << " ]";
//
//		return o;
//	}
//
//	// clusters output
//	std::ostream& operator<<(std::ostream& o, const Clusters& cs)
//	{
//		ClusterId cid = 1;
////----------di
////		//std::cout<<"the size of clusters is: "<<sizeof cs<<std::endl;
////----------
//		// interate for ||cs._clusters|| times.
//		BOOST_FOREACH(Cluster c, cs._clusters)
//		{
//			o << "c(" << cid++ << ")=";
//
//			BOOST_FOREACH(PointId pid, c)
//			{
//				o << cs._ps[pid];
//			}
//			o << std::endl;
//		}
//		return o;
//	}

};
