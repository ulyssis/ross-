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

unsigned int numCR = 500;
	Channel_onenode CaculateInit;	//a mediate vector for calculating bi

	PoolofNeighborhoods Pool_Neighbors;		//所走node的group的集合
	PoolofNeighborhoods Pool_Groups;		//basin nodes的group的集合

	//所有nodes的mother group的集合, 包括那些只有一个mother group的node
	//its element is node(index)'s mother ndoes
	PoolofNeighborhoods Pool_MotherGroups(1000,Neighbors());
	PoolofNeighborhoods Pool_MotherGroups_copy(1000,Neighbors());

	PoolofNeighborhoods backuppool_original_group(1000,Neighbors());//store original group in the phaseII

	robustness_vector ai_vector;
	robustness_vector bi_vector;

	Cluster c;   // a new cluster for storing nodes serving as clusterheads
	Cluster NodesAfterSearch(1000,0);
//	Cluster NodesSurrender(100,0);
	Cluster set_surrender;	//record the nodes which channged their manors, phaseII
//	Cluster NodesInPhaseII(100,0);

	Cluster overlapping_nodes;
	Cluster real_overlapping_nodes;

//	std::vector<int> Nodefinalized(100,0);
//	std::vector<int> Nodefinalized_AsHead(100,0);

	Channel_allnode v_all;
	Channel_allnode v_all_copy;
	Channel_allnode v_all_PR;
	Channel_allnode v_all_PR_copy;

	//phaseI 的遍历遍数
	int Num=0;

	Points PRGroup;

	//---new for SOC

	PoolofClusters Pool_manor(numCR,Cluster());	//store every node's manor（领地）in phaseI and phaseII
	PoolofClusters Pool_manor_copy(numCR,Cluster());

	Cluster set_of_y(numCR, 0);	//contain the value of number of common channels
	Cluster set_of_size(numCR, 0);	//contain size of each node's manor
	Cluster set_Q(numCR, 0);

	//globlal variation to pass argument
	unsigned int NumberOfCR;
	double RadiusOfCR;

    Cluster set_size;
//    set_size.reserve(100);


//decide the location of nodes, the number is num_points.
void randomInit	(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR, unsigned int network_geo_size)
	{
//	srand(time(NULL));	//as long as not been excuted all together, this part of program will have new initialization!
//	srand(6);
    		for (unsigned int j = 0; j < NumCR; j++)		//each node node within the possible scope
		{
		  Point p(GeoDim);
//		  std::cout<<"point "<<j<<"=(";


//		  if (NumCR < 100)
//			  {
//			  for (unsigned int i = 0; i < GeoDim; i++)
//				{
//					//points distributed in the plane of 100m X 100m
//				p[i]=network_geo_size *(rand() * (1.0) / RAND_MAX);
//	//				std::cout <<p(i) << ' ';
//				}
//			  }
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
		  else
		  {
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


		  //----------------------------

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
		    p(i)=network_geo_size*(rand() * (1.0) / RAND_MAX);
	    //				std::cout <<p(i) << ' ';
		    }
	    //		  std::cout<<")"<<"\t";
		PRGroup.push_back(p);
	    //		std::cout << std::endl;
		}
		NumberOfCR = NumCR;
	}


// Neighbors is a vector of int
// using this structure, I can formulate into a clustering scheme based on distance!

inline Neighbors Clusters::findNeighbors(PointId pid, double radius)
{
	Neighbors ne;

	std::cout<<"Node "<<pid<<"'s neighbors are (";
	for (unsigned int j=0; j < _ps.size(); j++)
	{
		if 	((pid != j ) && (_dis(pid, j) < radius))
		{
			ne.push_back(j);
//			std::cout << "dis(" << pid  << "," << j << ")" << _dis(pid, j) << "<" << threshold << std::endl;
			std::cout<<j<<"\t";

		}
	}
	//output ne:
	std::cout <<")"<<std::endl;

	Pool_Neighbors.push_back(ne);

//	unsigned int w=Pool_Neighbors.size();
//	std::cout<<"there are "<<w<<" neighbors...of course, this number is the same with _ps"<<"\n";
	return ne;
};

//------------------------------------------------------------------------------------------------
//	1. Initialize channels randomly on PR nodes and obtain the available channels on CR nodes.
//	2. separate this function into channelRandnomInit and caculate_ai_bi, if have time
//	3. parameters: - cardinality of PR and CR (ChDim)    1 X parameter
// 		       - radius of PR and CR (RadiusPR, RadiusCR)	2 X parameter
//------------------------------------------------------------------------------------------------

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
		RadiusOfCR = RadiusCR;
	}

//------------------------------------
//这个function是为了计算一个group中common channel的数量
//variable is a group of nodes, the last element in the group is head node!
inline unsigned int Clusters::CaculateNumCCC(Cluster group)
    {
	  unsigned int NumCCC;
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

	  return NumCCC;
    }




//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//---	phaseI: calculate Pool_manor and set_of_y after initialing available channels on nodes!
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------

#ifdef CLUSTER_SIZE_CONTROL
	void Clusters::ClusteringPhaseI(unsigned int optimalSize)
#endif
#ifndef CLUSTER_SIZE_CONTROL
	void Clusters::ClusteringPhaseI()
#endif
	{
    for(unsigned int k=0; k<NumOfNodes; k++)
		{
		unsigned int x_k=0, y_k=0, Qvalue;

		Neighbors ne1;
		ne1 = findNeighbors(k, RadiusOfCR);

		//do not output to save time
		//output channel infomation!
		for(unsigned int n_ch=0; n_ch<v_all[k].size(); n_ch++)
			{
			std::cout<<v_all[k][n_ch]<<"\t";
			}
		std::cout<<"\n";

		Cluster Qvalue_set;
		Cluster Qset;
		PoolofClusters Qset_series;

		unsigned int index=0;

		Cluster innerset; //store the number of common channels between i and its neighbors
		CaculateInit=v_all[k];
		//put the node itself into Qset, which is the initial and smallest group
		//Qset is the available manor waiting for being selected!
		Qset.push_back(k);
		Qset_series.push_back(Qset); // the first cluster in Qset_series.
		unsigned int tem_inner;

		unsigned int Qvalue_max;

		/* get sets on the base of k's neighborhood ne1
		 * Note: we DON'T obtain all the possible groups of k, in contrary, we sequentially put neighbors of k into Qset in ascending sequence (on the number of common channels), and put all Qset into Qset_series.
		 */
		while(ne1.size() > 0)	//the size of ne1 keeps decreasing
			{
			// clear innerset, to store new group(the original group after elimination)
			innerset.clear();

			/*
			 * find innerset which stores the pairwise number of common channels between k and its neighbors.
			 */
			for (unsigned int i = 0; i < ne1.size(); i++)
				{
				unsigned int inner = inner_product( CaculateInit.begin(),CaculateInit.end(),v_all[ne1[i]].begin(), 0);
				innerset.push_back(inner);
				}
			//get the biggest inner_product between
			unsigned int max_inner = *max_element(innerset.begin(),innerset.end());

			// add the neighbor which has max_inner common channels with i into Qset, and eliminate it later on
			unsigned int break_beacon = 0;
			unsigned j = 0;
			while(!break_beacon)
				{
				tem_inner = inner_product( CaculateInit.begin(),CaculateInit.end(),v_all[ne1[j]].begin(), 0);
				if(tem_inner==max_inner)
					{
					Qset.push_back(ne1[j]);
					//for each node, get series of possible manor, then in the following step, we should decide which manor is the best!
					Qset_series.push_back(Qset);
					//the node sharing the most common channels with k is chosen into Qset, then get deleted
					//from ne1
					std::vector<unsigned int>::iterator where = std::find(ne1.begin(), ne1.end(), ne1[j]);
					ne1.erase(where);
					break_beacon=1;
					}
				j++;
				}
			}

		/*
		 * check Qset_series
		 * if Qset_series contains Qset whose size is greater than or equal to optimalSize, then delete the Qset whose size is smaller than optimalSize.
		 * elseif there is no Qset in Qset_series, whose size is bigger than optimalSize
		 */

#ifdef CLUSTER_SIZE_CONTROL

		unsigned int counter_delete = 0;

		Cluster QsetSize_series;
		for(unsigned int q=0; q < Qset_series.size(); q++)
			{
			QsetSize_series.push_back(Qset_series[q].size());
			}
		unsigned int maxSetSize = *max_element(QsetSize_series.begin(), QsetSize_series.end());

		if (maxSetSize >=optimalSize)
			{
			for(unsigned int q=0; q < Qset_series.size(); q++)
				{
				std::cout << "size of Qset_series[" << q << "] is " << Qset_series[q].size() << "\n";
				// the first condition is to delete big clusters, the second is to delete single node clusters
				if(Qset_series[q].size() > optimalSize)// || Qset_series[q].size() < optimalSize-1)
					{
					counter_delete++;
								std::cout << "size of Qset_series[" << q << "] is " << Qset_series[q].size() << "\n";

//					PoolofClusters::iterator where = std::find(Qset_series.begin(), Qset_series.end(), Qset_series[q]);
//					Qset_series.erase(where);
								Qset_series.erase(Qset_series.begin()+q);
					}
				}
			}
		std::cout << "-----number of times to delete amnors with unproper sizes: " << counter_delete << "\n";
#endif


	/*
	 * find the manor which has the biggest Qvalue
	 */
		for(unsigned int m=0; m < Qset_series.size(); m++)
			{
//			std::cout << "Qset_series[" << m << "] is " << Qset_series[m].size() << "\n";
			//caculate Q, get the biggest one and record the index of Qset! and in the same time, x1 and y1 are obtained!
			x_k = Qset_series[m].size();	//cluster size
			y_k = CaculateNumCCC(Qset_series[m]);	//number of common channels

			Qvalue = x_k * y_k;
			Qvalue_set.push_back(Qvalue);
			}

		Qvalue_max = *max_element(Qvalue_set.begin(), Qvalue_set.end());

		//for(unsigned int m=Qset_series.size()-1; m != 0 ; m--)
		for(unsigned int m=0; m < Qset_series.size(); m++)
			{
			x_k=Qset_series[m].size();	//number of members
			y_k=CaculateNumCCC(Qset_series[m]);	//number of common channels
			Qvalue = x_k * y_k;
			if(Qvalue == Qvalue_max)
				{
				index=m;			//m is the position index in Qset_series, which means the best group for node k!
				//*****
				//Pool_manor is the most important variant!
				//*****
				Pool_manor[k] = Qset_series[index];
				set_of_y[k] = CaculateNumCCC(Pool_manor[k]);
				set_of_size[k]=Pool_manor[k].size();
				set_Q[k] = Pool_manor[k].size() * set_of_y[k];
				}
			}
		}
    Pool_manor_copy = Pool_manor;	 //another pool, for use of extract final clusters

#ifdef PRINT_ON_SCREAN
    //do not output to save time
    std::cout <<"==========ClusteringPhase I result:=========="<<std::endl;
    for(unsigned int t=0; t < Pool_manor.size(); t++)
	{

	std::cout <<"the manor of node ("<<t<<") is: (";
	for(unsigned int v=0; v < Pool_manor[t].size(); v++ )
	    {
	    std::cout << Pool_manor[t][v]<<"\t";
	    }

	std::cout <<")"<<std::endl<<"Q is "<<set_Q[t]<<"\t"<<"\t"<<"size is "<< Pool_manor[t].size()<<"\t"<<"NumCC is "<<set_of_y[t]<<std::endl;
	}
#endif
    std::cout <<"==========ClusteringPhase I result end!=========="<<std::endl;
    std::cout <<std::endl;

    }	//ClusteringPhaseI end!




//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//---	SOC phase II	---
//---	Updating cluster membership information	---
//---
//---	when a node receives info from neighbors, this node will check whether itself
//---	belongs to a better clique(the clique of its possible mother nodes)
//
//---	It is easy to say this principle, but there is a fact that, a node can possibly
//---	receive info more than one neighbor.
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------

//For a node, traversal every manor to see whether this manor contains it!!!!!
inline void Clusters::MotherGroups_Update()
    {
    unsigned int Flag_SurrenderedNode = 0;
    //PoolofNeighborhoods Pool_MotherGroups(100,Neighbors());
    //clear Pool_MotherGroups(100,Neighbors())
    for(unsigned int m=0; m<Pool_MotherGroups.size(); m++)
	{
	Pool_MotherGroups[m].clear();
	}
    //finding Pool_MotherGroups
    for(unsigned int h=0;h < Pool_manor.size();h++)
	{
	for(unsigned int i=0;i < Pool_manor[h].size();i++)	//对于每个group里面的node：Pool_manor[h][i]
	    {
	    unsigned int index=Pool_manor[h][i];	//index is "index of node."

	    if(!NodesAfterSearch[index])
		{
		//---	authetification begin
		for(unsigned int j=0; j < Pool_manor.size(); j++)	//遍历每个manor里面的node：Pool_manor[j][k]
		    {
		    //it is not accountable if included by a surrendered node!
		    for(unsigned int l=0; l<set_surrender.size(); l++)
			{
			//check whether the manor head is a surrendered node
			if(set_surrender[l]==j)
			    {
			    Flag_SurrenderedNode = 1;
			    }
			}
			//it is not accountable if included by a surrendered node!
			if(Flag_SurrenderedNode==0)
			    {
			     for(unsigned int k=0; k < Pool_manor[j].size(); k++)
				{
				if(Pool_manor[j][k]==index)			// to see, who has node "index"
				    {
				    //collect the nodes whose manor contains the node index!
				    Pool_MotherGroups[index].push_back(j);
				    //该node在查找完之后被labeled！
				    NodesAfterSearch[index]=1;
				    }
				}
			    }
			Flag_SurrenderedNode=0;
		    }
		}
	    }
	 }	//	1st part of phase II ends!
		//	for each node "index", the set of nodes(including node "index") whose manor conclude this node has been obtained!
	//finding Pool_MotherGroups
    }



void Clusters::ClusteringPhaseII()
    {
    //from 0 to Pool_MotherGroups.size():
    unsigned int beacon_Q=0;

    Cluster set_biggestQ;
    Cluster set_biggestsize;
    Cluster set_biggerQ_samesize;

    unsigned int biggestsize_biggestQ = 0;
    unsigned int index_biggestsize_biggestQ=0;

    for(unsigned int i=0; i < Pool_MotherGroups.size(); i++)
	{
	//the update of MotherGroup is necesseary!
	MotherGroups_Update();

//	//do not output to save time
	//output before node i is absorbed into other clusters
#ifdef PRINT_ON_SCREAN
	std::cout <<"node ("<<i<<")'s real possible mother nodes are : (";
    	for(unsigned int v=0; v < Pool_MotherGroups[i].size(); v++ )
    	    {
    	    std::cout << Pool_MotherGroups[i][v]<<"\t";
    	    }
    	std::cout <<")"<<std::endl;
#endif
	//	find out the set of nodes having biggest Q!
	unsigned int biggestQ = 0;
	//test
	unsigned int test_size=Pool_MotherGroups[i].size();


	if(Pool_MotherGroups[i].size()>0)
	    {
	    //compare which element in Pool_MotherGroups[i](node i's possible mum) has the biggest clique!
	    for(unsigned int j = 0; j < Pool_MotherGroups[i].size(); j++)
		{
		//test
		//unsigned int tem = Pool_MotherGroups[i][j];

		//search for biggest Q:
		if(biggestQ < set_Q[Pool_MotherGroups[i][j]])
		    {
		    biggestQ = set_Q[Pool_MotherGroups[i][j]];
		    }
		}


	    for(unsigned int j=0; j < Pool_MotherGroups[i].size(); j++)
		{
		if(set_Q[Pool_MotherGroups[i][j]] == biggestQ)	//at least equal to Pool_MotherGroups[i][j]'s Q
		    {
		    beacon_Q++;
		    set_biggestQ.push_back(Pool_MotherGroups[i][j]);
		    }
		}

	    if(set_Q[i]<=biggestQ)
		{
		//there are at least one(or one) mother node having the same biggest Q, including the node v itself!!!!
		if(beacon_Q > 1)
		    {
		    //try one by one!
		    //check Q first ---> if(have one biggest Q) ---> index_biggest_size_biggestQ update
		    //		|
		    //		---> if(have servel biggest Q) ---> check size ---> if ("<")   index_biggest_size_biggerQ update
		    for(unsigned int k=0; k < set_biggestQ.size(); k++)
			{
			if(biggestsize_biggestQ < set_of_size[set_biggestQ[k]])
			    {
			    biggestsize_biggestQ = set_of_size[set_biggestQ[k]];
			    index_biggestsize_biggestQ = set_biggestQ[k];
			    }
			}

		    for(unsigned int k=0; k < set_biggestQ.size(); k++)	//find the node with biggest index out of those having biggest size
			{
			if(set_of_size[set_biggestQ[k]] == biggestsize_biggestQ)
			    {
			    if(index_biggestsize_biggestQ <= set_biggestQ[k])
				{
				index_biggestsize_biggestQ = set_biggestQ[k];
				}
			    }
			}

		    //找出在存在最大Q的mum nodes的的情况下，size最大的mum node(在有几个并列的情况下就记录最大的index)
		    Pool_manor[i] = Pool_manor[index_biggestsize_biggestQ];

		    //从其他试图拥有它的manor(set_biggestQ)中剔除！
		    for(unsigned int j = 0; j < set_biggestQ.size(); j++)
			{
			if(set_biggestQ[j]!=index_biggestsize_biggestQ)
			    {
			    //剔除！
			    std::vector<unsigned int>::iterator where = std::find(Pool_manor[set_biggestQ[j]].begin(), Pool_manor[set_biggestQ[j]].end(), i);
			    Pool_manor[set_biggestQ[j]].erase(where);
			    }
			}

		    if(i!=index_biggestsize_biggestQ)
			{
			set_surrender.push_back(i);
			}
		    set_of_y[i] = set_of_y[index_biggestsize_biggestQ];
		    set_Q[i] = set_Q[index_biggestsize_biggestQ];
		    set_of_size[i] = set_of_size[index_biggestsize_biggestQ];
		    }


		//there is only one mother node having the biggest Q!
		if(beacon_Q == 1)
		    {
		    index_biggestsize_biggestQ=set_biggestQ[0];
		    if(i!= index_biggestsize_biggestQ)
			{
			set_surrender.push_back(i);
			}
		    //test
		    //unsigned int xxx=set_biggestQ[0];

		    Pool_manor[i] = Pool_manor[index_biggestsize_biggestQ];

		    set_of_y[i] = set_of_y[set_biggestQ[0]];
		    set_Q[i] = set_Q[set_biggestQ[0]];
		    set_of_size[i] = set_of_size[index_biggestsize_biggestQ];
		    }
		}
	    }//decide node i's manor in phaseII


	set_biggestQ.clear();
	beacon_Q = 0;
	index_biggestsize_biggestQ=0;
	biggestsize_biggestQ=0;

	biggestsize_biggestQ = 0;
	index_biggestsize_biggestQ=0;

	for(unsigned a=0; a<NodesAfterSearch.size(); a++)
	    {
	    NodesAfterSearch[a]=0;
	    }
	} //for all nodes, whether to surrender is clear!

    //do not output to save time
    std::cout <<"==========ClusteringPhase II results:=========="<<std::endl;

    //do not output to save time
    //---output set_surrender
    std::cout <<"There are "<<set_surrender.size()<<" nodes channged their manor :"<<std::endl;;
    for(unsigned int u=0; u < set_surrender.size(); u++)
	{
	std::cout << set_surrender[u] << "\t";
	}
    std::cout << "\n";

#ifdef PRINT_ON_SCREAN
    //do not output to save time
    //---output Pool_manor
    for(unsigned int t=0; t < Pool_manor.size(); t++)
	{
	std::cout <<"the manor of node ("<<t<<") is: (";
	for(unsigned int v=0; v < Pool_manor[t].size(); v++ )
	    {
	    std::cout << Pool_manor[t][v]<<"\t";
	    }

	std::cout << ")"<<"\n";
	}
#endif
    std::cout <<"==========ClusteringPhase II results end!=========="<<std::endl;

    std::cout <<std::endl;

    }	//phase II ends!
	//每个group的第一个element是这个group的head！


//	=========	Caculate CV of a vector of numbers	=========
    inline double Clusters::CaculateCV(Cluster group)
	{
    	if (!group.size()){
    		return 0;
    	}
    	else{
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
	}


    //	=========	Caculate CI of a vector of numbers	=========

    	double Clusters::CaculateCV(doubleVector group){
        	if (!group.size()){
        		return 0;
        	}
        	else{
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
    	}


//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//---	SOC phase III
//---	Finalizing cluster membership	---
//
//---	check every member of each node(i)'s manor, if the member's manor fail to include the node i,
//---	node (i) will eliminate this member
//--------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------

//	Traversal every node's manor, to check whether the member's manor include the node!

void Clusters::ClusteringPhaseIII()
    {
    unsigned int beacon_belong;
    for(unsigned int i=0; i < Pool_manor.size(); i++)
		{
		//store Pool_manor[i] before operating on it!
		Cluster tem_manor=Pool_manor[i];
		for(unsigned int j=0; j < Pool_manor[i].size(); j++)
			{
	//	    unsigned int sizeofmanor=Pool_manor[i].size();
			unsigned int tem_mem = Pool_manor[i][j];
			if(tem_mem!= i)
			{
			beacon_belong = 0;			//remeber to reset beacon_belong!!! @_@
			for(unsigned int k=0; k<Pool_manor[tem_mem].size(); k++)
				{
				if(i==Pool_manor[tem_mem][k])
				{
				beacon_belong = 1;		//make sure i is member of tem_men's manor!
				}
				}
			if(beacon_belong == 0)
				{
				//tem_mem's manor doesn't include node i, so this tem_mem should be elimiated from node i's manor!
				//Pool_manor[i] - tem_mem
				std::vector<unsigned int>::iterator where = std::find(tem_manor.begin(), tem_manor.end(), tem_mem);
				tem_manor.erase(where);	//
				}
			}
			}
		Pool_manor[i] = tem_manor;
		}// finalizing membership is actually finished, but there is another task: clarifying which is cluster head!





#ifdef PRINT_ON_SCREAN
    //do not output to save time
    //---output Pool_manor
    std::cout <<"==========ClusteringPhase III results:=========="<<std::endl;
    for(unsigned int t=0; t < Pool_manor.size(); t++)
	{
	std::cout <<"the manor of node ("<<t<<") is: (";
	for(unsigned int v=0; v < Pool_manor[t].size(); v++ )
	    {
	    std::cout << Pool_manor[t][v]<<"\t";
	    }

	std::cout << ")"<<"\n";
	}
#endif

    //extract only one unique clique!
    //if the index of a manor equals to the first element of that manor,
    //this manor is a original one! choosen!
    //then,
    //output Pool_Groups which is the finally formed clusters.
    Cluster c;
    std::cout<<"the final clustering output is :"<<"\n";
    	    unsigned int counter=0;
    	    for(unsigned int i=0; i < Pool_manor.size(); i++)
    		{
			unsigned int forward=Pool_manor[i].front();
			if(i==forward)
				{
				std::cout<<"cluster ("<< i <<") contains : ";
				c.push_back(Pool_manor[i].front());

				for(unsigned int j=0; j < Pool_manor[i].size(); j++)
				{
				std::cout<<Pool_manor[i][j]<<"\t";
				counter++;
				}
				std::cout<<"\n";
				}
    		}
    	    //form Pool_Groups
    	    for(unsigned int k=0; k < c.size(); k++)
				{
    	    	Pool_Groups.push_back(Pool_manor[c[k]]);
				}

    	    Cluster repeated;
    	    for(unsigned int k=0; k < Pool_Groups.size(); k++)
    		{
		for(unsigned int l=0; l<Pool_Groups[k].size(); l++)
		    {
		    unsigned int index=Pool_Groups[k][l];

		    for(unsigned int m=0; m < Pool_Groups.size(); m++)
			{
			for(unsigned int n=0; n < Pool_Groups[m].size(); n++)
			    {
			    if(index==Pool_Groups[m][n])
				{
				repeated.push_back(index);
				}
			    }
			}
		    }
    		}
    	 for(unsigned int m=0; m < repeated.size(); m++)
    	     {
	     std::cout<<repeated[m]<<"\t";
	     }
    	std::cout<<"\n";
    	std::cout<<"****************"<<"\n";


    	std::cout<<"the number of involved nodes is :"<<counter<<"\n";
        std::cout <<"==========ClusteringPhase III resultes end!=========="<<std::endl;
        std::cout<<std::endl;



#ifdef SURVIVAL_CLUSTERS
//    	std::cout<<"*****SURVIVAL_CLUSTERS**********"<<"\n";

    /*
     * newPU is the number of added PUs every time
     */
    unsigned int oldPU =30;
    unsigned int newPU =10;// number of PRs added each time, the total number is decided as nuwPU*(a number)
    unsigned int ChDim =10;
    unsigned int RadiusPR = 20;
    unsigned int network_geo_size = 100;

    std::ofstream survival_file;
    survival_file.open("survival.txt", 	std::ios::app);



    /*
     * the total number of PUs is controlled here. currently, add PUs for 10 times.
     */
    while(PRGroup.size() < oldPU + newPU*30 ){
	//add new PUs in the network
//        survival_file << "PRGroup.size()="<<PRGroup.size()<<"\n";

	    for (unsigned int j = 0; j < newPU; j++)		//each node node within the possible scope
		{
		Point p(2);
		for (unsigned int i = 0; i < 2; i++)
		    {
		    //points distributed in the plane of 100m X 100m
		    p[i]=network_geo_size *(rand() * (1.0) / RAND_MAX);
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
	    unsigned int num_unclustered=0;
//	    survival << "Pool_Groups.size() = " << Pool_Groups.size() << "\n";

	    for(unsigned int i=0; i< Pool_Groups.size(); i++)
		{
//		survival << "Pool_Groups[" << i << "].size()= " << Pool_Groups[i].size() << "\n";
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
	    survival_file << num_unclustered << "\t";
    }
	survival_file << "\n";
	survival_file.close();

#endif
















//===the final results output:===

        //====//output #CCC and #OCC!
            std::cout<<"\n";
            std::cout<<"number of ccc and outwards cc are as follows: "<<"\n";
            double SumNumCCC=0;
            Cluster set_SumNumCCC;

            double SumNumCCC_allsize = 0;
            Cluster set_SumNumCCC_allsize;

//            double SumNumOutCC=0;
//            Cluster set_SumNumOutCC;

            unsigned int SumOCC_onecluster=0;
            unsigned int SumOCC_onenode=0;

            for(unsigned int i=0;i<Pool_Groups.size();i++)
        	{
        	    std::cout<<"inner ccc and outwards cc of Group("<<Pool_Groups[i].front()<<")= "<<"\t";
        	    unsigned int NumCCC=CaculateNumCCC(Pool_Groups[i]);
        		if(Pool_Groups[i].size()>1)
        			{
        			SumNumCCC += NumCCC; //SumNumCCC is the total number of CCC in the whole network!
        			set_SumNumCCC.push_back(NumCCC);
        			}
        		SumNumCCC_allsize += NumCCC;
        		set_SumNumCCC_allsize.push_back(NumCCC);

        	    //sum of products of size and #CCC
        	    prod_size_ccc += NumCCC * Pool_Groups[i].size();
        	    std::cout<<NumCCC<<"\t";
        	    unsigned int inner=0;

//        	    for(unsigned int j=0; j<Pool_Groups[i].size(); j++)	//Pool_Groups[i][j] is a member of group i, we want to find the number of cc between it and other nodes!
//        		{
//        		Cluster tem_neighbor= Pool_Neighbors[Pool_Groups[i][j]];	//node Pool_Groups[i][j]'s neighborhood
//        		Cluster tem_group = Pool_Groups[i];				//the group where node Pool_Groups[i][j] is
////***		problem!! you have to find the real cluster where node i belong!
//
//        		//tem_tem_neighbor is unchanged, unchanged, unchanged, unchanged!
//        		Cluster tem_tem_neighbor=tem_neighbor;
//
////        		std::vector<unsigned int>::iterator where = std::find(tem_tem_neighbor.begin(), tem_tem_neighbor.end(), Pool_Groups[i][j]);
////        		tem_group.erase(where);
//
//        		//elimilate Pool_Groups[i][j]'s neighbors locating in the same group with Pool_Groups[i][j].
//        		for(unsigned int k=0; k<tem_group.size(); k++)
//        		    {
//        		    for(unsigned m=0; m<tem_tem_neighbor.size(); m++)
//        			{
//        			if (tem_tem_neighbor[m]==tem_group[k])
//        			    {
//        			    std::vector<unsigned int>::iterator where = std::find(tem_neighbor.begin(), tem_neighbor.end(), tem_tem_neighbor[m]);
//        			    tem_neighbor.erase(where);
//        			    }
//        			}
//        		    }
//
//        		tem_tem_neighbor=tem_neighbor;
//        		//elimilate Pool_Groups[i][j]'s neighbors which are clusterhead of other clusters
////			for(unsigned int k=0; k<c.size(); k++)
////			    {
////			    for(unsigned m=0; m<tem_tem_neighbor.size(); m++)
////				{
////				if (tem_tem_neighbor[m]==c[k])
////				    {
////				    if(Pool_Groups[c[k]].size()>1)
////					{
////					std::vector<unsigned int>::iterator where = std::find(tem_neighbor.begin(), tem_neighbor.end(), tem_tem_neighbor[m]);
////					tem_neighbor.erase(where);
////					}
////				    }
////
////				}
////			    }
//
//
////        		for(unsigned int l=0; l< tem_neighbor.size(); l++)
////        		    {
////        		    inner=inner+inner_product(v_all_copy[tem_neighbor[l]].begin(), v_all_copy[tem_neighbor[l]].end(),v_all_copy[Pool_Groups[i][j]].begin(), 0);
////        		    }
//			// For each border node, the #OCC is the number of channels via which outward connection can be made!
//			Channel_onenode v_tem(10,0);
//			Channel_onenode v_onebordernode(10,0);
//
//			//tem_neighbor is the set of border nodes of node Pool_Groups[i][j]!!!
//			for(unsigned int l=0; l< tem_neighbor.size(); l++)
//			    {
//			    for(unsigned int m=0; m< 10; m++)
//				{
//				v_tem[m] = v_all_copy[Pool_Groups[i][j]][m]*v_all_copy[tem_neighbor[l]][m];
//				}
//			    for(unsigned int n=0; n< 10; n++)
//				{
//				v_onebordernode[n]=v_onebordernode[n]||v_tem[n];
//				}
//			    }
//			//all the usable channels on which conncections can be bulit!
//			unsigned int SumOCC_onenode = std::accumulate( v_onebordernode.begin(), v_onebordernode.end(), 0 );
//			//the sum of usable channels in cluster Pool_Groups[i]!
//			SumOCC_onecluster = SumOCC_onecluster + SumOCC_onenode;
//        		}
//        		set_SumNumOutCC.push_back(SumOCC_onecluster);
//
//        		//output #OCC
//        		std::cout<<SumOCC_onecluster<<"\t";
//        		std::cout<<"\n";
//        		SumOCC_onecluster=0;
//        		SumOCC_onenode=0;
//        		inner=0;
        	}

	    for(unsigned int i=0; i<100;i++)
		{
		size_distr.push_back(0);
		CCC_distr.push_back(0);
//		OCC_distr.push_back(0);
		}


	    //---here to record the distribution of #CCC
	    //set_SumNumCCC
	    for(unsigned int i=0; i < set_SumNumCCC.size(); i++)
		{
		CCC_distr[set_SumNumCCC[i]]=CCC_distr[set_SumNumCCC[i]]+1;
		}

            average_CCC = SumNumCCC/set_SumNumCCC.size();
            average_CCC_allsize = SumNumCCC_allsize/set_SumNumCCC_allsize.size();
            std::cout<<"Average number of inner CCC is "<<average_CCC<<"\n";

            //the average usable channels in each cluster!
//            unsigned int SumOCC = std::accumulate( set_SumNumOutCC.begin(), set_SumNumOutCC.end(), 0 );
//            //---here to record the distribution of #OCC
//            for(unsigned int i=0; i < set_SumNumOutCC.size(); i++)
//        	{
//        	OCC_distr[set_SumNumOutCC[i]]=OCC_distr[set_SumNumOutCC[i]]+1;
//        	}
//
//            average_OCC=SumOCC/(Pool_Groups.size()*1.0);
//            std::cout<<"Average number of outwards cc per node is "<<average_OCC<<"\n";

            number_clusters=Pool_Groups.size();
            std::cout<<"number of clusters in this topology is: "<<number_clusters <<std::endl;


	    for(unsigned int i=0;i<Pool_Groups.size();i++)
			{
			set_size.push_back(Pool_Groups[i].size());
			}
	    //---here to record the distribution of cluster size
	    for(unsigned int j=0; j<set_size.size(); j++)
			{
			size_distr[set_size[j]]=size_distr[set_size[j]]+1;
			}

//	    //replace the above 2 blocs
//	    for(unsigned int i=0;i<Pool_Groups.size();i++)
//			{
//	    	unsigned int size= Pool_Groups[i].size();
//			size_distr[size]=size_distr[size]+1;
//			}



	    average_size=std::accumulate(set_size.begin(),set_size.end(),0)/(Pool_Groups.size()*1.0);
	    std::cout<<"average size of clusters is: "<<average_size<<std::endl;
            //========output CV of cluster size
	    cv_size=CaculateCV(set_size);
	    std::cout<<"CV of cluster size is: "<<cv_size<<std::endl;

            //========output CV of CCC and OCC
	    cv_CCC=CaculateCV(set_SumNumCCC);
	    std::cout<<"CV of CCC is: "<<cv_CCC<<std::endl;
//	    cv_OCC=CaculateCV(set_SumNumOutCC);
//	    std::cout<<"CV of OCC is: "<<cv_OCC<<std::endl;
	    std::cout<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << "\n";
//	    set_size.clear();
    }

/*
 * v_all contains initial channel availability before adding new PUs
 * 22.07.2014
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
//	Pool_Groups_copy.clear();

	for(unsigned int i=0; i<Pool_MotherGroups.size(); i++)
	    {
	    Pool_MotherGroups[i].clear();
	    }
	for(unsigned int i=0; i<Pool_MotherGroups_copy.size(); i++)
	    {
	    Pool_MotherGroups_copy[i].clear();
	    }
	for(unsigned int i=0; i<backuppool_original_group.size(); i++)
	    {
	    backuppool_original_group[i].clear();
	    }
	ai_vector.clear();
	bi_vector.clear();
	c.clear();


//	Pool_CompetitiveMotherGroups.clear();	//overlapping node的mother group的集合


	for(unsigned int i=0; i<NodesAfterSearch.size(); i++)
	    {
	    NodesAfterSearch[i]=0;
	    }
//	for(unsigned int i=0; i<NodesSurrender.size(); i++)
//	    {
//	    NodesSurrender[i]=0;
//	    }
	set_surrender.clear();
//	for(unsigned int i=0; i<NodesInPhaseII.size(); i++)
//	    {
//	    NodesInPhaseII[i]=0;
//	    }
	overlapping_nodes.clear();
	real_overlapping_nodes.clear();

//	for(unsigned int i=0; i<Nodefinalized.size(); i++)
//	    {
//	    Nodefinalized[i]=0;
//	    }
//	for(unsigned int i=0; i<Nodefinalized_AsHead.size(); i++)
//	    {
//	    Nodefinalized_AsHead[i]=0;
//	    }

	v_all.clear();
	v_all_copy.clear();
	v_all_PR.clear();
	v_all_PR_copy.clear();
	PRGroup.clear();

	for(unsigned int i=0; i<Pool_manor.size(); i++)
	    {
	    Pool_manor[i].clear();
	    }

	for(unsigned int i=0; i<Pool_manor_copy.size(); i++)
	    {
	    Pool_manor_copy[i].clear();
	    }
	for(unsigned int i=0; i<set_of_y.size(); i++)
	    {
	    set_of_y[i]=0;
	    }
	for(unsigned int i=0; i<set_of_size.size(); i++)
	    {
	    set_of_size[i]=0;
	    }
	for(unsigned int i=0; i<set_Q.size(); i++)
	    {
	    set_Q[i]=0;
	    }


	NumberOfCR=0;
	RadiusOfCR=0;
	//phaseI 的遍历遍数
	int Num=0;



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
