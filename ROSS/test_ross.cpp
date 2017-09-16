/*
 * feature of size control can be changed by modifying macro in cluster.h
 */



#include "distance.h"
#include <iostream>
#include "clusters.h"
#include <time.h>

int main()
    {

	for (unsigned int loop =1; loop<2; loop++){
		std::ofstream myfile;
		std::string fileName;
		std::string seed_index = SSTR(loop);
		fileName = "result_" + seed_index + ".txt";
		const char *str = fileName.c_str();
//		std::ifstream myfile(str);
		myfile.open (str);

		using namespace Clustering;

		Clustering::Points ps;
		Cluster vector_numOverlappingNode;
		doubleVector vector_prod_size_ccc;
		float sum_prod_size_ccc =0;
		doubleVector vector_average_CCC;
		doubleVector vector_average_CCC_allsize;
		doubleVector vector_average_OCC;
		doubleVector vector_average_size;
		doubleVector vector_average_size_allsize;
		doubleVector vector_average_ClaimingCluster;
		doubleVector vector_number_clusters;
		Cluster vector_UpdateRound;

		float sum_average_CCC=0;
		float sum_average_CCC_allsize=0;
		float sum_average_OCC=0;
		float sum_number_clusters=0;
		float sum_average_size=0;
		float sum_cv_CCC=0;
		float sum_cv_OCC=0;
		float sum_average_neighborSize=0;
		float sum_cv_neighborSize=0;
		float sum_numOverlappingNode =0;
		float sum_average_ClaimingCluster =0;
		float sum_average_num_channel_pernode =0;
		float sum_UpdateRound =0;
		std::vector<float> distr(200,0);
		std::vector<float> accum_numCCC(200,0);
		std::vector<float> accum_numOCC(200,0);
		std::vector<float> clustersize_cdf;
		std::vector<float> numCCC_cdf;
		std::vector<float> numOCC_cdf;

		//store data in every run to cacalate the confidence interval
		std::vector<float> container;

		unsigned int numLoops=50;
		unsigned int numCR =loop*100;
		for(unsigned int seed=0; seed<numLoops; seed++)
			{

			srand(seed);
			//void randomInit	(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR, unsigned int network_geo_size)

			/*
			 * when macros 'SURVIVAL_CLUSTERS' or 'SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME' are used, remeber to set network_geo_size as the same
			 * with the last parameter in the following function.
			 * unsigned int network_geo_size = 100;
			 * in clusters.cpp
			 */
			Clustering::randomInit(ps, 2, numCR, 30, 50);
//			Clustering::randomInit(ps, 2, numCR, 10, 30);
//

			Clustering::Clusters clusters(ps); //private member _ps is initialized here!
			clusters.computeDistance(); //private member _dis is initialized here!

			/*
			 * void channelRandomInit(unsigned int ChDim, double RadiusPR, double RadiusCR);
			 * number of channel; broadcast radius of PR; broadcast radius of CR;
			 */

			Channel_allnode v_all;
			v_all = clusters.channelRandomInit(10,20,10);
//			v_all = clusters.channelRandomInit(10,10,10);

	#ifdef PRINT_OUT_Q_MATRIX
			clusters.legitimateClusters(ps, numCR, seed, v_all);
	#endif

	#ifdef CLUSTER_SIZE_CONTROL
			clusters.ClusteringPhaseI(6);
	#endif


	#ifndef CLUSTER_SIZE_CONTROL
			clusters.ClusteringPhaseI();
	#endif
			clusters.OutputClustersAfterPhaseI();
#ifndef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME
			clusters.ClusteringPhaseII();
#endif
#ifdef SURVIVAL_CLUSTERS_TEST_CENTRALIZED_SCHEME
    clusters.ClusteringPhaseII(seed);
#endif

			//accumulation!
			sum_average_CCC = sum_average_CCC + clusters.average_CCC;
			sum_average_CCC_allsize = sum_average_CCC_allsize + clusters.average_CCC_allsize;
			sum_average_OCC = sum_average_OCC + clusters.average_OCC;
	//	    //to caculate standard deviation of 50 average_OCC(s)
	//	    container.push_back(clusters.average_OCC);
			sum_number_clusters= sum_number_clusters + clusters.number_clusters;
			vector_number_clusters.push_back(clusters.number_clusters);
			sum_average_size=sum_average_size+clusters.average_size;

			sum_average_neighborSize = sum_average_neighborSize+clusters.average_neighborSize;
			sum_cv_neighborSize = sum_cv_neighborSize+clusters.cv_neighborSize;

			sum_numOverlappingNode+=clusters.numOverlappingNode;
#ifdef GREEDY
			sum_UpdateRound+=clusters.UpdateRound;
#endif

			vector_numOverlappingNode.push_back(clusters.numOverlappingNode/numCR);
			vector_prod_size_ccc.push_back(clusters.prod_size_ccc);
			sum_prod_size_ccc += clusters.prod_size_ccc;
			clusters.prod_size_ccc =0; //clear to 0

			vector_average_CCC.push_back(clusters.average_CCC);
			vector_average_CCC_allsize.push_back(clusters.average_CCC_allsize);
			vector_average_OCC.push_back(clusters.average_OCC);
			vector_average_size.push_back(clusters.average_size);
			vector_average_size_allsize.push_back(clusters.average_size_allsize);
			sum_average_ClaimingCluster+= clusters.average_ClaimingCluster;
			sum_average_num_channel_pernode += clusters.average_num_channel_pernode;
			vector_average_ClaimingCluster.push_back(clusters.average_ClaimingCluster);
#ifdef GREEDY
			vector_UpdateRound.push_back(clusters.UpdateRound);
#endif
			for(unsigned int i=0; i<distr.size(); i++)
			{
			distr[i]=distr[i]+clusters.size_distr[i];
			}

			for(unsigned int i=0; i<accum_numCCC.size(); i++)
			{
			accum_numCCC[i]=accum_numCCC[i]+clusters.CCC_distr[i];
			}

			for(unsigned int i=0; i<accum_numOCC.size(); i++)
			{
			accum_numOCC[i]=accum_numOCC[i]+clusters.OCC_distr[i];
			}

			clusters.clearall();
			//clear ps
			ps.clear();
			std::cout<<"seed is "<< seed <<"\n";

			}
	//	vector_numOverlappingNode.clear();
		Clustering::Clusters temp(ps);


		std::cout<<"\n final output of ROSS is:"<<std::endl;
		std::cout<<"average_CCC(size>1) = "<<sum_average_CCC/numLoops<<"\n"; // \sum_{num_channels}/num_all_clusters
		std::cout<<"average_CCC(size>=1) = "<<sum_average_CCC_allsize/numLoops<<"\n";
		std::cout<<"average_OCC = "<<sum_average_OCC/numLoops<<"\n";
		std::cout<<"number_clusters = "<<sum_number_clusters/numLoops<<"\n";
		std::cout<<"CI of number_clusters = "<<temp.CaculateCV(vector_number_clusters)<<"\n";
		std::cout<<"average_size = "<<sum_average_size/numLoops<<"\n";
		std::cout<<"average number of channels per node = "<< sum_average_num_channel_pernode/numLoops<<"\n";
		std::cout<<"CI of average cluster size(size>1) = "<<temp.CaculateCV(vector_average_size) <<"\n";
		std::cout<<"CI of average cluster size(size>=1) = "<<temp.CaculateCV(vector_average_size_allsize) <<"\n";
		std::cout<<"CI of average_CCC (size>1)" << temp.CaculateCV(vector_average_CCC) <<"\n";
		std::cout<<"CI of average_CCC (size>=1)" << temp.CaculateCV(vector_average_CCC_allsize) <<"\n";
		std::cout<<"CI of average_OCC = "<< temp.CaculateCV(vector_average_OCC) <<"\n";
		std::cout<<"Average of products of size and number of ICC in the network is " << sum_prod_size_ccc/numLoops <<"\n";
		std::cout<<"CI of products of size and number of ICC in the network is " << temp.CaculateCV(vector_prod_size_ccc) <<"\n";
		std::cout<<"average number of neighbors over 50 runs =" <<sum_average_neighborSize/numLoops <<"\n";
		std::cout<<"average percentage of overlapping nodes is = " << sum_numOverlappingNode/numCR/numLoops << "\n";
		std::cout<<"CI of percentage of overlapping nodes is " << temp.CaculateCV(vector_numOverlappingNode) << "\n";
		std::cout<<"Average Number of claiming clusters is = " << sum_average_ClaimingCluster/numCR*1.0 << "\n";
		std::cout<< "CI of average number of claiming clusters is " << temp.CaculateCV(vector_average_ClaimingCluster) << "\n";
		std::cout<< "Number of control messages in greedy phase is " << sum_UpdateRound/numLoops << "\n\n";


		std::cout<<"---------only data bellow-----:"<< "\n";
		std::cout<<sum_average_CCC/numLoops<<"\n";
		std::cout<<sum_average_CCC_allsize/numLoops<<"\n";
		std::cout<<sum_average_OCC/numLoops<<"\n";
		std::cout<<sum_number_clusters/numLoops<<"\n";
		std::cout<<temp.CaculateCV(vector_number_clusters)<<"\n";
		std::cout<<sum_average_size/numLoops<<"\n";
		std::cout<<sum_average_num_channel_pernode/numLoops<<"\n";
		std::cout<<temp.CaculateCV(vector_average_size)<<"\n";
		std::cout<<temp.CaculateCV(vector_average_size_allsize)<<"\n";
		std::cout<<temp.CaculateCV(vector_average_CCC)<<"\n";
		std::cout<<temp.CaculateCV(vector_average_CCC_allsize) <<"\n";
		std::cout<<temp.CaculateCV(vector_average_OCC)<<"\n";
		std::cout<< sum_prod_size_ccc/numLoops << "\n";
		std::cout<< temp.CaculateCV(vector_prod_size_ccc) << "\n";
		std::cout<< sum_average_neighborSize/numLoops <<"\n";
		std::cout<< sum_numOverlappingNode/numCR/numLoops <<"\n";
		std::cout<< temp.CaculateCV(vector_numOverlappingNode) << "\n";
		std::cout<< sum_average_ClaimingCluster/numLoops << "\n";
		std::cout<< temp.CaculateCV(vector_average_ClaimingCluster) << "\n";
		std::cout<< sum_UpdateRound/numLoops << "\n";


		//---------
		myfile<< numCR <<":\n";
		myfile<<sum_average_CCC/numLoops<<"\n";
		myfile<<sum_average_CCC_allsize/numLoops<<"\n";
		myfile<<sum_average_OCC/numLoops<<"\n";
		myfile<<sum_number_clusters/numLoops<<"\n";
		myfile<<temp.CaculateCV(vector_number_clusters)<<"\n";
		myfile<<sum_average_size/numLoops<<"\n";
		myfile<<sum_average_num_channel_pernode/numLoops<<"\n";
		myfile<<temp.CaculateCV(vector_average_size)<<"\n";
		myfile<<temp.CaculateCV(vector_average_size_allsize)<<"\n";
		myfile<<temp.CaculateCV(vector_average_CCC)<<"\n";
		myfile<<temp.CaculateCV(vector_average_CCC_allsize) <<"\n";
		myfile<<temp.CaculateCV(vector_average_OCC)<<"\n";
		myfile<< sum_prod_size_ccc/numLoops << "\n";
		myfile<< temp.CaculateCV(vector_prod_size_ccc) << "\n";
		myfile<< sum_average_neighborSize/numLoops <<"\n";
		myfile<< sum_numOverlappingNode/numCR/numLoops <<"\n";
		myfile<< temp.CaculateCV(vector_numOverlappingNode) << "\n";
		myfile<< sum_average_ClaimingCluster/numLoops << "\n";
		myfile<< temp.CaculateCV(vector_average_ClaimingCluster) << "\n";
		myfile<< sum_UpdateRound/numLoops << "\n";
		myfile<< "\n";

		//-----------


		temp.clearall();

		//	CDF of size of clusters
			for(unsigned int i=0; i<distr.size(); i++)
			{
			distr[i]=distr[i]/numLoops;
			}

			unsigned int j=distr.size();
			while(!distr[j-1])
				{
			if(distr.back()==0.0)
				{
				distr.pop_back();
				}
			j--;
			}
			//output the number of nodes grouped into clusters with certain size
			std::cout<< "//output the number of nodes grouped into clusters with certain size"<< "\n";
			for(unsigned int k=0; k<distr.size(); k++)
			{
			std::cout<<distr[k]*k<<"\n";
			myfile<<distr[k]*k<<"\n";
			}
			myfile<< "\n";
			std::cout<<std::endl;

//			//output cdf of number of nodes according to clusters size
//			clustersize_cdf.push_back(0);
//			for(unsigned int k=0; k<distr.size(); k++)
//			{
//			float pre_cum=clustersize_cdf.back();
//			clustersize_cdf.push_back(pre_cum+(distr[k+1]*(k+1)));
//			std::cout<<clustersize_cdf[k]<<"\n";
//			}
//			std::cout<<std::endl;


			//output standard deviation of 50 average_OCC(s)
	//	    for(unsigned int l=0; l<container.size(); l++)
	//		{
	//		std::cout<<container[l]<<",";
	//		}

			std::cout<<std::endl;

	////	CDF of number of CCC
	//	    for(unsigned int i=0; i<accum_numCCC.size(); i++)
	//		{
	//		accum_numCCC[i]=accum_numCCC[i]/50;
	//		}
	//
	//	    unsigned int jj=accum_numCCC.size();
	//	    while(!accum_numCCC[jj-1])
	//		{
	//		if(accum_numCCC.back()==0.0)
	//		    {
	//		    accum_numCCC.pop_back();
	//		    }
	//		jj--;
	//		}
	//	    float firstnum=accum_numCCC[0];
	//	    numCCC_cdf.push_back(firstnum);
	//
	//	    std::cout<<"there output the series of accumulative numbers: "<<"\n";
	//	    for(unsigned int k=0; k<accum_numCCC.size(); k++)
	//		{
	//		float pre_cum=numCCC_cdf.back();
	//		numCCC_cdf.push_back(pre_cum+(accum_numCCC[k+1]));
	//
	//		std::cout<<numCCC_cdf[k]<<"\n";
	//		}
	////	CDF of number of OCC
	//	    for(unsigned int i=0; i<accum_numOCC.size(); i++)
	//		{
	//		accum_numOCC[i]=accum_numOCC[i]/50;
	//		}
	//
	//	    unsigned int jjj=accum_numOCC.size();
	//	    while(!accum_numOCC[jjj-1])
	//		{
	//		if(accum_numOCC.back()==0.0)
	//		    {
	//		    accum_numOCC.pop_back();
	//		    }
	//		jjj--;
	//		}
	//
	//	    float firstnum2=accum_numOCC[0];
	//	    numOCC_cdf.push_back(firstnum2);
	//
	//		std::cout<<"there output the series of accumulative numbers: "<<"\n";
	//	    for(unsigned int k=0; k<accum_numOCC.size(); k++)
	//		{
	//		float pre_cum=numOCC_cdf.back();
	//		numOCC_cdf.push_back(pre_cum+(accum_numOCC[k+1]));
	//
	//		std::cout<<numOCC_cdf[k]<<"\n";
	//		}
			//clear containers!
			myfile.close();
		}

}

