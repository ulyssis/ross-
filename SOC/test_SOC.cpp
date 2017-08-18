/*
 * remember to modify the numCR in cluster.cpp file
 */

#include "soc.h"
#include "distance.h"
#include <iostream>



int main()
    {

 	using namespace Clustering;

	Clustering::Points ps;
	Clustering::Channel_allnode v_all;

	doubleVector vector_prod_size_ccc;
	float sum_prod_size_ccc =0;
	doubleVector vector_average_CCC;
	doubleVector vector_average_CCC_allsize;
	doubleVector vector_average_OCC;
	doubleVector vector_average_size;



	float sum_average_CCC=0;
	float sum_average_CCC_allsize =0;
	float sum_average_OCC=0;
	float sum_number_clusters=0;
	float sum_average_size=0;
	float sum_cv_size=0;
	float sum_cv_CCC=0;
	float sum_cv_OCC=0;
	std::vector<float> distr(100,0);
	std::vector<float> accum_numCCC(100);
	std::vector<float> accum_numOCC(100,0);
	std::vector<float> clustersize_cdf;
	std::vector<float> numCCC_cdf;
	std::vector<float> numOCC_cdf;

	unsigned int numLoops=50;
	unsigned int numCR =100;
	for(unsigned int seed=0; seed < numLoops; seed++)
	    {
	    srand(seed);
	    //void randomInit(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR, unsigned int network_geo_size)
	    Clustering::randomInit(ps, 2, numCR, 30, 50);
//	    Clustering::randomInit(ps, 2, numCR, 10, 30);

	    Clustering::Clusters clusters(ps); //private member _ps is initialized here!
	    clusters.computeDistance(); //private member _dis is initialized here!

	    //void channelRandomInit(unsigned int ChDim, double RadiusPR, double RadiusCR);
	    clusters.channelRandomInit(10,20,10);
//	    clusters.channelRandomInit(10,10,10);

	    // Maximum edge biclique computation
#ifdef CLUSTER_SIZE_CONTROL
	    clusters.ClusteringPhaseI(2);
#endif
#ifndef CLUSTER_SIZE_CONTROL
	    clusters.ClusteringPhaseI();
#endif


	    clusters.ClusteringPhaseII();

	    clusters.ClusteringPhaseIII();

	    sum_average_CCC = sum_average_CCC + clusters.average_CCC;
	    sum_average_CCC_allsize = sum_average_CCC_allsize + clusters.average_CCC_allsize;
	    sum_average_OCC = sum_average_OCC + clusters.average_OCC;
	    sum_number_clusters= sum_number_clusters + clusters.number_clusters;
	    sum_average_size=sum_average_size+clusters.average_size;
	    sum_cv_size=sum_cv_size+clusters.cv_size;
	    sum_cv_CCC=sum_cv_CCC+clusters.cv_CCC;
	    sum_cv_OCC=sum_cv_OCC+clusters.cv_OCC;

	    vector_average_CCC.push_back(clusters.average_CCC);
	    vector_average_CCC_allsize.push_back(clusters.average_CCC_allsize);
	    vector_average_OCC.push_back(clusters.average_OCC);
	    vector_average_size.push_back(clusters.average_size);
	    vector_prod_size_ccc.push_back(clusters.prod_size_ccc);
	    sum_prod_size_ccc += clusters.prod_size_ccc;
	    clusters.prod_size_ccc =0; //clear to 0

	    //store the number of clusters with certain sizes!
	    for(unsigned int i=0; i<distr.size(); i++)
		{
		distr[i]=distr[i]+clusters.size_distr[i];
		}

	    //store the number of clusters with certain number of CCC!
	    for(unsigned int i=0; i<accum_numCCC.size(); i++)
		{
		accum_numCCC[i]=accum_numCCC[i]+clusters.CCC_distr[i];
		}

	    //store the number of clusters with certain number of OCC!
//	    for(unsigned int i=0; i<accum_numOCC.size(); i++)
//		{
//		accum_numOCC[i]=accum_numOCC[i]+clusters.OCC_distr[i];
//		}

	    clusters.clearall();
	    //clear ps
	    ps.clear();

	    }


	Clustering::Clusters temp(ps);

	std::cout<<"\n the averaged values of numLoops experiments are:"<<std::endl;
	std::cout<<"average_CCC(size>1) = "<<sum_average_CCC/numLoops<<"\n";
	std::cout<<"average_CCC(size>=1) = "<<sum_average_CCC_allsize/numLoops<<"\n";
	std::cout<<sum_average_OCC/numLoops<<"\n";
	std::cout<<"number_clusters = "<<sum_number_clusters/numLoops<<"\n";
	std::cout<<"\n";
	std::cout<<sum_average_size/numLoops<<"\n";
	std::cout<<"\n";
	std::cout<<temp.CaculateCV(vector_average_size)<<"\n";
	std::cout<<"\n";
	std::cout<<temp.CaculateCV(vector_average_CCC)<<"\n";
	std::cout<<temp.CaculateCV(vector_average_CCC_allsize)<<"\n";//......
	std::cout<<temp.CaculateCV(vector_average_OCC)<<"\n";
	std::cout<<sum_prod_size_ccc/numLoops <<"\n"; //"Average of products of size and number of ICC in the network is " <<
	std::cout<< temp.CaculateCV(vector_prod_size_ccc) <<"\n";

	temp.clearall();

//	std::cout<<"***************************"<<std::endl;
////	std::cout<<"final output of ROSS is:"<<std::endl;
////	std::cout<<"average_OCC = "<<sum_average_CCC/numLoops<<"\n";
////	std::cout<<"average_CCC = "<<sum_average_OCC/numLoops<<"\n";
////	std::cout<<"number_clusters = "<<sum_number_clusters/numLoops<<"\n";
////	std::cout<<"average_size = "<<sum_average_size/numLoops<<"\n";
////	std::cout<<"cv_size = "<<sum_cv_size/numLoops<<"\n";
////	std::cout<<"cv_CCC = "<<sum_cv_CCC/numLoops<<"\n";
////	std::cout<<"cv_OCC = "<<sum_cv_OCC/numLoops<<"\n";
//
//	std::cout<<std::endl;
//	std::cout<<std::endl;
//	std::cout<<sum_average_CCC/numLoops<<"\n";
//	std::cout<<sum_average_OCC/numLoops<<"\n";
//	std::cout<<sum_number_clusters/numLoops<<"\n";
//	std::cout<<sum_average_size/numLoops<<"\n";
//	std::cout<<sum_cv_size/numLoops<<"\n";
//	std::cout<<sum_cv_CCC/numLoops<<"\n";
//	std::cout<<sum_cv_OCC/numLoops<<"\n";

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
	    std::cout<< "xxxxxx  distr.size() = "<< distr.size() << "\n";

	    //output the number of nodes grouped into clusters with certain size
	    std::cout<< "//output the number of nodes grouped into clusters with certain size"<< "\n";
	    for(unsigned int k=0; k<distr.size(); k++)
		{
		std::cout<<distr[k]*k<<"\n";
		}
	    std::cout<<std::endl;

//	    //output cdf of number of nodes according to clusters size
//	    clustersize_cdf.push_back(0);
//	    for(unsigned int k=0; k<distr.size(); k++)
//		{
//		float pre_cum=clustersize_cdf.back();
//		clustersize_cdf.push_back(pre_cum+(distr[k+1]*(k+1)));
//		std::cout<<clustersize_cdf[k]<<"\n";
//		}
//	    std::cout<<std::endl;

//	    //	CDF of number of CCC
//	    	    for(unsigned int i=0; i<accum_numCCC.size(); i++)
//	    		{
//	    		accum_numCCC[i]=accum_numCCC[i]/numLoops;
//	    		}
//
//	    	    unsigned int jj=accum_numCCC.size();
//	    	    while(!accum_numCCC[jj-1])
//	    		{
//	    		if(accum_numCCC.back()==0.0)
//	    		    {
//	    		    accum_numCCC.pop_back();
//	    		    }
//	    		jj--;
//	    		}
//	    	    float firstnum=accum_numCCC[0];
//	    	    numCCC_cdf.push_back(firstnum);
//
//	    	    std::cout<<"there output the series of accumulative numbers of #CCC: "<<"\n";
//	    	    for(unsigned int k=0; k<accum_numCCC.size(); k++)
//	    		{
//	    		float pre_cum=numCCC_cdf.back();
//	    		numCCC_cdf.push_back(pre_cum+(accum_numCCC[k+1]));
//
//	    		std::cout<<numCCC_cdf[k]<<"\n";
//	    		}
//	    //	CDF of number of OCC
//	    	    for(unsigned int i=0; i<accum_numOCC.size(); i++)
//	    		{
//	    		accum_numOCC[i]=accum_numOCC[i]/numLoops;
//	    		}
//
//	    	    unsigned int jjj=accum_numOCC.size();
//	    	    while(!accum_numOCC[jjj-1])
//	    		{
//	    		if(accum_numOCC.back()==0.0)
//	    		    {
//	    		    accum_numOCC.pop_back();
//	    		    }
//	    		jjj--;
//	    		}
//
//	    	    float firstnum2=accum_numOCC[0];
//	    	    numOCC_cdf.push_back(firstnum2);
//
//	    		std::cout<<"there output the series of accumulative numbers of #OCC: "<<"\n";
//	    	    for(unsigned int k=0; k<accum_numOCC.size(); k++)
//	    		{
//	    		float pre_cum=numOCC_cdf.back();
//	    		numOCC_cdf.push_back(pre_cum+(accum_numOCC[k+1]));
//
//	    		std::cout<<numOCC_cdf[k]<<"\n";
//	    		}
}
