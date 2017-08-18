#include "distance.h"
#include <iostream>
#include "clusters.h"
#include<time.h>

int main()
    {

 	using namespace Clustering;

	Clustering::Points ps;
	Clustering::Channel_allnode v_all;

	float sum_average_CCC=0;
	float sum_average_OCC=0;
	float sum_number_clusters=0;
	float sum_average_size=0;
	float sum_cv_size=0;
	float sum_cv_CCC=0;
	float sum_cv_OCC=0;
	unsigned int SumUpdateRound=0;
	std::vector<float> distr(100,0);
	std::vector<float> accum_numCCC(100,0);
	std::vector<float> accum_numOCC(100,0);
	std::vector<float> clustersize_cdf;
	std::vector<float> numCCC_cdf;
	std::vector<float> numOCC_cdf;


	//store data in every run to cacalate the confidence interval
	std::vector<float> container;

	for(unsigned int seed=0; seed<50; seed++)
	    {
	    srand(seed);
	    //void randomInit	(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR)
//	    Clustering::randomInit(ps, 2, 100, 70);
	    Clustering::randomInit(ps, 2, 100, 30);

	    Clustering::Clusters clusters(ps); //private member _ps is initialized here!
	    clusters.computeDistance(); //private member _dis is initialized here!

	    //void channelRandomInit(unsigned int ChDim, double RadiusPR, double RadiusCR);
	    clusters.channelRandomInit(10,10,10);


	    clusters.ClusteringPhaseI();
	    clusters.OutputClustersAfterPhaseI();
	    clusters.ClusteringPhaseII();

	    //accumulation!
	    sum_average_CCC = sum_average_CCC + clusters.average_CCC;
	    sum_average_OCC = sum_average_OCC + clusters.average_OCC;
//	    //to caculate standard deviation of 50 average_OCC(s)
//	    container.push_back(clusters.average_OCC);
	    sum_number_clusters= sum_number_clusters + clusters.number_clusters;
	    sum_average_size=sum_average_size+clusters.average_size;
	    sum_cv_size=sum_cv_size+clusters.cv_size;
	    sum_cv_CCC=sum_cv_CCC+clusters.cv_CCC;
	    sum_cv_OCC=sum_cv_OCC+clusters.cv_OCC;

	    SumUpdateRound=SumUpdateRound+clusters.UpdateRound;

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

	    }

//	std::cout<<"final output of ROSS is:"<<std::endl;
//	std::cout<<"average_OCC = "<<sum_average_CCC/50<<"\n";
//	std::cout<<"average_CCC = "<<sum_average_OCC/50<<"\n";
//	std::cout<<"number_clusters = "<<sum_number_clusters/50<<"\n";
//	std::cout<<"average_size = "<<sum_average_size/50<<"\n";
//	std::cout<<"cv_size = "<<sum_cv_size/50<<"\n";
//	std::cout<<"cv_CCC = "<<sum_cv_CCC/50<<"\n";
//	std::cout<<"cv_OCC = "<<sum_cv_OCC/50<<"\n";

	std::cout<<std::endl;
	std::cout<<sum_average_CCC/50<<"\n";
	std::cout<<sum_average_OCC/50<<"\n";
	std::cout<<sum_number_clusters/50<<"\n";
	std::cout<<sum_average_size/50<<"\n";
	std::cout<<sum_cv_size/50<<"\n";
	std::cout<<sum_cv_CCC/50<<"\n";
	std::cout<<sum_cv_OCC/50<<"\n";
	std::cout<<SumUpdateRound/50<<"\n";


	//	CDF of size of clusters
	    for(unsigned int i=0; i<distr.size(); i++)
		{
		distr[i]=distr[i]/50;
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
	    for(unsigned int k=0; k<distr.size(); k++)
		{
		std::cout<<distr[k]*k<<"\t";
		}
	    std::cout<<std::endl;

	    //output cdf of number of nodes according to clusters size
	    clustersize_cdf.push_back(0);
	    for(unsigned int k=0; k<distr.size(); k++)
		{
		float pre_cum=clustersize_cdf.back();
		clustersize_cdf.push_back(pre_cum+(distr[k+1]*(k+1)));
		std::cout<<clustersize_cdf[k]<<"\n";
		}
	    std::cout<<std::endl;
	    //output standard deviation of 50 average_OCC(s)
//	    for(unsigned int l=0; l<container.size(); l++)
//		{
//		std::cout<<container[l]<<",";
//		}

	    std::cout<<std::endl;

//	CDF of number of CCC
	    for(unsigned int i=0; i<accum_numCCC.size(); i++)
		{
		accum_numCCC[i]=accum_numCCC[i]/50;
		}

	    unsigned int jj=accum_numCCC.size();
	    while(!accum_numCCC[jj-1])
		{
		if(accum_numCCC.back()==0.0)
		    {
		    accum_numCCC.pop_back();
		    }
		jj--;
		}
	    float firstnum=accum_numCCC[0];
	    numCCC_cdf.push_back(firstnum);

	    std::cout<<"there output the series of accumulative numbers: "<<"\n";
	    for(unsigned int k=0; k<accum_numCCC.size(); k++)
		{
		float pre_cum=numCCC_cdf.back();
		numCCC_cdf.push_back(pre_cum+(accum_numCCC[k+1]));

		std::cout<<numCCC_cdf[k]<<"\n";
		}
//	CDF of number of OCC
	    for(unsigned int i=0; i<accum_numOCC.size(); i++)
		{
		accum_numOCC[i]=accum_numOCC[i]/50;
		}

	    unsigned int jjj=accum_numOCC.size();
	    while(!accum_numOCC[jjj-1])
		{
		if(accum_numOCC.back()==0.0)
		    {
		    accum_numOCC.pop_back();
		    }
		jjj--;
		}

	    float firstnum2=accum_numOCC[0];
	    numOCC_cdf.push_back(firstnum2);

		std::cout<<"there output the series of accumulative numbers: "<<"\n";
	    for(unsigned int k=0; k<accum_numOCC.size(); k++)
		{
		float pre_cum=numOCC_cdf.back();
		numOCC_cdf.push_back(pre_cum+(accum_numOCC[k+1]));

		std::cout<<numOCC_cdf[k]<<"\n";
		}
	    //clear containers!

}

