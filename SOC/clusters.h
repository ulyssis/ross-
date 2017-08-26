//#define CLUSTER_SIZE_CONTROL
//#define PRINT_ON_SCREAN

//#define _DOUBLE_CHECK_UNACCURATE_SPECTRUM_SENSING
//#define FALSE_NEGATIVE_SPECTRUM_SENSING 0.1

//#define SURVIVAL_CLUSTERS

#include <vector>
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>
#include <fstream>
#include <sstream>
//#include "distance.h"

namespace Clustering{

// a single point is represented by vector of doubles
	typedef boost::numeric::ublas::vector<double> Point;
//a group of nodes.
	typedef std::vector<Point> Points;
	typedef std::vector<double> doubleVector;

	typedef unsigned int ClusterId;
	typedef unsigned int PointId;	
//ChannelId
	typedef unsigned int ChannelId;

// a cluster is a set of pointid
	typedef std::vector<PointId> Cluster;
// a Neighbors is a set of pointid
	typedef std::vector<PointId> Neighbors;

	typedef std::vector<Neighbors> PoolofNeighborhoods;
	typedef std::vector<Cluster> PoolofClusters;

// several ints represent a set of channels avaible on one node.
	typedef std::vector<unsigned int> Channel_onenode;
	typedef std::vector<Channel_onenode> Channel_allnode;
	typedef std::vector<unsigned int> robustness_vector;

	  //ai_vector and bi_vector store caculated values of ai and bi



	void randomInit	(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR, unsigned int network_geo_size);

	class Clusters
	{
	public:
		Clusters (Points & ps) : _ps(ps)	//constructor function;
							//Initialization Lists;
							//set _ps to ps;
		{
//			_pointId_to_clusterId.resize(_ps.size(), 0);
			NumOfNodes=_ps.size();
			//size_distr initialization!
			    for(unsigned int k=0; k<100; k++)
				{
				size_distr.push_back(0);
				}
		};

		// assign each point to a new cluster
		//7-12, this a empty function, the initial purpose should be distributing the nodes.
//		void uniformPartition();

		void computeDistance()
		{

			unsigned int size = _ps.size();
			_dis.resize(size, size, false);			//formulate vector to a matrix of size X size
			for (unsigned int i=0; i < size; i++)
			{
				for (unsigned int j=i+1; j < size; j++)
				{
					Point x =_ps[i];
					Point y =_ps[j];
					double d = sqrt((y[0]-x[0])*(y[0]-x[0])+(y[1]-x[1])*(y[1]-x[1]));
					//return d;

				    _dis(j, i) = _dis(i, j) = d;
//				    std::cout << "(" << i << ", " << j << ")=" << d << " ";
				}
//				std::cout << std::endl;
			}
		std::cout<<"\n"<<"void computeDistance() is OK"<<"\n"<<std::endl;
		};

		void channelRandomInit(unsigned int ChDim, double RadiusPR, double RadiusCR);//this function initializes the channels on nodes, and caculate the value of ai and bi on each node.
		unsigned int survival_check(Cluster cluster, unsigned int RadiusCR, unsigned int ChDim);
		void cluster_with_accurate_spectrum(Cluster cluster, unsigned int index, float false_positive, unsigned int ChDim);


	//	void GetRobustness();

		Neighbors findNeighbors(PointId pid, double threshold);

		// cacalate the number of common control channel in a group of nodes
		unsigned int CaculateNumCCC(Cluster group);
		double CaculateCV(Cluster group);
		double CaculateCV(doubleVector doubleVector); //vector of double
		void MotherGroups_Update();

		//void ClusteringPhaseI(unsigned int num_nodes);
		void ClusteringPhaseI();
		void ClusteringPhaseI(unsigned int optimalSize);


//		check whether every node is rectuited into at least one cluster
		unsigned int CheckVisitComplete();

//		output clusters after phaseI
		void OutputClusters();

		void ClusteringPhaseII();
		void ClusteringPhaseIII();

		double RadiusPR;
		double RadiusCR;

		void clearall();

		float average_CCC;
		float average_CCC_allsize;
		float average_OCC;
		float number_clusters;
		float average_size;
		float cv_size;
		float cv_CCC;
		float cv_OCC;
		std::vector<float> size_distr;
		std::vector<float> CCC_distr;
		std::vector<float> OCC_distr;
		double prod_size_ccc;

	protected:
		// the collection of points we are working on
		Points& _ps;
		
		// the scale of _ps
		unsigned int NumOfNodes;

		// mapping point_id -> clusterId
		std::vector<ClusterId> _pointId_to_clusterId;

		// the collection of clusters
		std::vector<Cluster> _clusters;

		// simarity_matrix
		boost::numeric::ublas::matrix<double> _dis;



		friend	
			std::ostream& operator << (std::ostream& o, const Clusters& c);
	};

}
