#include <vector>
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/foreach.hpp>
#include <fstream>
#include <sstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

namespace Clustering{

	// a single point is represented by vector of doubles
	typedef boost::numeric::ublas::vector<double> Point;
	//a group of nodes.
	typedef std::vector<Point> Points;

	typedef unsigned int ClusterId;
	typedef unsigned int PointId;	
        //ChannelId
	typedef unsigned int ChannelId;

	// a cluster is a set of pointid
	typedef std::vector<PointId> Cluster;
	// a Neighbors is a set of pointid
	typedef std::vector<PointId> Neighbors;

	typedef std::vector<Neighbors> PoolofNeighborhoods;

	// several ints represent a set of channels avaible on one node.
	typedef std::vector<unsigned int> Channel_onenode;
	typedef std::vector<Channel_onenode> Channel_allnode;
	typedef std::vector<unsigned int> robustness_vector;

	  //ai_vector and bi_vector store caculated values of ai and bi



	void randomInit	(Points & ps, unsigned int GeoDim, unsigned int NumCR, unsigned int NumPR);

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
		void legitimateClusters(Points & ps, unsigned int clusterSize);
		void legitimateClusters(Points & ps, unsigned int numCR, unsigned int seed);
		void findNeighbors(double threshold);

		// cacalate the number of common control channel in a group of nodes
		unsigned int checkUniqueness(Cluster group);
		unsigned int inOneCluster(Cluster group);
		unsigned int CaculateNumCCC(Cluster group);
		double CaculateCV(Cluster group);
		//void ClusteringPhaseI(unsigned int num_nodes);
		void ClusteringPhaseI();

//		check whether every node is rectuited into at least one cluster
		unsigned int CheckVisitComplete();

//		output clusters after phaseI
		void OutputClustersAfterPhaseI();

		void ClusteringPhaseII();
		void clearall();

		float average_CCC;
		float average_OCC;
		float number_clusters;
		float average_size;
		float cv_size;
		float cv_CCC;
		float cv_OCC;
		std::vector<float> size_distr;
		std::vector<float> CCC_distr;
		std::vector<float> OCC_distr;
		
				// simarity_matrix
		boost::numeric::ublas::matrix<double> _dis;

	protected:
		// the collection of points we are working on
		Points& _ps;
		
		// the scale of _ps
		unsigned int NumOfNodes;

		// mapping point_id -> clusterId
		std::vector<ClusterId> _pointId_to_clusterId;

		// the collection of clusters
		std::vector<Cluster> _clusters;








		friend	
			std::ostream& operator << (std::ostream& o, const Clusters& c);
	};

}
