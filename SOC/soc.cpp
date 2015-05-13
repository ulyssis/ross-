//the criterion of ai can be used here!
//#include  "stdafx.h"
#include "soc.h"
#include <boost/foreach.hpp>

/*
namespace Clustering
{
	void DBSCAN::run_cluster() //claim member function run_clusters() out of the class of DBSCAN, so '::'is used here.
	{
		Cluster Neighbor_same_ai; //in case of neighboring node has the same value of a_i, use this vector to deal with it.
		ClusterId cid = 1;
		// foreach pid
		for (PointId pid = 0; pid < _ps.size(); pid++)
		{
			// not already visited
			if (!_visited[pid])
			    {
				_visited[pid] = true;
//此处执行一次pid=0,调用findNeighbors
				// get the neighbors
				Neighbors ne = findNeighbors(pid, _eps);


				unsigned int x=0,y=0,z=0;

				//If value of ai on this node is smaller than its neibors......
				for (unsigned int i = 0; i < ne.size(); i++)
				  {
/////////////////////  Bug  ////
				    x=x+((ai_vector[i]>ai_vector[pid])||(ai_vector[i]=ai_vector[pid]));
				    y=y+(ai_vector[i]>ai_vector[pid]);
				  }
				if(!x)
				  {
				  //pid is a basin node!
				  std::cout << "Point i=" << pid << " is a basin node " << std::endl;// = true;
				  // Add pid to current cluster
				  Cluster c;              // a new cluster
				  c.push_back(pid);       //this new cluster only contain node pid, which is the cluster head of this cluster.
				  for(unsigned int m=0; m < ne.size(); m++)
				      {
				  c.push_back(m);
				      }
				  _pointId_to_clusterId[pid]=cid;//collect nodes which are cluster heads
				  _clusters.push_back(c); //put the newly found cluster into _clusters.
				  cid++;
				  }
				else if(!y)                //x!=0 and y=0, there are neighboring nodes having identical ai with pid.
				  {for (unsigned int i = 0; i < ne.size(); i++)
				      {
				       if (ai_vector[pid]=ai_vector[i])
				      {
				       Neighbor_same_ai.push_back(i);
				      }
				      }
				      for(unsigned int j=0;j<Neighbor_same_ai.size();j++)
				      {
				      z=z+(bi_vector[pid]<bi_vector[j]);
				      }
				      if(!z)
				      {
				       //pid is a basin node!
				       std::cout << "Point i=" << pid << " is a basin node " << std::endl;// = true;
				       /// Add p to current cluster
	                                  Cluster c;              // a new cluster
	                                  c.push_back(pid);       //this new cluster only contain node pid, which is the cluster head of this cluster.
	                                  for(unsigned int m=0; m < ne.size(); m++)
	                                     {
	                                    c.push_back(m);
	                                     }
	                                  _pointId_to_clusterId[pid]=cid;//collect nodes which are cluster heads
	                                  _clusters.push_back(c); //put the newly found cluster into _clusters.
	                                  cid++;
				       }
				    }
				else{
				  std::cout<<"Point pid is not a basin node"<<std::endl;
				    }
				/*
				 if (ne.size() < _minPts)
				{
					_noise[pid] = true;
				}
				*/

/*
					// go to neighbors of P
					for (unsigned int i = 0; i < ne.size(); i++)
					{
						PointId nPid = ne[i];

						// not already visited
						if (!_visited[nPid])
						{
							_visited[nPid] = true;

							// go to neighbors
							Neighbors ne1 = findNeighbors(nPid, _eps);

							// enough support
							if (ne1.size() >= _minPts)
							{
								std::cout << "\t Expanding to pid=" << nPid << std::endl;	
								// join
								BOOST_FOREACH(Neighbors::value_type n1, ne1)
								{
									// join neighbord
									ne.push_back(n1); 
									std::cerr << "\tPushback pid=" << n1 << std::endl;
								}
								std::cout << std::endl;
							}
						}

						// not already assigned to a cluster
						if (!_pointId_to_clusterId[nPid])
						{
							std::cout << "\tadding pid=" << nPid << std::endl;
							c.push_back(nPid);
							_pointId_to_clusterId[nPid]=cid;
						}
					}


				}

			} // if (!visited

		} // for
	}
*/

