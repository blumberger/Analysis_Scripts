#ifndef CLUSTER_HEADER
#define CLUSTER_HEADER

#include <vector>
#include <iostream>
#include <set>
#include <cmath>

//		// The category of each point
//		std::vector<int> point_category(npts, -1);

namespace clust {

	/*
	 * Will cluster using the DBSCAN algorithm.
	 *
	 * Inputs:
	 *	pos (vec<vec<double>>) => Shape (natom, 3). The position coordinates
	 *	epislon (double)	   => The epislon parameter (how close points should be)
	 */
	std::vector<std::set<int>> naive(std::vector<std::vector<double>> const &pos, const double epsilon) {

		const int npts = pos.size();
		const double ep2 = epsilon*epsilon;
		std::vector<std::set<int>> all_cluster_inds;
		// atoms that haven't been visited by DBSCAN yet.
		std::set<int> visited_points;
		std::set<int> unvisited_points;
		for (auto i=0; i<npts; i++) unvisited_points.insert(i);

		auto pti = *unvisited_points.begin();
		std::set<int> cluster_buffer;
		while (unvisited_points.size() > 0) {
			auto &pos1 = pos[pti];
	
			cluster_buffer.insert(pti);
			visited_points.insert(pti);
			for (auto jpt: unvisited_points) {
				// Ignore categorised points again
				auto &pos2 = pos[jpt];
				double dist2 = std::pow(pos1[0]-pos2[0], 2) + 
							  std::pow(pos1[1]-pos2[1], 2) +
							  std::pow(pos1[2]-pos2[2], 2);
				
				if (dist2 < ep2) cluster_buffer.insert(jpt);
			}
			for (auto tmp: cluster_buffer) unvisited_points.erase(tmp);
	
			// A cluster of 1
			if (cluster_buffer.size() == 1) {
				all_cluster_inds.push_back(cluster_buffer);
				cluster_buffer = {};
				pti = *unvisited_points.begin();
			}
			// Decide how we choose the next point
			else {
				bool noPtsRemaining = true;
				for (auto i: cluster_buffer) {
					if (visited_points.find(i) == visited_points.end()) {
						noPtsRemaining = false;
						pti = i;
						break;
					}
				}
				if (noPtsRemaining == true) {
					all_cluster_inds.push_back(cluster_buffer);
				cluster_buffer = {};
					pti = *unvisited_points.begin();
				}
			}
		}


		return all_cluster_inds;
	}
}

#endif
