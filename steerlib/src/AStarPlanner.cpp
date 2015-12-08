//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//

#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
#include <functional>
#include <queue>
#include <math.h>
#include "planning/AStarPlanner.h"

#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MANHATTAN false

namespace SteerLib
{
	AStarPlanner::AStarPlanner() {}

	AStarPlanner::~AStarPlanner() {}

	bool AStarPlanner::canBeTraversed(int id, int ob)
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x, z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		x_range_min = (x < ob) ? 0 : (x - ob);
		x_range_max = MIN(x+ob, gSpatialDatabase->getNumCellsX() - 1);

		z_range_min = (z < ob) ? 0 : (z - ob);
		z_range_max = MIN(z+ob, gSpatialDatabase->getNumCellsZ() - 1);

		for (int i = x_range_min; i <= x_range_max; i += GRID_STEP)
		{
			for (int j = z_range_min; j <= z_range_max; j += GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords(i, j);
				traversal_cost += gSpatialDatabase->getTraversalCost(index);
			}
		}

		if (traversal_cost > COLLISION_COST)
			return false;
		return true;
	}

	
	
	float heuristic_function(const Util::Point& p1, const Util::Point& p2) {
		
		if (MANHATTAN) {
			return (abs(p2.x - p1.x) + abs(p2.y - p1.y) + abs(p2.z - p1.z));
		}
		else {
			return Util::distanceBetween(p1, p2);
		}
		

	}

	int search(const std::vector<AStarPlannerNode>& set, const Util::Point& point, SteerLib::GridDatabase2D* _gSpatialDatabase)
     {
       
		int point_index = _gSpatialDatabase->getCellIndexFromLocation(point);
        for (size_t i = 0; i < set.size(); ++i) {
            int curr_index = _gSpatialDatabase->getCellIndexFromLocation(set[i].point);
            if (curr_index == point_index)
                 return i;
        }
         return -1;
     }


	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D* _gSpatialDatabase, bool append_to_path,bool diagonal, int ob )
	{
		gSpatialDatabase = _gSpatialDatabase;

		bool path = false;
		float weight = 1;
		std::vector<AStarPlannerNode> closed;
		std::vector<AStarPlannerNode> open;
		
		open.push_back(AStarPlannerNode(start, weight*heuristic_function(start, goal), 0, -1));

		while (open.size() > 0) {
			int cDex = 0;

			double min_f = open[0].f;

			for (int i = 1; i < open.size(); ++i) {
				if (open[i].f < min_f) {
					min_f = open[i].f;
					cDex = i;
				}
			}


			int goal_grid_index = _gSpatialDatabase->getCellIndexFromLocation(goal);
            int curr_grid_index = _gSpatialDatabase->getCellIndexFromLocation(open[cDex].point);

            if (goal_grid_index == curr_grid_index) {
				path = true;
				break;
			}

			
			closed.push_back(open[cDex]);
			open.erase(open.begin() + cDex);

			AStarPlannerNode& curr = closed.back();

			
			std::vector<Util::Point> prevNode;
			
				//int minx = MAX(gSpatialDatabase->getOriginX(), curr.point.x - 1);
				//int maxx = MIN(curr.point.x + 1, gSpatialDatabase->getNumCellsX() + gSpatialDatabase->getOriginX());

				//int minz = MAX(gSpatialDatabase->getOriginZ(), curr.point.z - 1);
				//int maxz = MIN(curr.point.z + 1, gSpatialDatabase->getNumCellsZ() + gSpatialDatabase->getOriginZ());
				int startIndex = gSpatialDatabase->getCellIndexFromLocation(curr.point);
				unsigned int x, z;
				gSpatialDatabase->getGridCoordinatesFromIndex(startIndex, x, z);

				int minx = (x == 0) ? 0 : (x - 1);
				int maxx = MIN(x + 1, gSpatialDatabase->getNumCellsX() - 1);

				int minz = (z == 0) ? 0 : (z - 1);
				int maxz = MIN(z + 1, gSpatialDatabase->getNumCellsZ() - 1);


				for (int i = minx; i <= maxx; i++) {
					for (int j = minz; j <= maxz; j++) {
						/*if (!(i == curr.point.x && j == curr.point.z)) {
							int index = gSpatialDatabase->getCellIndexFromLocation(i, j);
							if (canBeTraversed(index))
								prevNode.push_back(Util::Point(i, 0, j));
						}*/
						
						if (!(i == x && j == z)) {
							int index = gSpatialDatabase->getCellIndexFromGridCoords(i, j);
							if (canBeTraversed(index, ob)) {
								Util::Point p;
								gSpatialDatabase->getLocationFromIndex(index, p);
								prevNode.push_back(p);
							}
						}
					}
				}
			

			//add
			for (int i = 0; i < prevNode.size(); ++i) {
				int closedDex = search(closed,prevNode[i], _gSpatialDatabase);
				if (closedDex != -1) {
					continue;
				}

				double g = curr.g + 1;
										  
				double f = g + weight * heuristic_function(prevNode[i], goal);


				int o = search(open, prevNode[i],_gSpatialDatabase);
			

				if (o == -1) {
					
					open.push_back(AStarPlannerNode(prevNode[i], f, g, closed.size() - 1));
				}
				else if (g < open[o].g) {
					
					open[o].f = f;
					open[o].g = g;
					open[o].parent_index = closed.size() - 1;
				}
			}
		}

		if (path==true) {
			if (!append_to_path) {
				agent_path.clear();
			}
			int gDex = search(open,goal,_gSpatialDatabase);
			assert(gDex != -1);

		
			std::vector<Util::Point> temp;
			const AStarPlannerNode* current = &open[gDex];
			while (current->point != start) {
				temp.push_back(current->point);
				current = &closed[current->parent_index];
			}
			temp.push_back(start);

			
			for (std::vector<Util::Point>::reverse_iterator iter = temp.rbegin(); iter != temp.rend(); iter++) {
				agent_path.push_back(*iter);
			}
		}

		return path;
	}
}