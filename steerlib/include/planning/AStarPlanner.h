//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//


#ifndef __STEERLIB_A_STAR_PLANNER_H__
#define __STEERLIB_A_STAR_PLANNER_H__


#include <vector>
#include <stack>
#include <set>
#include <map>
#include "SteerLib.h"

namespace SteerLib
{

	class STEERLIB_API AStarPlannerNode {
	public:
		double f;
		double g;
		size_t parent_index;
		Util::Point point;

		AStarPlannerNode(Util::Point _point, double _f, double _g, size_t _parent_index)
		{
			f = _f;
			point = _point;
			g = _g;
			parent_index = _parent_index;
		}

		bool operator<(const AStarPlannerNode& other) const
		{
			return this->f < other.f;
		}

		bool operator>(const AStarPlannerNode& other) const
		{
			return this->f > other.f;
		}

		bool operator==(const AStarPlannerNode& other) const
		{
			return ((this->point.x == other.point.x) &&
				(this->point.z == other.point.z));
		}
	};


	class STEERLIB_API AStarPlanner {
	public:
		AStarPlanner();
		~AStarPlanner();


		bool canBeTraversed(int id, int ob);


		Util::Point getPointFromGridIndex(int id);

		bool computePath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D* _gSpatialDatabase, bool append_to_path = true, bool diagonal=true, int ob = 0);

	private:
		SteerLib::GridDatabase2D* gSpatialDatabase;

		std::vector<Util::Point> getSuccessors(const Util::Point& p);
	};


}


#endif