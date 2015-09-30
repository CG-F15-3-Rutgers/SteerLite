//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function drawCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================

	// Robustness: make sure there is at least two control point: start and end points

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	
	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	std::cout<<"THIS METHOD HAS BEEN CALLED" << "\n";
	
	std::vector<CurvePoint> cPoints = getControPoints(); //vector to hold control points
	std::vector<CurvePoint> sortedControlPoints = getControPoints(); //vector to hold the sorted points
	
	int numP = controlPoints.size(); //number of control points
	int pos = 0; //position in vector of lowest point
	int count = cPoints.size(); //size of Cpoints
	
	
	for (int j = 0; j<count;j++){ //print all of the unsorted points
		std::cout<<"Control Point at " << j << " is: " << "Time: " << cPoints[j].time << " Position: " << cPoints[j].position <<"\n";
	}
	std::cout<<"\n";
	
	while(count>0){ //run until we have erased all of the points
		Util::CurvePoint lowestPoint = cPoints[0]; //set the first point to the lowest
		
		//std::cout<<"breakpoint1";
		for(int i = 0; i<count;i++){ //go through the entire array and save the lowest time
			//std::cout<<"breakpoint2";
			if(cPoints[i].time<lowestPoint.time){//new lowest point
				pos = i;
				lowestPoint=cPoints[i];
				//std::cout<<"The lowest point is: " << lowestPoint.time << " , " << lowestPoint.position;
			}else{
				//not lowest
			}
		}
		sortedControlPoints[numP - count]=lowestPoint; //populate sortedControlPoints with the lowest time
		//std::cout<<"breakpoint2.5";
		//std::cout<<"pos is: " << pos;
		cPoints.erase(cPoints.begin()+pos);	//erase 1 entry from cPoints
		pos=0;
		//std::cout<<"count is: "<<count;
		//std::cout<<"breakpoint3";
		count--;
	}
	//std::cout<<"breakpoint4"<<"\n";
	
	for (int k = 0; k<sortedControlPoints.size();k++){ //print all of the sorted points
		std::cout<<"Sorted Control Point at : " << k << " is: "<< "Time: " << sortedControlPoints[k].time << " Position: " << sortedControlPoints[k].position <<"\n";
	}
	
	return;
}


// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function checkRobust is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================


	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function findTimeInterval is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================


	return true;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function useHermiteCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================


	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Hermite curve

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function useCatmullCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================


	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Catmull-Rom curve
	
	// Return result
	return newPosition;
}