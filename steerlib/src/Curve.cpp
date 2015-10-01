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

	if(checkRobust()==false){
		std::cerr << "Error, checkRobust failed in drawCurve!";
		return;
	}
	
	float time = 0;
	Point startP = controlPoints[0].position;
	Point endP;
	float tWindow = 0.5;
	float end = controlPoints[controlPoints.size()-1].time;
	while (tWindow <= end){
		if(calculatePoint(endP,time)){
			DrawLib::drawLine(startP,endP,curveColor,curveThickness);
			startP =endP;
		}else{
			std::cerr<<"Error can't find the point";
		}
		time = time + tWindow;
	}

	// Robustness: make sure there is at least two control point: start and end points

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	
	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
bool compareFunction(CurvePoint cp1, CurvePoint cp2) {
 	return (cp1.time < cp2.time);
 }
 
 // Sort controlPoints vector in ascending order: min-first
 void Curve::sortControlPoints()
 {
 
 	for(std::vector<CurvePoint>::iterator it = controlPoints.begin(); it < controlPoints.end(); ++it) {
 		std::cout << it->time << std::endl;
 	}
 
 	std::sort(controlPoints.begin(), controlPoints.end(), compareFunction);
 	for(std::vector<CurvePoint>::iterator it = controlPoints.begin(); it < controlPoints.end(); ++it) {
 		std::cout << it->time << std::endl;
 	}
 }


// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust()){
		std::cout<<"error with check ROBUST"<<"\n"<<"\n";
		return false;
		
	}
	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time)){
		std::cout<<"error with findTimeInterval"<<"\n"<<"\n";
		return false;
	}

	// Calculate position at t = time on curve
	//if (type == hermiteCurve)
	//{
	outputPoint = useHermiteCurve(nextPoint, time);
	//}
	//else if (type == catmullCurve)
	//{
	//	outputPoint = useCatmullCurve(nextPoint, time);
	//}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{	
	std::vector<CurvePoint> rob_curves = getControPoints(); // vector
	
	if (rob_curves.size() < 2) {	//if there are not two points... bad!
		return false;
	}
	else {
		return true;				//if there are... good!
	}
	return true;					//shouldn't happen
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{

	/*if(!((controlPoints[0].time<=time&&controlPoints.back().time>=time))){
		std::cout<<"error with findTimeInterval if statement";
		return false;
	}*/
	std::vector<CurvePoint>::iterator it;
	for(it=controlPoints.begin();it!=controlPoints.end();++it){
		
		if(it->time > time){
			nextPoint = it - controlPoints.begin();
			return true;
		}
	}
	return false;
}

<<<<<<< HEAD
float h1 (float t) {
 	return (2 * t * t * t - 3 * t * t + 1);
 }
 
 float h2 (float t) {
 	return (-2 * t * t * t + 3 * t * t);
 }
 
 float h3 (float t) {
 	return (t * t * t - 2 * t * t + t);
 }
 
 float h4 (float t) {
 	return (t * t * t - t * t);
 } 
 
float dist (Point p0, Point p1){
	float diffx = pow((p0.x - p1.x), 2);
	float diffy = pow((p0.y - p1.y), 2);
	float diffz = pow((p0.z - p1.z), 2);
=======
// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	std::vector<CurvePoint> curveP = getControPoints();		//vector of controlPoints
	Point nextP = curveP[nextPoint].position;
	Point prevP = curveP[nextPoint-1].position;				//get the nextPoint and currentPoint's time and positions
	float nextTime = curveP[nextPoint].time;
	float prevTime = curveP[nextPoint-1].time;
	
	intervalTime = nextTime-prevTime;						//interval time is just the later time - earlier time
	normalTime = ((time-prevTime)/intervalTime);			//normal time is interval time normalized
	
	
	
	
	
	
	
	
	
	
	//==========================================================================================
	//			WE
	//					STILL
	//								NEED
	//											TO
	//													IMPLEMENT
	//																	THIS
	//																			STUFF
	//===========================================================================================
	
	
	
	
	
	
	
	
	
	
	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Hermite curve
>>>>>>> parent of b9e553d... added hermite curve

	return sqrt(diffx + diffy + diffz);
}

 // Implement Hermite curve
 Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
 {
 	Point newPosition;
 	float normalTime, intervalTime;
 
 	//=========================================================================
 
 	if (nextPoint == 0) {
 		newPosition = controlPoints[nextPoint].position;
 		return newPosition;
 	}
 	// Calculate time interval, and normal time required for later curve calculations
 	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
 	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;
 	// Calculate position at t = time on Hermite curve
 	Point p0 = controlPoints[nextPoint - 1].position;
 	Point p1 = controlPoints[nextPoint].position;
 	Vector v0 = controlPoints[nextPoint - 1].tangent;
 	Vector v1 = controlPoints[nextPoint].tangent;
 	float _h1 = h1(normalTime);
 	float _h2 = h2(normalTime);
 	float _h3 = h3(normalTime);
 	float _h4 = h4(normalTime);
	
	float _dist = dist(p0, p1);

	newPosition.x = p0.x * _h1 + p1.x * _h2 + v0.x * _h3 * _dist + v1.x * _h4 * _dist;
	newPosition.y = p0.y * _h1 + p1.y * _h2 + v0.y * _h3 * _dist + v1.y * _h4 * _dist;
	newPosition.z = p0.z * _h1 + p1.z * _h2 + v0.z * _h3 * _dist + v1.z * _h4 * _dist;
 	// Return result
 	return newPosition;
 }

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	std::vector<CurvePoint> curveP = getControPoints();		//vector of controlPoints
	int size = curveP.size();
	float nextTime = curveP[nextPoint].time;
	float prevTime = curveP[nextPoint-1].time;
	
	intervalTime = nextTime-prevTime;						//interval time is just the later time - earlier time
	normalTime = ((time-prevTime)/intervalTime);			//normal time is interval time normalized
	
	if(nextPoint<=1){
		newPosition = useHermiteCurve(nextPoint,time);
		return newPosition;
	}
	if(nextPoint>=size-1){
		newPosition = useHermiteCurve(nextPoint,time);
		return newPosition;
	}
	
	Point point3 = curveP[nextPoint+1].position;
	Point point2 = curveP[nextPoint].position;
	Point point1 = curveP[nextPoint-1].position;
	Point point0 = curveP[nextPoint-2].position;//get the nextPoint and currentPoint's time and positions
	
	float t = normalTime;
	float t2 = normalTime*normalTime;
	float t3 = normalTime*normalTime*normalTime;
	
	/*newPosition.position.x = ((-t3 + 2*t2-t)*(point0.position.x) + (3*t3-5*t2+2)*(point1.position.x) + (-3*t3+4*t2+t)* (point2.position.x) + (t3-t2)*(point4.position.x))/2;
	newPosition.position.y = ((-t3 + 2*t2-t)*(point0.position.y) + (3*t3-5*t2+2)*(point1.position.y) + (-3*t3+4*t2+t)* (point2.position.y) + (t3-t2)*(point4.position.y))/2;
	newPosition.position.z = ((-t3 + 2*t2-t)*(point0.position.z) + (3*t3-5*t2+2)*(point1.position.z) + (-3*t3+4*t2+t)* (point2.position.z) + (t3-t2)*(point4.position.z))/2;
	*/
	
	newPosition.x = ((2 * point1.x) + (-point0.x + point2.x) * t + (2*point0.x - 5*point1.x + 4*point2.x - point3.x) * t2 + (-point0.x + 3*point1.x - 3*point2.x + point3.x) * t3) * 0.5f;
	newPosition.y = ((2 * point1.y) + (-point0.y + point2.y) * t + (2*point0.y - 5*point1.y + 4*point2.y - point3.y) * t2 + (-point0.y + 3*point1.y - 3*point2.y + point3.y) * t3) * 0.5f;
	newPosition.z = ((2 * point1.z) + (-point0.z + point2.z) * t + (2*point0.z - 5*point1.z + 4*point2.z - point3.z) * t2 + (-point0.z + 3*point1.z - 3*point2.z + point3.z) * t3) * 0.5f;
	
	return newPosition;
	
	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Catmull-Rom curve
	
	// Return result
	return newPosition;
}