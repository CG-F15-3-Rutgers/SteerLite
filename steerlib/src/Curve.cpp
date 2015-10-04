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
 		//std::cerr << "Error, checkRobust failed in drawCurve!";
 	}
 	
 	/*float timeInterval, startT, endT, dT, time;
 	Point endPoint;
 	Point startPoint;
 	int curr;
 	
 	for(int i = 1; i<controlPoints.size(); i++){	//for each control point in controlPoints
 		
 	
 		startPoint = controlPoints[i-1].position;	//get the start point, start time, and end time
 		startT = controlPoints[i-1].time;
 		endT = controlPoints[i].time;
 		timeInterval = endT-startT;					//calculate the delta which is just the interval/window
 		dT = timeInterval/(float)window;
 		curr = i;
 		
 		for(time =startT; time <= endT; time=time+dT){ //go through the time adding the delta time constant to draw either the hermite or the catmull curve
 			
 			if(type==hermiteCurve){
 				endPoint = useHermiteCurve(curr,time);
 			}else{
 				endPoint = useCatmullCurve(curr,time);
 			}
 			
 			DrawLib::drawLine(startPoint,endPoint,curveColor,curveThickness); //draw function
 		}
 
 
 	}*/
	
	Point startPoint = controlPoints[0].position;
 	Point endPoint;
	
	float time = 0;
 	float endTime = controlPoints[controlPoints.size() - 1].time;
 	float startTime = controlPoints[0].time;
	
	
	float timeWin = 0.5;
 	while (time <= endTime) {			//while time is less than or equal to the endtime
 		if (calculatePoint(endPoint, time)) { //if calculate point spits back true we can keep drawing the line
 			DrawLib::drawLine(startPoint, endPoint, curveColor, curveThickness);
 			startPoint = endPoint;
 		} else {			//if it spits out false we have problems -.-
 			//std::cerr<<"Failed to find next Point at time " << time << std::endl;
 		}
 		time = time + timeWin; //go to next time window
	}
	
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
 	if (!findTimeInterval(nextPoint, time)){
			
		//std::cout<<"failed in time interval"<<"\n";
 		return false;
	}
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
	//outputPoint = useHermiteCurve(nextPoint, time);
	//outputPoint = useCatmullCurve(nextPoint, time);
 
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
	
	std::vector<CurvePoint> cPoints = getControPoints();	//vector
 	std::vector<CurvePoint>::iterator iter;					//vector iterator
 	iter = cPoints.begin();
 	for(int a=0;a<cPoints.size();a++){
 		
 		if(time<cPoints[a].time){		//if time is less than the a'th spot in cPoints, set the distance, and return true
 			nextPoint = std::distance(cPoints.begin(),iter);
 			return true;
 			}
 		iter++;
 	}
	//std::cout<<"error in findTimeInterval"<<"\n"<<"\n";
 	return false;		//if the time is outside of the time interval we return false
	
	
 }
 
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
	
	
	float t = normalTime;
	float t2 = normalTime*normalTime;
	float t3 = normalTime*normalTime*normalTime;
	
	float x0,x1,x2,x3;
	float y0,y1,y2,y3;
	float z0,z1,z2,z3;
 	float x, y, z;
		
 	x0 = (2*t3 - 3*t2 + 1)*prevP.x;		//hermite math for x
	x1 = (-2*t3+3*t2)*nextP.x;
	x2 = (t3 -2*t2+t)*curveP[nextPoint-1].tangent.x;
	x3 = (t3 - t2)*curveP[nextPoint].tangent.x;
	
	x = x0 + x1 + x2 + x3;
	
	y0 = (2*t3 - 3*t2 + 1)*prevP.y;		//hermite math for y
	y1 = (-2*t3+3*t2)*nextP.y;
	y2 = (t3 -2*t2+t)*curveP[nextPoint-1].tangent.y;
	y3 = (t3 - t2)*curveP[nextPoint].tangent.y;
	
	y = y0 + y1 + y2 + y3;
		
	z0 = (2*t3 - 3*t2 + 1)*prevP.z;		//hermite math for x
	z1 = (-2*t3+3*t2)*nextP.z;
	z2 = (t3 -2*t2+t)*curveP[nextPoint-1].tangent.z;
	z3 = (t3 - t2)*curveP[nextPoint].tangent.z;
	
	z = z0 + z1 + z2 + z3;
 
 	newPosition = Point(x,y,z);
	
	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Hermite curve

	// Return result
	return newPosition;
	
 }
 
 // Implement Catmull-Rom curve
 Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
 {
 	Point newPosition;
 	float normalTime, intervalTime;
 
 	
 	int size = controlPoints.size();
 	float nextTime = controlPoints[nextPoint].time;
 	float prevTime = controlPoints[nextPoint-1].time;
 	intervalTime = nextTime-prevTime;						//interval time is just the later time - earlier time
 	normalTime = ((time-prevTime)/intervalTime);			//normal time is interval time normalized
 	
	if(nextPoint==0){
		newPosition = controlPoints[nextPoint].position;
		return newPosition;
	}
	
 	if(nextPoint<=1){
 		newPosition = useHermiteCurve(nextPoint,time);
 		return newPosition;
 	}
 	if(nextPoint>=size-1){
 		newPosition = useHermiteCurve(nextPoint,time);
 		return newPosition;
 	}
 	
 	Point point3 = controlPoints[nextPoint+1].position;
 	Point point2 = controlPoints[nextPoint].position;
 	Point point1 = controlPoints[nextPoint-1].position;
 	Point point0 = controlPoints[nextPoint-2].position;//get the nextPoint and currentPoint's time and positions
 	
 	float t = normalTime;
 	float t2 = normalTime*normalTime;
 	float t3 = normalTime*normalTime*normalTime;

 	
	newPosition.x = ((2 * point1.x) + (-point0.x + point2.x) * t + (2*point0.x - 5*point1.x + 4*point2.x - point3.x) * t2 + (-point0.x + 3*point1.x - 3*point2.x + point3.x) * t3) * 0.5f;
	newPosition.x = ((2 * point1.x) + (-point0.x + point2.x) * t + (2*point0.x - 5*point1.x + 4*point2.x - point3.x) * t2 + (-point0.x + 3*point1.x - 3*point2.x + point3.x) * t3) * 0.5f;	//catmull-rom math
 	newPosition.y = ((2 * point1.y) + (-point0.y + point2.y) * t + (2*point0.y - 5*point1.y + 4*point2.y - point3.y) * t2 + (-point0.y + 3*point1.y - 3*point2.y + point3.y) * t3) * 0.5f;
 	newPosition.z = ((2 * point1.z) + (-point0.z + point2.z) * t + (2*point0.z - 5*point1.z + 4*point2.z - point3.z) * t2 + (-point0.z + 3*point1.z - 3*point2.z + point3.z) * t3) * 0.5f;
 	
	
 	
 	// Calculate time interval, and normal time required for later curve calculations
 
 	// Calculate position at t = time on Catmull-Rom curve
 	
 	// Return result
 	return newPosition;
 }