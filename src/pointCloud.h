#pragma once
#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H
#include "basicStructure.h"
#include <vector>
#include <string>
class PointCloud {
public:
	std::vector<Point> points;
	std::vector<double> normal;
	BoundingBox bb;

	PointCloud(std::string filePath);
	PointCloud(std::string filePath, std::vector<double> &normal);
	PointCloud(std::vector<Point> p);
};

#endif
