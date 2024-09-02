#include "pointCloud.h"
#include <iostream>
#include <fstream>

PointCloud::PointCloud(std::string filePath)
{
	std::ifstream fin;
	fin.open(filePath, std::ios::in);
	if (!fin.is_open()) {
		std::cout << filePath << "文件读取错误" << std::endl;
		return;
	}
	points.clear();
	double x, y, z;
	double maxX , maxY, maxZ;
	double minX, minY, minZ;
	if (fin >> x && fin >> y && fin >> z) {
		maxX = x;
		minX = x;
		maxY = y;
		minY = y;
		maxZ = z;
		minZ = z;
		points.push_back(Point(x, y, z));
	}
	while (fin >> x && fin >> y && fin >> z) {
		points.push_back(Point(x, y, z));
		maxX = maxX < x ? x : maxX;
		minX = minX > x ? x : minX;
		maxY = maxY < y ? y : maxY;
		minY = minY > y ? y : minY;
		maxZ = maxZ < z ? z : maxZ;
		minZ = minZ > z ? z : minZ;
	}
	fin.close();
	bb.pxMax = maxX;
	bb.pxMin = minX;
	bb.pyMax = maxY;
	bb.pyMin = minY;
	bb.pzMax = maxZ;
	bb.pzMin = minZ;

	bb.scale = maxX - minX;
	double scale = maxY - minY;
	bb.scale = bb.scale < scale ? scale : bb.scale;
	scale = maxZ - minZ;
	bb.scale = bb.scale < scale ? scale : bb.scale;
	bb.scale *= BSCALE;
	//bb.blx = 0.5 * ((1 + BSCALE) * minX + (1 - BSCALE) * maxX);
	//bb.bly = 0.5 * ((1 + BSCALE) * minY + (1 - BSCALE) * maxY);
	//bb.blz = 0.5 * ((1 + BSCALE) * minZ + (1 - BSCALE) * maxZ);
	bb.blx = 0.5 * (minX + maxX - bb.scale);
	bb.bly = 0.5 * (minY + maxY - bb.scale);
	bb.blz = 0.5 * (minZ + maxZ - bb.scale);

}

PointCloud::PointCloud(std::string filePath, std::vector<double> &normal)
{
	std::ifstream fin;
	fin.open(filePath, std::ios::in);
	if (!fin.is_open()) {
		std::cout << filePath << "文件读取错误" << std::endl;
		return;
	}
	points.clear();
	double x, y, z, nx, ny, nz;
	double maxX , maxY, maxZ;
	double minX, minY, minZ;
	normal.clear();
	if (fin >> x && fin >> y && fin >> z && fin >> nx && fin >> ny && fin >> nz) {
		maxX = x;
		minX = x;
		maxY = y;
		minY = y;
		maxZ = z;
		minZ = z;
		points.push_back(Point(x, y, z));
		normal.push_back(nx);
		normal.push_back(ny);
		normal.push_back(nz);
	}
	while (fin >> x && fin >> y && fin >> z && fin >> nx && fin >> ny && fin >> nz) {
		points.push_back(Point(x, y, z));
		maxX = maxX < x ? x : maxX;
		minX = minX > x ? x : minX;
		maxY = maxY < y ? y : maxY;
		minY = minY > y ? y : minY;
		maxZ = maxZ < z ? z : maxZ;
		minZ = minZ > z ? z : minZ;
		normal.push_back(nx);
		normal.push_back(ny);
		normal.push_back(nz);
	}
	fin.close();
	bb.pxMax = maxX;
	bb.pxMin = minX;
	bb.pyMax = maxY;
	bb.pyMin = minY;
	bb.pzMax = maxZ;
	bb.pzMin = minZ;

	bb.scale = maxX - minX;
	double scale = maxY - minY;
	bb.scale = bb.scale < scale ? scale : bb.scale;
	scale = maxZ - minZ;
	bb.scale = bb.scale < scale ? scale : bb.scale;
	bb.scale *= BSCALE;
	/*bb.blx = 0.5 * ((1 + BSCALE) * minX + (1 - BSCALE) * maxX);
	bb.bly = 0.5 * ((1 + BSCALE) * minY + (1 - BSCALE) * maxY);
	bb.blz = 0.5 * ((1 + BSCALE) * minZ + (1 - BSCALE) * maxZ);*/
	bb.blx = 0.5 * (minX + maxX - bb.scale);
	bb.bly = 0.5 * (minY + maxY - bb.scale);
	bb.blz = 0.5 * (minZ + maxZ - bb.scale);
}

PointCloud::PointCloud(std::vector<Point> p)
{
	points.clear();
	double maxX, maxY, maxZ;
	double minX, minY, minZ;
	maxX = p[0].x;
	minX = p[0].x;
	maxY = p[0].y;
	minY = p[0].y;
	maxZ = p[0].z;
	minZ = p[0].z;
	for (int i = 0; i < p.size(); i++) {
		points.push_back(p[i]);
		maxX = maxX < p[i].x ? p[i].x : maxX;
		minX = minX > p[i].x ? p[i].x : minX;
		maxY = maxY < p[i].y ? p[i].y : maxY;
		minY = minY > p[i].y ? p[i].y : minY;
		maxZ = maxZ < p[i].z ? p[i].z : maxZ;
		minZ = minZ > p[i].z ? p[i].z : minZ;
	}
	bb.pxMax = maxX;
	bb.pxMin = minX;
	bb.pyMax = maxY;
	bb.pyMin = minY;
	bb.pzMax = maxZ;
	bb.pzMin = minZ;

	bb.scale = maxX - minX;
	double scale = maxY - minY;
	bb.scale = bb.scale < scale ? scale : bb.scale;
	scale = maxZ - minZ;
	bb.scale = bb.scale < scale ? scale : bb.scale;
	bb.scale *= BSCALE;
	//bb.blx = 0.5 * ((1 + BSCALE) * minX + (1 - BSCALE) * maxX);
	//bb.bly = 0.5 * ((1 + BSCALE) * minY + (1 - BSCALE) * maxY);
	//bb.blz = 0.5 * ((1 + BSCALE) * minZ + (1 - BSCALE) * maxZ);
	bb.blx = 0.5 * (minX + maxX - bb.scale);
	bb.bly = 0.5 * (minY + maxY - bb.scale);
	bb.blz = 0.5 * (minZ + maxZ - bb.scale);

}

