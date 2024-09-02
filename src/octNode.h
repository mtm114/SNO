#pragma once
#ifndef OCTNODE_H
#define OCTNODE_H
#include <vector>
#include "basicStructure.h"

class OctNode
{
public:
	int nid;
	int lnid;
	bool isleaf = true;
	std::vector<Point> inPoints;
	std::vector<int> inPointId;
	std::vector<int> around;
	int fa = -1;
	int child[8];
	int depth;
	Point center;

	OctNode(int i, int d);
	OctNode(int i, int d, Point c);
	OctNode(int i, int d, int f, Point c);
	OctNode(int i, int d, Point c, const std::vector<Point> &points, const std::vector<int> &pointId);
	void set_lnid(int i);
	void set_fa(int f);
	void add_point(Point p, int pid);
	void add_around(int ni);
	void creat_childNode(std::vector<OctNode *> &nodes);
	bool isNeedDivide(int maxDepth);
	
private:

};

#endif // !OCTNODE_H
