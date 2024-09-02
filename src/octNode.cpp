#include "octNode.h"
OctNode::OctNode(int i, int d)
{
	nid = i;
	depth = d;
}
OctNode::OctNode(int i, int d, Point c)
{
	nid = i;
	depth = d;
	center = c;
}
OctNode::OctNode(int i, int d, int f, Point c)
{
	nid = i;
	depth = d;
	fa = f;
	center = c;
}
OctNode::OctNode(int i, int d, Point c, const std::vector<Point> &points, const std::vector<int> &pointId)
{
	nid = i;
	depth = d;
	center = c;

	for (int i = 0; i < points.size(); i++) {
		inPoints.push_back(points[i]);
		inPointId.push_back(pointId[i]);
	}
}
void OctNode::set_fa(int f)
{
	fa = f;
}
void OctNode::set_lnid(int i)
{
	lnid = i;
}
void OctNode::add_point(Point p, int pid)
{
	inPoints.push_back(p);
	inPointId.push_back(pid);
}
void OctNode::add_around(int ni)
{
	around.push_back(ni);
}
void OctNode::creat_childNode(std::vector<OctNode *> &nodes)
{
	int beginId = nodes.size();
	double l = DEPLEN[depth - 1];
	isleaf = false;
	for (int i = 0; i < 8; i++) {
		int dx = 2 * (i % 2) - 1;
		int dy = 2 * ((i >> 1) % 2) - 1;
		int dz = 2 * ((i >> 2) % 2) - 1;
		Point c = Point(center.x + dx * l / 4, center.y + dy * l / 4, center.z + dz * l / 4);
		nodes.push_back(new OctNode(beginId + i, depth + 1, nid, c));
		child[i] = beginId + i;
	}
	for (int i = 0; i < inPoints.size(); i++) {
		Point p = inPoints[i];
		int id = 0;
		id += p.x - center.x > 0 ? 1 : 0;
		id += p.y - center.y > 0 ? 2 : 0;
		id += p.z - center.z > 0 ? 4 : 0;
		nodes[beginId + id]->add_point(p, inPointId[i]);
	}
	inPoints.clear();
	inPoints.shrink_to_fit();
	inPointId.clear();
	inPointId.shrink_to_fit();
}

bool OctNode::isNeedDivide(int maxDepth)
{
	if (depth < APTDEPTH) {
		return true;
	}
	if (depth < maxDepth && inPoints.size() > NPN) {
		return true;
	}
	return false;
}




