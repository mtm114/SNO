#include "octTree.h"
#include <queue>
#include <iostream>
#include <unordered_set>
#include <stack>
#include <ctime>
#include <fstream>
#include <Eigen/SVD>  
#include <Eigen/Eigen>
typedef std::pair<int, double> pad;
OctTree::OctTree(PointCloud *p)
{
	pc = p;
	unitPointCloud();
	std::cout << "build Oct" << std::endl;
	buildTree();
	get_outNode();
	find_Ngbr();
	find_around();
	std::cout << "calc_B" << std::endl;
	calc_B();
	std::cout << "calc_A" << std::endl;
	calc_A();
	find_near_points();
	calc_P();
}


void OctTree::unitPointCloud()
{
	Point origin = Point(pc->bb.blx, pc->bb.bly, pc->bb.blz);
	double scale = pc->bb.scale;
	pxMax = (pc->bb.pxMax - origin.x) / scale;
	pxMin = (pc->bb.pxMin - origin.x) / scale;
	pyMax = (pc->bb.pyMax - origin.y) / scale;
	pyMin = (pc->bb.pyMin - origin.y) / scale;
	pzMax = (pc->bb.pzMax - origin.z) / scale;
	pzMin = (pc->bb.pzMin - origin.z) / scale;
	points.clear();
	for (int i = 0; i < pc->points.size(); i++) {
		points.push_back((pc->points[i] - origin) * (1 / scale));
	}
}
void OctTree::buildTree()
{
	unsigned nv = points.size();
	ANNpointArray dataPts = annAllocPts(nv, 3);
	for (int i = 0; i < nv; i++) {
		dataPts[i][0] = points[i].x;
		dataPts[i][1] = points[i].y;
		dataPts[i][2] = points[i].z;
	}

	if (kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nv, 3);
	std::queue<int> que;
	std::vector<int> pid(points.size());
	for (int i = 0; i < points.size(); i++) {
		pid[i] = i;
	}

	calc_avgdis();

	pdep.clear();
	pdep.resize(points.size(), 0);
	for (int i = 0; i < points.size(); i++) {
		int logdis = round(-log2(avgdis[i])) + 1;
		if (logdis < MAXDEPTH) {
			pdep[i] = logdis;
		}
		
	}

	nodeArr.clear();
	nodeArr.push_back(new OctNode(0, 1, Point(0.5, 0.5, 0.5), points, pid));
	que.push(0);
	while (!que.empty()) {
		int nodeid = que.front();
		que.pop();
		int bid = nodeArr.size();
		nodeArr[nodeid]->creat_childNode(nodeArr);
		for (int i = 0; i < 8; i++) {
			if (nodeArr[bid + i]->isNeedDivide(MAXDEPTH)) {
				que.push(bid + i);
			}
		}
	}

	for (int i = 0; i < nodeArr.size(); i++) {
		OctNode *nodei = nodeArr[i];
		int dep = nodei->depth;
		if (!nodei->isleaf || dep >= MAXDEPTH) {
			continue;
		}
		if (nodei->inPointId.size() > 0) {
			bool needDiv = false;
			for (int j = 0; j < nodei->inPointId.size(); j++) {
				if (pdep[nodei->inPointId[j]] > dep) {
					needDiv = true;
					break;
				}
			}
			if (needDiv) {
				nodei->creat_childNode(nodeArr);
			}
		}
	}
	
	pointNodeId.clear();
	pointNodeId.resize(points.size());
	for (int i = 0; i < nodeArr.size(); i++) {
		for (int j = 0; j < nodeArr[i]->inPointId.size(); j++) {
			pointNodeId[nodeArr[i]->inPointId[j]] = i;
			pdep[nodeArr[i]->inPointId[j]] = nodeArr[i]->depth;
		}
	}

	

	for (int i = 0; i < nodeArr.size(); i++) {
		OctNode *nodei = nodeArr[i];
		if (!nodei->isleaf || nodei->inPointId.size() > 0 || nodei->depth >= MAXDEPTH) {
			continue;
		}
		Point o = nodei->center;
		ANNpoint tp = annAllocPt(3);
		tp[0] = o.x;
		tp[1] = o.y;
		tp[2] = o.z;
		ANNidxArray nnIdx = new ANNidx[8];
		ANNdistArray dists = new ANNdist[8];
		kdTree->annkSearch(tp, 8, nnIdx, dists);
		int pd = 0;
		for (int k = 0; k < 8; k++) {
			pd += nodeArr[pointNodeId[nnIdx[k]]]->depth;
		}
		pd = round(pd / 8.0);
		if (nodei->depth < pd && dists[0] < 3 * DEPLEN[nodei->depth] * DEPLEN[nodei->depth]) {
			nodei->creat_childNode(nodeArr);
		}
	}


	leafNodeArr.clear();
	int lid = 0;
	for (int i = 0; i < nodeArr.size(); i++) {
		if (nodeArr[i]->isleaf) {
			leafNodeArr.push_back(nodeArr[i]);
			nodeArr[i]->set_lnid(lid);
			lid++;

			depthNode[nodeArr[i]->depth - 1].push_back(nodeArr[i]->nid);
		}
	}


	for (int i = 0; i < leafNodeArr.size(); i++) {
		if (leafNodeArr[i]->inPointId.size() == 0) {
			nopNode.push_back(leafNodeArr[i]->nid);
		}
		//nopNode.push_back(leafNodeArr[i]->nid);
	}


}

void insertpad(std::vector<pad> &parr, pad pa, int kn)
{
	int k = parr.size() + 1;
	k = k > kn ? kn : k;
	int insp = k - 1;
	for (int i = 0; i < k - 1; i++) {
		if (parr[i].second > pa.second) {
			insp = i;
			break;
		}
	}
	if (parr.size() < kn) {
		parr.push_back(pad(0, 0));
	}
	for (int i = k - 1; i >= insp; i--) {
		if (i == insp) {
			parr[i] = pa;
		}
		else {
			parr[i] = parr[i - 1];
		}
	}
}

bool is_inbound(Point i, double *bound)
{
	if (i.x < bound[0] || i.x > bound[1]) {
		return false;
	}
	if (i.y < bound[2] || i.y > bound[3]) {
		return false;
	}
	if (i.z < bound[4] || i.z > bound[5]) {
		return false;
	}
	return true;
}

bool bound_inset(double *bound1, double * bound2) {
	if (bound1[1] < bound2[0] || bound1[0] > bound2[1]) {
		return false;
	}
	if (bound1[3] < bound2[2] || bound1[2] > bound2[3]) {
		return false;
	}
	if (bound1[5] < bound2[4] || bound1[4] > bound2[5]) {
		return false;
	}
	return true;
}

inline double Fo(Point a, Point b, int depth)
{
	double d = DEPLEN[depth - 1];
	double x = abs(a.x - b.x) / d;
	double y = abs(a.y - b.y) / d;
	double z = abs(a.z - b.z) / d;
	if (x > 1.5 || y > 1.5 || z > 1.5) {
		return 0;
	}
	double xv = 0.5 * pow(1.5 - x, 2);
	double yv = 0.5 * pow(1.5 - y, 2);
	double zv = 0.5 * pow(1.5 - z, 2);
	if (x < 0.5) {
		xv = (0.75 - pow(x, 2));
	}
	if (y < 0.5) {
		yv = (0.75 - pow(y, 2));
	}
	if (z < 0.5) {
		zv = (0.75 - pow(z, 2));
	}
	return xv * yv * zv;
}

inline double B_func(double x, int d)
{
	double dep = DEPLEN[0.5, d - 1];
	x = x / dep;
	if (abs(x) > 1.5) {
		return 0;
	}
	if (abs(x) < 0.5) {
		return 0.75 - pow(x, 2);
	}
	else
	{
		return 0.5 * pow(1.5 - abs(x), 2);
	}
}

inline double DB_func(double x, int d)
{
	double dep = DEPLEN[0.5, d - 1];
	x = x / dep;
	if (abs(x) > 1.5) {
		return 0;
	}
	if (abs(x) < 0.5) {
		return -2 * x / dep;
	}
	else if (x > 0)
	{
		return (x - 1.5) / dep;
	}
	else {
		return (x + 1.5) / dep;
	}
}

inline double guassInt_BB(double x1, int d1, double x2, int d2)
{
	if (x1 > x2) {
		std::swap(x1, x2);
		std::swap(d1, d2);
	}
	double dep1 = DEPLEN[d1];
	double dep2 = DEPLEN[d2];
	if (x1 + 3 * dep1 < x2 - 3 * dep2) {
		return 0;
	}
	double I[8] = { x1 - 3 * dep1, x1 - dep1, x1 + dep1, x1 + 3 * dep1,
					x2 - 3 * dep2, x2 - dep2, x2 + dep2, x2 + 3 * dep2 };

	std::sort(I, I + 8);
	double s = 0;
	for (int i = 0; i < 7; i++) {
		if (I[i] < x2 - 3 * dep2 || I[i + 1] > x1 + 3 * dep1) {
			continue;
		}
		double a = I[i];
		double b = I[i + 1];
		double l = (b - a) / 2;
		double p0 = (a + b) / 2;
		double p1 = p0 - l * sqrt(3.0 / 5);
		double p2 = p0 + l * sqrt(3.0 / 5);
		s += l * ((5.0 / 9) * B_func(p1 - x1, d1) * B_func(p1 - x2, d2) +
			(8.0 / 9) * B_func(p0 - x1, d1) * B_func(p0 - x2, d2) +
			(5.0 / 9) * B_func(p2 - x1, d1) * B_func(p2 - x2, d2));
	}
	return s;
}

inline double guassInt_DBB(double x1, int d1, double x2, int d2)
{
	int format = 0;
	if (x1 > x2) {
		std::swap(x1, x2);
		std::swap(d1, d2);
		format = 1;
	}
	double dep1 = DEPLEN[d1];
	double dep2 = DEPLEN[d2];
	if (x1 + 3 * dep1 < x2 - 3 * dep2) {
		return 0;
	}
	double I[8] = { x1 - 3 * dep1, x1 - dep1, x1 + dep1, x1 + 3 * dep1,
					x2 - 3 * dep2, x2 - dep2, x2 + dep2, x2 + 3 * dep2 };
	std::sort(I, I + 8);
	double s = 0;
	for (int i = 0; i < 7; i++) {
		if (I[i] < x2 - 3 * dep2 || I[i + 1] > x1 + 3 * dep1) {
			continue;
		}
		double a = I[i];
		double b = I[i + 1];
		double l = (b - a) / 2;
		double p0 = (a + b) / 2;
		double p1 = p0 - l * sqrt(1.0 / 3);
		double p2 = p0 + l * sqrt(1.0 / 3);
		if (!format) {
			s += l * (DB_func(p1 - x1, d1) * B_func(p1 - x2, d2) +
				DB_func(p2 - x1, d1) * B_func(p2 - x2, d2));
		}
		else {
			s += l * (B_func(p1 - x1, d1) * DB_func(p1 - x2, d2) +
				B_func(p2 - x1, d1) * DB_func(p2 - x2, d2));
		}
	}
	return s;
}

inline double guassInt_DBDB(double x1, int d1, double x2, int d2)
{
	if (x1 > x2) {
		std::swap(x1, x2);
		std::swap(d1, d2);
	}
	double dep1 = DEPLEN[d1];
	double dep2 = DEPLEN[d2];
	if (x1 + 3 * dep1 < x2 - 3 * dep2) {
		return 0;
	}
	double I[8] = { x1 - 3 * dep1, x1 - dep1, x1 + dep1, x1 + 3 * dep1,
					x2 - 3 * dep2, x2 - dep2, x2 + dep2, x2 + 3 * dep2 };
	std::sort(I, I + 8);
	double s = 0;
	for (int i = 0; i < 7; i++) {
		if (I[i] < x2 - 3 * dep2 || I[i + 1] > x1 + 3 * dep1) {
			continue;
		}
		double a = I[i];
		double b = I[i + 1];
		double l = (b - a) / 2;
		double p0 = (a + b) / 2;
		double p1 = p0 - l * sqrt(1.0 / 3);
		double p2 = p0 + l * sqrt(1.0 / 3);
		s += l * (DB_func(p1 - x1, d1) * DB_func(p1 - x2, d2) +
			DB_func(p2 - x1, d1) * DB_func(p2 - x2, d2));
	}
	return s;
}
void OctTree::set_gamma(double gam)
{
	gamma = gam;
}

void OctTree::get_outNode()
{
	outNode.clear();
	for (int i = 0; i < nopNode.size(); i++) {
		OctNode *node = nodeArr[nopNode[i]];
		Point o = node->center;
		if (o.x < pxMin || o.x > pxMax || o.y < pyMin || o.y > pyMax || o.z < pzMin || o.z > pzMax) {
			outNode.push_back(i);
		}
	}
}
void OctTree::get_allpoints(std::vector<int> &pid, int nid)
{
	OctNode *node = nodeArr[nid];
	if (node->isleaf) {
		for (int i = 0; i < node->inPointId.size(); i++) {
			pid.push_back(node->inPointId[i]);
		}
	}
	else {
		for (int i = 0; i < 8; i++) {
			get_allpoints(pid, node->child[i]);
		}
	}
}

void OctTree::calc_avgdis()
{
	avgdis.clear();
	avgdis.resize(points.size(), 0);

	for (int i = 0; i < points.size(); i++) {
		int K = 1;
		Point p = points[i];
		ANNpoint tp = annAllocPt(3);
		tp[0] = p.x;
		tp[1] = p.y;
		tp[2] = p.z;
		ANNidxArray nnIdx = new ANNidx[K + 1];
		ANNdistArray dists = new ANNdist[K + 1];
		kdTree->annkSearch(tp, K + 1, nnIdx, dists);
		for (int j = 0; j < K; j++) {
			avgdis[i] += dists[j + 1];
		}
		avgdis[i] = sqrt(avgdis[i] / K);
	}
}

void OctTree::get_knear(Point p, int k, std::vector<int> &pid)
{
	std::vector<pad> pais;
	std::stack<int> sta;
	sta.push(0);
	while (!sta.empty()) {
		OctNode *node = nodeArr[sta.top()];
		sta.pop();
		if (node->isleaf) {
			for (int i = 0; i < node->inPointId.size(); i++) {
				Point np = node->inPoints[i];
				double dis = pow(p.x - np.x, 2) + pow(p.y - np.y, 2) + pow(p.z - np.z, 2);
				if (pais.size() < k) {
					insertpad(pais, pad(node->inPointId[i], dis), k);
				}
				else {
					if (pais[k - 1].second > dis) {
						insertpad(pais, pad(node->inPointId[i], dis), k);
					}
				}
			}
		}
		else {
			std::vector<pad> npis;
			for (int i = 0; i < 8; i++) {
				OctNode *nodei = nodeArr[node->child[i]];
				Point oi = node->center;
				double d = DEPLEN[node->depth];
				double bound[6] = { oi.x - d, oi.x + d, oi.y - d, oi.y + d, oi.z - d, oi.z + d };
				
				double dis = 0;
				if (is_inbound(p, bound)) {
					dis = 0;
				}
				else {
					if (p.x < oi.x) {
						dis += pow(p.x - bound[0], 2);
					}
					else {
						dis += pow(p.x - bound[1], 2);
					}
					if (p.y < oi.y) {
						dis += pow(p.y - bound[2], 2);
					}
					else {
						dis += pow(p.y - bound[3], 2);
					}
					if (p.z < oi.z) {
						dis += pow(p.z - bound[4], 2);
					}
					else {
						dis += pow(p.z - bound[5], 2);
					}
				}
				insertpad(npis, pad(nodei->nid, dis), 8);
			}
			for (int i = 7; i >= 0; i--) {
				if (pais.size() == k && pais[k - 1].second < npis[i].second) {
					continue;
				}
				sta.push(npis[i].first);
			}
		}
	}
	pid.clear();
	for (int i = 0; i < pais.size(); i++) {
		pid.push_back(pais[i].first);
	}
}

void OctTree::find_around()
{
	for (int i = 0; i < leafNodeArr.size(); i++) {
		OctNode *nodei = leafNodeArr[i];
		Point oi = nodei->center;
		double d = DEPLEN[nodei->depth];
		double bound[6] = { oi.x - 3 * d, oi.x + 3 * d, oi.y - 3 * d, oi.y + 3 * d, oi.z - 3 * d, oi.z + 3 * d };

		std::stack<int> sta;
		sta.push(0);
		while (!sta.empty()) {
			int nid = sta.top();
			sta.pop();
			OctNode *nodej = nodeArr[nid];
			if (nodej->isleaf) {
				nodei->add_around(nodej->nid);
			}
			else {
				for (int k = 0; k < 8; k++) {
					OctNode *nodek = nodeArr[nodej->child[k]];
					Point ok = nodek->center;
					double dk = DEPLEN[nodek->depth];
					double boundk[6] = { ok.x - 3 * dk, ok.x + 3 * dk, ok.y - 3 * dk, ok.y + 3 * dk, ok.z - 3 * dk, ok.z + 3 * dk };
					if (bound_inset(boundk, bound)) {
						sta.push(nodek->nid);
					}
				}
			}
		}
	}
}


void OctTree::find_Ngbr()
{
	Ngbr.clear();
	Ngbr.resize(points.size(), std::vector<int>());

	int lN = leafNodeArr.size();
	ANNpointArray dataPts = annAllocPts(lN, 3);
	for (int i = 0; i < lN; i++) {
		Point o = leafNodeArr[i]->center;
		dataPts[i][0] = o.x;
		dataPts[i][1] = o.y;
		dataPts[i][2] = o.z;
	}

	if (kdTree_o) delete kdTree_o;
	kdTree_o = new ANNkd_tree(dataPts, lN, 3);

	for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		ANNpoint tp = annAllocPt(3);
		tp[0] = p.x;
		tp[1] = p.y;
		tp[2] = p.z;
		ANNidxArray nnIdx = new ANNidx[8];
		ANNdistArray dists = new ANNdist[8];
		kdTree_o->annkSearch(tp, 8, nnIdx, dists);
		for (int k = 0; k < 8; k++) {
			Ngbr[i].push_back(leafNodeArr[nnIdx[k]]->nid);
		}
	}

	/*for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		double d = DEPLEN[pdep[i] - 1];
		double bound[6] = { p.x - d, p.x + d, p.y - d, p.y + d, p.z - d, p.z + d };
		std::stack<int> sta;
		sta.push(0);
		while (!sta.empty()) {
			int nid = sta.top();
			sta.pop();
			OctNode *node = nodeArr[nid];
			if (node->isleaf) {
				if (node->depth >= pdep[i] && isNgbr(p, node->center, node->depth)) {
					Ngbr[i].push_back(node->nid);
				}
			}
			else {
				for (int j = 0; j < 8; j++) {
					OctNode *nodej = nodeArr[node->child[j]];
					if (is_intersect(nodej, bound)) {
						sta.push(nodej->nid);
					}
				}
			}
		}
	}*/
}



void OctTree::calc_B()
{
	std::vector<Td> tripletList_B[3];
	std::vector<Td> tripletList_Bw[3];
	wds.resize(points.size(), 0);
	for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		std::unordered_set<int> aroSet;
		std::vector<double> alpha(Ngbr[i].size(), 0);
		//double d = DEPLEN[pdep[i] - 1];
		for (int j = 0; j < Ngbr[i].size(); j++) {
			OctNode *nodej = nodeArr[Ngbr[i][j]];
			for (int k = 0; k < nodej->around.size(); k++) {
				if (!aroSet.count(nodej->around[k])) {
					aroSet.insert(nodej->around[k]);
				}
			}
			Point oj = nodej->center;
			//alpha[j] = pow(2, pdep[i] - 1) * (1 - abs(p.x - oj.x) / d) * (1 - abs(p.y - oj.y) / d) * (1 - abs(p.z - oj.z) / d);
			/*alpha[j] = pow(2, pdep[i] - 1);*/
			alpha[j] = 1;

		}
		for (auto it : aroSet) {
			OctNode *nodes = nodeArr[it];
			Point ps = nodes->center;
			int deps = nodes->depth;
			double wn[3] = { 0, 0, 0 };
			double kpsr = 0;
			for (int j = 0; j < Ngbr[i].size(); j++) {
				OctNode *nodej = nodeArr[Ngbr[i][j]];
				Point pj = nodej->center;
				int depj = nodej->depth;
				double DBBx = guassInt_DBB(ps.x, deps, pj.x, depj);
				double BBx = guassInt_BB(ps.x, deps, pj.x, depj);
				double DBBy = guassInt_DBB(ps.y, deps, pj.y, depj);
				double BBy = guassInt_BB(ps.y, deps, pj.y, depj);
				double DBBz = guassInt_DBB(ps.z, deps, pj.z, depj);
				double BBz = guassInt_BB(ps.z, deps, pj.z, depj);
				wn[0] += alpha[j] * DBBx * BBy * BBz;
				wn[1] += alpha[j] * BBx * DBBy * BBz;
				wn[2] += alpha[j] * BBx * BBy * DBBz;
				kpsr += sigma_g * Fo(ps, pj, depj);

				for (int t = 0; t < nodes->inPointId.size(); t++) {
					int ptid = nodes->inPointId[t];
					Point pt = nodes->inPoints[t];
					wds[ptid] += alpha[j] * Fo(pt, nodej->center, depj);
				}
			}
			//std::cout << wn[0] << std::endl;
			tripletList_B[0].push_back(Td(nodes->lnid, i, wn[0]));
			tripletList_B[1].push_back(Td(nodes->lnid, i, wn[1]));
			tripletList_B[2].push_back(Td(nodes->lnid, i, wn[2]));
		}
	}
	for (int i = 0; i < tripletList_B[0].size(); i++) {
		int row = tripletList_B[0][i].row();
		int col = tripletList_B[0][i].col();
		//std::cout <<  wds[col] << std::endl;
		if (abs(wds[col]) < 1e-16) {
			//std::cout << "ERROR" << i << " " << col << " " << wds[col] << std::endl;
			wds[col] = 1;
		}
		tripletList_Bw[0].push_back(Td(row, col, tripletList_B[0][i].value() / wds[col]));
		tripletList_Bw[1].push_back(Td(row, col, tripletList_B[1][i].value() / wds[col]));
		tripletList_Bw[2].push_back(Td(row, col, tripletList_B[2][i].value() / wds[col]));
	}
	B[0].resize(leafNodeArr.size(), points.size());
	B[0].setFromTriplets(tripletList_Bw[0].begin(), tripletList_Bw[0].end());

	B[1].resize(leafNodeArr.size(), points.size());
	B[1].setFromTriplets(tripletList_Bw[1].begin(), tripletList_Bw[1].end());

	B[2].resize(leafNodeArr.size(), points.size());
	B[2].setFromTriplets(tripletList_Bw[2].begin(), tripletList_Bw[2].end());

	B_T[0] = B[0].transpose();
	B_T[1] = B[1].transpose();
	B_T[2] = B[2].transpose();
}


void OctTree::calc_A()
{
	std::vector<Td> tripletList_A;
	for (int i = 0; i < leafNodeArr.size(); i++) {
		OctNode *nodei = leafNodeArr[i];
		Point pi = nodei->center;
		double depi = nodei->depth;
		for (int j = 0; j < nodei->around.size(); j++) {
			OctNode *nodej = nodeArr[nodei->around[j]];
			Point pj = nodej->center;
			int depj = nodej->depth;
			double DBDBx = guassInt_DBDB(pi.x, depi, pj.x, depj);
			double BBx = guassInt_BB(pi.x, depi, pj.x, depj);
			double DBDBy = guassInt_DBDB(pi.y, depi, pj.y, depj);
			double BBy = guassInt_BB(pi.y, depi, pj.y, depj);
			double DBDBz = guassInt_DBDB(pi.z, depi, pj.z, depj);
			double BBz = guassInt_BB(pi.z, depi, pj.z, depj);
		
			
			double w = DBDBx * BBy * BBz;
			w += BBx * DBDBy * BBz;
			w += BBx * BBy * DBDBz;
			//std::cout << w << std::endl;
			if (abs(w) > 1e-6) {
				tripletList_A.push_back(Td(nodei->lnid, nodej->lnid, w));
			}
		}

	}
	A.resize(leafNodeArr.size(), leafNodeArr.size());
	A.setFromTriplets(tripletList_A.begin(), tripletList_A.end());

	A_solver.compute(A);
}



void OctTree::calc_P()
{
	std::vector<Td> tripletList_p;
	for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		OctNode *nodei = nodeArr[pointNodeId[i]];
		for (int j = 0; j < nodei->around.size(); j++) {
			OctNode *nodej = nodeArr[nodei->around[j]];
			double v = Fo(p, nodej->center, nodej->depth);
			if (abs(v) > 1e-12) {
				tripletList_p.push_back(Td(i, nodej->lnid, v));
			}
		}
	}
	P.resize(points.size(), leafNodeArr.size());
	P.setFromTriplets(tripletList_p.begin(), tripletList_p.end());

	std::vector<Td> tripletList_p1;
	std::vector<Td> tripletList_p2;
	int cnt = 0;
	int outn = outNode.size();
	for (int i = 0; i < nopNode.size(); i++) {
		OctNode *nodei = nodeArr[nopNode[i]];
		Point ci = nodei->center;
		if (cnt < outn && outNode[cnt] == i) {
			for (int j = 0; j < nodei->around.size(); j++) {
				OctNode *nodej = nodeArr[nodei->around[j]];
				double v = Fo(ci, nodej->center, nodej->depth);
				if (abs(v) > 1e-12) {
					tripletList_p1.push_back(Td(i, nodej->lnid, v));
					tripletList_p2.push_back(Td(cnt, nodej->lnid, v));
				}
			}
			cnt++;
		}
		else {
			for (int j = 0; j < nodei->around.size(); j++) {
				OctNode *nodej = nodeArr[nodei->around[j]];
				double v = Fo(ci, nodej->center, nodej->depth);
				if (abs(v) > 1e-12) {
					tripletList_p1.push_back(Td(i, nodej->lnid, v));
				}
			}
		}
	}
	P1.resize(nopNode.size(), leafNodeArr.size());
	P1.setFromTriplets(tripletList_p1.begin(), tripletList_p1.end());

	P2.resize(outn, leafNodeArr.size());
	P2.setFromTriplets(tripletList_p2.begin(), tripletList_p2.end());

	P_T = P.transpose();
	P1_T = P1.transpose();
}


void OctTree::find_near_points()
{
	int N = points.size();
	int nopN = nopNode.size();
	near_points.clear();
	near_points.resize(nopNode.size(), std::vector<double>());
	near_points_id.clear();
	near_points_id.resize(nopNode.size(), std::vector<int>());

	for (int i = 0; i < nopNode.size(); i++) {
		OctNode *nodei = nodeArr[nopNode[i]];
		Point o = nodei->center;
		ANNpoint tp = annAllocPt(3);
		tp[0] = o.x;
		tp[1] = o.y;
		tp[2] = o.z;
		ANNidxArray nnIdx = new ANNidx[k_n];
		ANNdistArray dists = new ANNdist[k_n];
		kdTree->annkSearch(tp, k_n, nnIdx, dists);
		for (int k = 0; k < k_n; k++) {
			int pid = nnIdx[k];
			near_points_id[i].push_back(pid);
			Point p = points[pid];
			double norm = sqrt(pow(o.x - p.x, 2) + pow(o.y - p.y, 2) + pow(o.z - p.z, 2));
			near_points[i].push_back((o.x - p.x) / norm);
			near_points[i].push_back((o.y - p.y) / norm);
			near_points[i].push_back((o.z - p.z) / norm);
		}
	}

	int cnt = 0;
	ANNpointArray dataPts = annAllocPts(nopNode.size(), 3);
	for (int i = 0; i < nopNode.size(); i++) {
		OctNode *nodei = nodeArr[nopNode[i]];
		if (nodei->depth == APTDEPTH) {
			Point o = nodei->center;
			dataPts[cnt][0] = o.x;
			dataPts[cnt][1] = o.y;
			dataPts[cnt][2] = o.z;
			cnt++;
		}
	}

	if (kdTree_no) delete kdTree_no;
	kdTree_no = new ANNkd_tree(dataPts, cnt, 3);
	for (int i = 0; i < points.size(); i++) {
		Point p = points[i];
		ANNpoint tp = annAllocPt(3);
		tp[0] = p.x;
		tp[1] = p.y;
		tp[2] = p.z;
		ANNidxArray nnIdx = new ANNidx[k_p];
		ANNdistArray dists = new ANNdist[k_p];
		kdTree_no->annkSearch(tp, k_p, nnIdx, dists);
		for (int k = 0; k < k_p; k++) {
			int npid = nnIdx[k];
			near_points_id[npid].push_back(i);
			OctNode* nodek = nodeArr[nopNode[npid]];
			Point o = nodek->center;
			double norm = sqrt(pow(o.x - p.x, 2) + pow(o.y - p.y, 2) + pow(o.z - p.z, 2));
			double al = DEPLEN[APTDEPTH - 1];
			al = norm;

			near_points[npid].push_back((o.x - p.x) / norm);
			near_points[npid].push_back((o.y - p.y) / norm);
			near_points[npid].push_back((o.z - p.z) / norm);
		}
	}
}

void func_callback_w_uv(const size_t N, const std::vector<double>& x, double& f, std::vector<double>& g, void* user_supply)
{
	auto ptr_this = static_cast<OctTree*>(user_supply);
	ptr_this->energy_evaluation_w_uv(x, f, g);
}

void OctTree::optimize_LBFS_w_uv(std::vector<double>& x)
{
	int N = points.size();
	if (2 * N != x.size()) {
		std::cout << "x size error" << std::endl;
	}
	HLBFGS solver;
	int it = 0;
	solver.set_number_of_variables(2 * N);
	solver.set_verbose(true);
	solver.set_func_callback(func_callback_w_uv, 0, 0, 0, 0);
	solver.optimize_without_constraints(&(x[0]), 10000, this, &it);

	clear_near_points();
	set_gamma(0);
	sig = 0.005;
	solver.set_Info(it, 1);
	solver.optimize_without_constraints(&(x[0]), 10000, this, &it);
	sig = 0.002;
	solver.set_Info(it, 1);
	solver.optimize_without_constraints(&(x[0]), 10000, this, &it);
}
void OctTree::clear_near_points()
{
	for (int i = 0; i < nopNode.size(); i++) {
		while (near_points_id[i].size() > k_n) {
			near_points_id[i].pop_back();
			near_points[i].pop_back();
			near_points[i].pop_back();
			near_points[i].pop_back();
		}
	}
}
void OctTree::energy_evaluation_w_uv(const std::vector<double>& x, double& f, std::vector<double>& g)
{
	int N = points.size();
	int leafN = leafNodeArr.size();
	int nopN = nopNode.size();
	f = 0;
	g.clear();
	g.resize(2 * N, 0);
	Eigen::VectorXd normal[3];
	normal[0].resize(N);
	normal[1].resize(N);
	normal[2].resize(N);
	for (int i = 0; i < N; i++) {
		normal[0](i) = sin(x[2 * i]) * cos(x[2 * i + 1]);
		normal[1](i) = sin(x[2 * i]) * sin(x[2 * i + 1]);
		normal[2](i) = cos(x[2 * i]);
	}
	Eigen::VectorXd tempbv = B[0] * normal[0] + B[1] * normal[1] + B[2] * normal[2];
	Eigen::VectorXd alp = A_solver.solve(tempbv);
	Eigen::VectorXd mu = P * alp;
	Eigen::VectorXd mu1 = P1 * alp;
	std::cout << "mean: " << mu.mean()  << std::endl;
	mu1 = mu1 - Eigen::VectorXd::Constant(nopN, mu.mean());
	Eigen::VectorXd mu1_unit = mu1 / sig;
	Eigen::VectorXd Fdmu = Eigen::VectorXd::Constant(nopN, 0);
	int oi = 0;
	for (int i = 0; i < nopN; i++) {
		double ln = 0;
		int dep = nodeArr[nopNode[i]]->depth;
		for (int n = 0; n < near_points_id[i].size(); n++) {
			int pid = near_points_id[i][n];
			ln += near_points[i][3 * n] * normal[0](pid) +
				near_points[i][3 * n + 1] * normal[1](pid) + near_points[i][3 * n + 2] * normal[2](pid);
		}
		double V = pow(2, nodeArr[nopNode[i]]->depth - APTDEPTH);
		V = pow(DEPLEN[dep - 1], 3);
		ln *= V;
		double p = 0.5 * erf(-mu1_unit(i) / sqrt(2));
		if (i == outNode[oi]) {
			p = -0.5;
		}
		f += ln * p;
		double dpdmu = (-1 / (sqrt(2 * PI) * sig)) * exp(-0.5 * pow(mu1_unit(i), 2));
		if (i == outNode[oi]) {
			dpdmu = 0;
			oi++;
		}
		Fdmu(i) = ln * dpdmu;

		for (int n = 0; n < near_points_id[i].size(); n++) {
			int pid = near_points_id[i][n];
			double sinuj = sin(x[2 * pid]);
			double cosuj = cos(x[2 * pid]);
			double sinvj = sin(x[2 * pid + 1]);
			double cosvj = cos(x[2 * pid + 1]);
			g[2 * pid] += V * p * (near_points[i][3 * n] * cosuj * cosvj +
				near_points[i][3 * n + 1] * cosuj * sinvj + near_points[i][3 * n + 2] * -sinuj);
			g[2 * pid + 1] += V * p * (near_points[i][3 * n] * -sinuj * sinvj + near_points[i][3 * n + 1] * sinuj * cosvj);
		}
	}

	std::cout << "F: " << f << " ";

	double f0 = f;
	double fmu = mu.mean();
	//f += gamma * sf * sf;
	Eigen::VectorXd Gdmu = Eigen::VectorXd::Constant(N, 0);
	double gamma_ = gamma / N;
	for (int i = 0; i < N; i++) {
		double d = mu(i) - fmu;
		f += gamma_ * d * d;
		Gdmu(i) = 2 * gamma_ * d;
	}
	std::cout << "G: " << f - f0 << " ";

	Eigen::VectorXd Edx = P1_T * Fdmu + P_T * Gdmu;
	Eigen::VectorXd y = A_solver.solve(Edx);
	//std::cout << "A solver iterations2: " << A_solver.iterations() << std::endl;
	Eigen::VectorXd Edn[3] = { B_T[0] * y, B_T[1] * y, B_T[2] * y };
	for (int j = 0; j < N; j++) {
		double sinuj = sin(x[2 * j]);
		double cosuj = cos(x[2 * j]);
		double sinvj = sin(x[2 * j + 1]);
		double cosvj = cos(x[2 * j + 1]);
		g[2 * j] += Edn[0](j) * cosuj * cosvj + Edn[1](j) * cosuj * sinvj + Edn[2](j) * -sinuj;
		g[2 * j + 1] += Edn[0](j) * -sinuj * sinvj + Edn[1](j) * sinuj * cosvj;
	}
}

void OctTree::energy_evaluation_w_print_p(const std::vector<double>& x, std::vector<double>& p)
{
	int N = points.size();
	int leafN = leafNodeArr.size();
	int nopN = nopNode.size();
	p.resize(nopN);
	Eigen::VectorXd normal[3];
	normal[0].resize(N);
	normal[1].resize(N);
	normal[2].resize(N);
	for (int i = 0; i < N; i++) {
		normal[0](i) = sin(x[2 * i]) * cos(x[2 * i + 1]);
		normal[1](i) = sin(x[2 * i]) * sin(x[2 * i + 1]);
		normal[2](i) = cos(x[2 * i]);
	}
	Eigen::VectorXd tempbv = B[0] * normal[0] + B[1] * normal[1] + B[2] * normal[2];
	Eigen::VectorXd alp = A_solver.solve(tempbv);
	Eigen::VectorXd mu = P * alp;
	Eigen::VectorXd mu1 = P1 * alp;
	mu1 = mu1 - Eigen::VectorXd::Constant(nopN, mu.mean());
	Eigen::VectorXd mu1_unit = mu1 / sig;

	for (int i = 0; i < nopN; i++) {
		p[i] = 0.5 * erf(-mu1_unit(i) / sqrt(2)) + 0.5;
	}
}