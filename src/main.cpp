#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include "pointCloud.h"
#include "octTree.h"

void writeResult(std::string filePath, std::vector<Point> points, std::vector<double> normal)
{
	std::ofstream fout;
	fout.open(filePath, std::ios::out);
	if (!fout.is_open()) {
		if (!fout.is_open()) {
			std::cout << filePath << "File write error" << std::endl;
			return;
		}
	}
	/* ply
	fout << "ply" << std::endl;
	fout << "format ascii 1.0" << std::endl;
	fout << "element vertex " << points.size() << std::endl;
	fout << "property double x" << std::endl;
	fout << "property double y" << std::endl;
	fout << "property double z" << std::endl;
	fout << "property double nx" << std::endl;
	fout << "property double ny" << std::endl;
	fout << "property double nz" << std::endl;
	fout << "end_header" << std::endl;
	*/
	for (int i = 0; i < points.size(); i++) {
		fout << points[i].x << " ";
		fout << points[i].y << " ";
		fout << points[i].z << " ";
		fout << normal[3 * i] << " ";
		fout << normal[3 * i + 1] << " ";
		fout << normal[3 * i + 2] << std::endl;
	}
	fout.close();
}

void writeResult(std::string filePath, std::vector<Point> points)
{
	std::ofstream fout;
	fout.open(filePath, std::ios::out);
	if (!fout.is_open()) {
		if (!fout.is_open()) {
			std::cout << filePath << "File write error" << std::endl;
			return;
		}
	}
	
	for (int i = 0; i < points.size(); i++) {
		fout << points[i].x << " ";
		fout << points[i].y << " ";
		fout << points[i].z << std::endl;
	}
	fout.close();
}
void writeOctreeNode(std::string filePath, std::vector<Point> points, OctTree &octree, std::vector<double> &normal)
{
	std::ofstream fout;
	fout.open(filePath, std::ios::out);
	if (!fout.is_open()) {
		if (!fout.is_open()) {
			std::cout << filePath << "File write error" << std::endl;
			return;
		}
	}
	std::vector<double> p;
	octree.energy_evaluation_w_print_p(normal, p);
	int N = octree.points.size();
	fout << "ply" << std::endl;
	fout << "format ascii 1.0" << std::endl;
	fout << "element vertex " << points.size() << std::endl;
	fout << "property double x" << std::endl;
	fout << "property double y" << std::endl;
	fout << "property double z" << std::endl;
	fout << "property uint8 red" << std::endl;
	fout << "property uint8 green" << std::endl;
	fout << "property uint8 blue" << std::endl;
	fout << "end_header" << std::endl;
	for (int i = 0; i < points.size(); i++) {
		fout << points[i].x << " ";
		fout << points[i].y << " ";
		fout << points[i].z << " ";
		if (i < N) {
			fout << 255 << " ";
			fout << 0 << " ";
			fout << 0 << std::endl;
		}
		else {
			fout << (int)((p[i - N]) * 255) << " ";
			fout << (int)((1 - p[i - N]) * 255) << " ";
			fout << (int)((p[i - N]) * 255) << std::endl;
		}
	}

	fout.close();
}

int main(int argc, char *argv[])
{
	std::vector<double> normal;
	std::vector<Point> points;
	//std::string fileName = argv[1];
	//std::string dataPath = "..\\..\\data\\";
	//std::string resultPath = "..\\..\\result\\";
	std::string fileName = "input";
	std::string dataPath = "..\\data\\";
	std::string resultPath = "..\\result\\";
	std::string type = ".xyz";
	PointCloud pointcloud(dataPath + fileName + type);
	int N = pointcloud.points.size();

	Eigen::VectorXd Nx = Eigen::VectorXd::Constant(N, 0);
	Eigen::VectorXd Ny = Eigen::VectorXd::Constant(N, 0);
	Eigen::VectorXd Nz = Eigen::VectorXd::Constant(N, 0);
	
	std::vector<double> uv(2 * N);
	std::vector<double> nor(3 * N);
	srand(time(0));
	//initialize
	for (int i = 0; i < N; i++) {
		uv[2 * i] = rand();
		uv[2 * i + 1] = rand();
	}
	for (int i = 0; i < N; i++) {
		nor[3 * i] = sin(uv[2 * i]) * cos(uv[2 * i + 1]);
		nor[3 * i + 1] = sin(uv[2 * i]) * sin(uv[2 * i + 1]);
		nor[3 * i + 2] = cos(uv[2 * i]);
	}

	OctTree octree(&pointcloud);

	if (DEBUG) {
		for (int i = 0; i < pointcloud.points.size(); i++) {
			points.push_back(octree.points[i]);
		}

		for (int i = 0; i < octree.nopNode.size(); i++) {
			points.push_back(octree.nodeArr[octree.nopNode[i]]->center);
		}
		//writeOctreeNode(resultPath + fileName + "_grid_p_start" + ".ply", points, octree, uv);
	}

	octree.optimize_LBFS_w_uv(uv);

	for (int i = 0; i < N; i++) {
		nor[3 * i] = sin(uv[2 * i]) * cos(uv[2 * i + 1]);
		nor[3 * i + 1] = sin(uv[2 * i]) * sin(uv[2 * i + 1]);
		nor[3 * i + 2] = cos(uv[2 * i]);
	}
	if (DEBUG) {
		writeOctreeNode(resultPath + fileName + "_grid_p_end" + ".ply", points, octree, uv);
	}

	writeResult(resultPath + fileName + "_res" + type, pointcloud.points, nor);
	return 0;
}

