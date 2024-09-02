#pragma once
#ifndef BASIC_STRUCTURE_H
#define BASIC_STRUCTURE_H
#define PI 3.141592654
constexpr double BSCALE = 1.2;
constexpr int MAXDEPTH = 10;
constexpr int APTDEPTH = 6;
constexpr int NPN = 1;
constexpr double DEPLEN[13] = { 1.0, 0.5, 0.25, 0.125, 0.0625,
								0.03125, 0.015625, 0.0078125, 0.00390625,
								0.001953125, 0.0009765625,
								0.00048828125, 0.000244140625 };
constexpr double sigma_g = 0.02;
constexpr int k_n = 7;
constexpr int k_p = 2;
constexpr bool DEBUG = true;
struct Point {
	double x;
	double y;
	double z;
	Point() {
		x = 0; y = 0; z = 0;
	}
	Point(double xx, double yy, double zz) {
		x = xx;
		y = yy;
		z = zz;
	}
	const Point operator+(const Point &point) const {
		return Point(x + point.x, y + point.y, z + point.z);
	}

	const Point operator-(const Point &point) const {
		return Point(x - point.x, y - point.y, z - point.z);
	}

	const Point operator*(const double &a) const {
		return Point(x * a, y * a, z * a);
	}

};

struct BoundingBox {
	double blx;		// bottom-left x
	double bly;		// bottom-left y
	double blz;		// bottom-left z
	double pxMax;
	double pxMin;
	double pyMax;
	double pyMin;
	double pzMax;
	double pzMin;
	double scale;
};

#endif


