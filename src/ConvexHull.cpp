/*
 * Convex hull algorithm - Library (C++)
 *
 * Copyright (c) 2017 Project Nayuki
 * https://www.nayuki.io/page/convex-hull-algorithm
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING.txt and COPYING.LESSER.txt).
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include "Geometry.h"

using std::vector;




vector<Point> makeConvexHullPresorted(const vector<Point> &points) {
	if (points.size() <= 1)
		return vector<Point>(points);

	// Andrew's monotone chain algorithm. Positive y coordinates correspond to "up"
	// as per the mathematical convention, instead of "down" as per the computer
	// graphics convention. This doesn't affect the correctness of the result.

	vector<Point> upperHull;
	for (const Point &p : points) {
		while (upperHull.size() >= 2) {
			const Point &q = *(upperHull.cend() - 1);  // Same as .back()
			const Point &r = *(upperHull.cend() - 2);
			if ((q.x - r.x) * (p.y - r.y) >= (q.y - r.y) * (p.x - r.x))
				upperHull.pop_back();
			else
				break;
		}
		upperHull.push_back(p);
	}
	upperHull.pop_back();

	vector<Point> lowerHull;
	for (vector<Point>::const_reverse_iterator it = points.crbegin(); it != points.crend(); ++it) {
		const Point &p = *it;
		while (lowerHull.size() >= 2) {
			const Point &q = *(lowerHull.cend() - 1);  // Same as .back()
			const Point &r = *(lowerHull.cend() - 2);
			if ((q.x - r.x) * (p.y - r.y) >= (q.y - r.y) * (p.x - r.x))
				lowerHull.pop_back();
			else
				break;
		}
		lowerHull.push_back(p);
	}
	lowerHull.pop_back();

	if (!(upperHull.size() == 1 && upperHull == lowerHull))
		upperHull.insert(upperHull.end(), lowerHull.cbegin(), lowerHull.cend());
	return upperHull;
}


vector<Point> makeConvexHull(const vector<Point> &points)
{
	vector<Point> newPoints = points;
	std::sort(newPoints.begin(), newPoints.end());
	return makeConvexHullPresorted(newPoints);
}

vector<Point> CombineHull(const vector<Point> &points_1, const vector<Point> &points_2)
{
	vector<Point> newPoints = points_1;
	newPoints.insert (newPoints.end(), points_2.begin(), points_2.end());
	std::sort(newPoints.begin(), newPoints.end());
	return makeConvexHullPresorted(newPoints);
}

double HullArea(const vector<Point> &points)
{
double A = 0.0;
int n = points.size();

    if ( n < 3) return 0.0;
    for (int i = 0; i < n - 1; i++)
    {
        A += points[i].x*points[i+1].y - points[i+1].x*points[i].y;
    }
    //  Then close the polygon, last point is the first one.
    A += points[n-1].x*points[0].y - points[0].x*points[n-1].y;
    A = fabs(A);
	return A/2;
}

Rectangle Hull2BoundingBox(const vector<Point> &points)
{
double xmin, xmax, ymin, ymax;
int n = points.size();
Rectangle BBox;

    xmax = xmin = points[0].x;
    ymax = ymin = points[0].y;
    for (int i = 1; i < n - 1; i++)
    {
        xmin = fmin(xmin, points[i].x);
        ymin = fmin(ymin, points[i].y);
        xmax = fmax(xmax, points[i].x);
        ymax = fmax(ymax, points[i].y);
    }
    BBox.A = Point(xmin, ymin);
    BBox.B = Point(xmax, ymax);
	return BBox;
}

