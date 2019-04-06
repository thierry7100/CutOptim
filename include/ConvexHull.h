#ifndef CONVEXHULL_H_INCLUDED
#define CONVEXHULL_H_INCLUDED


vector<Point> makeConvexHull(const vector<Point> &points);
vector<Point> CombineHull(const vector<Point> &points_1, const vector<Point> &points_2);
double HullArea(const vector<Point> &points);
Rectangle Hull2BoundingBox(const vector<Point> &points);


#endif // CONVEXHULL_H_INCLUDED
