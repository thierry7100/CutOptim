#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
using namespace std;

const double my_Default_Precision = 1e-4;
const double angle_precision = 1e-6;

enum PositionType {
	TopLeft = 0,
	TopCenter = 1,
	TopRight = 2,
	CenterLeft = 4,
	CenterCenter = 5,
	CenterRight = 6,
    BottomLeft = 8,
    BottomCenter = 9,
    BottomRight = 10
};


class Point
{
public:
	double x;
	double y;
	Point (double newx=0,double newy=0)
	{
		x=newx;
		y=newy;
	}
	/*
		Overloaded == operator to check for equality between 2 objects of class point
	*/
	inline friend bool operator== (const Point& p1,const Point& p2)
	{
		return (p1.x==p2.x && p1.y==p2.y);
	}
	/*
		Overloaded != operator to check for non-equality between 2 objects of class point
	*/
	inline friend bool operator!= (const Point& p1,const Point& p2)
	{
		return (!(p1.x==p2.x && p1.y==p2.y));
	}

	/*
		Overloaded != operator to check if a point is "lower"
		Return true if X is lower. if X is equal return true if Y is lower
	*/
	inline friend bool operator< (const Point& p1,const Point& p2)
	{
        if ( p1.x != p2.x ) return(p1.x < p2.x);
        return(p1.y < p2.y);
	}
	/*
		Overloaded != operator to check if a point is "lower or equal"
		Return true if X is lower. if X is equal return true if Y is lower or equal
	*/
	inline friend bool operator<= (const Point& p1,const Point& p2)
	{
        if ( p1.x != p2.x ) return(p1.x < p2.x);
        return(p1.y <= p2.y);
	}
	/*
		Overloaded != operator to check if a point is "greater"
		Return true if X is greater. if X is equal return true if Y is greater
	*/
	inline friend bool operator> (const Point& p1,const Point& p2)
	{
        if ( p1.x != p2.x ) return(p1.x > p2.x);
        return(p1.y > p2.y);
	}

	/*
		Overloaded != operator to check if a point is "greater or equal"
		Return true if X is greater. if X is equal return true if Y is greater
	*/
	inline friend bool operator>= (const Point& p1,const Point& p2)
	{
        if ( p1.x != p2.x ) return(p1.x > p2.x);
        return(p1.y >= p2.y);
	}

	/*
		Overloaded + operator to returna point with coordinates which are the difference between the 2 points
	*/
	friend Point operator+ (const Point& p1,const Point& p2)
	{
		return (Point(p1.x+p2.x , p1.y+p2.y));
	}
	/*
		Overloaded + operator to returna point with coordinates which are the sum of the 2 points
	*/
	inline friend Point operator- (const Point& p1,const Point& p2)
	{
		return (Point(p1.x-p2.x , p1.y-p2.y));
	}
	/*
		Overloaded * operator multiplication by a scalar
	*/
	inline friend Point operator* (double Coeff,const Point& p2)
	{
		return (Point(Coeff*p2.x , Coeff*p2.y));
	}
	/*
		Overloaded ostream << operator to check for print object of class point to STDOUT
	*/
	friend ostream& operator<<(ostream& output,const Point& p)
	{
		output<<"("<<p.x<<","<<p.y<<")";
		return output;
	}

	inline double sqrDistance(Point &p2)
	{
        return( (x - p2.x)*(x - p2.x) + (y - p2.y)*(y - p2.y));
	}

	inline double Distance(Point &p2)
	{
        return( sqrt((x - p2.x)*(x - p2.x) + (y - p2.y)*(y - p2.y)));
	}

	inline void align_grid() { x = round(x/my_Default_Precision)*my_Default_Precision; y = round(y/my_Default_Precision)*my_Default_Precision;}

	inline void Transform(double angle, Point Centre, Point Translation)
	{
        double cos_a = cos(angle);
        double sin_a = sin(angle);
        double newx = (x - Centre.x)*cos_a - (y-Centre.y)*sin_a + Centre.x + Translation.x;
        double newy = (x - Centre.x)*sin_a + (y-Centre.y)*cos_a + Centre.y + Translation.y;
        newx = round(newx/my_Default_Precision)*my_Default_Precision;
        newy = round(newy/my_Default_Precision)*my_Default_Precision;
        x = newx;
        y = newy;
	}
};

class Rectangle
{
    public:
    Point A, B;

    double get_area() { return fabs((B.x - A.x)*(B.y - A.y)); }
};

/**
*   This simple class deals with lines.
*   A line is represented by its equation a.x + b.y + c = 0
*/
class Line
{
    protected:
        double a, b, c;

    public:
        Line(double la, double lb, double lc)
        {
            a = la;
            b = lb;
            c = lc;
        }
        Line()
        {
            a = 0;
            b = 0;
            c = 0;
        }

        /**
        *   Return the square of distance from the point p to the line itself.
        */
        inline double sqrDistancePoint(Point &p) const
        {
            return((a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b));
        }

        /**
        * Return the square of distance from the point p to the line itself, same as previous one but add the square root computation
        */

        inline double DistancePoint(Point &p) const
        {
            return(fabs(a * p.x + b * p.y + c)/sqrt(a*a + b*b));
        }

        /**
        *   Compute the intersection poi t between 2 lines.
        *   If lines are parallel, return the Point (nan, nan)
        */
        inline Point Intersect(const Line &Line2) const
        {
            double det = Line2.a * b - a*Line2.b;
            if ( fabs(det) < my_Default_Precision )       //  Line are parallel
            {
                return(Point(nan(""), nan("")));
            }
            return(Point((Line2.b*c - Line2.c*b)/det, (a*Line2.c - Line2.a*c)/det));
        }

        /**
        *   Return true if the line intersects with Line2
        *   If the intersection exists, update the intersection Point
        */
        inline bool hasIntersect(const Line *Line2, Point &ResIntersect) const
        {
            double det = Line2->a * b - a*Line2->b;
            if ( fabs(det) < my_Default_Precision )       //  Line are parallel
            {
                return(0);          //  Parallel lines, no intersection
            }
            ResIntersect = Point((Line2->b*c - Line2->c*b)/det, (a*Line2->c - Line2->a*c)/det);
            return 1;
        }
};

/**
*   This class is derived from line, but adds two ends
*/
class Segment:Line
{
    protected:
       double xA, xB,yA, yB;
    public:

        double xm, xM, ym, yM;

        Segment()
        {
        }

        Segment(Point &A, Point &B)
        {
            xA = A.x;
            xB = B.x;
            yA = A.y;
            yB = B.y;
            xm = fmin(xA, xB) - my_Default_Precision;
            xM = fmax(xA, xB) + my_Default_Precision;
            ym = fmin(yA, yB) - my_Default_Precision;
            yM = fmax(yA, yB) + my_Default_Precision;
            a = A.y - B.y;
            b = B.x - A.x;
            c = A.x * B.y - B.x * A.y;
        }
        inline Point getPointA()
        {
            return(Point(xA, yA));
        }
        inline Point getPointB()
        {
            return(Point(xB, yB));
        }
        inline void setPointA(Point &A)
        {
            xA = A.x;
            yA = A.y;
            a = yA - yB;
            b = xB - xA;
            c = xA*yB - xB*yA;
        }
        inline void setPointB(Point &B)
        {
            xB = B.x;
            yB = B.y;
            a = yA - yB;
            b = xB - xA;
            c = xA*yB - xB*yA;
        }

        // Compute segment size
        inline double SegmentSize()
        {
            return(sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB)));
        }

        /**
         * Return true if point p is in segment. p should be on the line
        */
        inline int InSegment(Point &p) const
        {
            if ( p.x < xm ) return 0;     //  Impossible lower than xmin
            if ( p.x > xM ) return 0;     //  Impossible greater than xmax
            if ( p.y < ym) return 0;     //  Impossible lower than ymin
            if ( p.y > yM ) return 0;     //  Impossible greater than ymax
            return(1);      //  OK
        }

        /**
        *  Return true if point p is in segment. p should be on the line, segment ends are not taken into account
        *   So return false if the point is an end.
        */
        inline int InSegmentNoEnd(Point &p) const
        {
            if ( p.x < xm ) return 0;     //  Impossible lower than xmin
            if ( p.x > xM ) return 0;     //  Impossible greater than xmax
            if ( p.y < ym) return 0;     //  Impossible lower than ymin
            if ( p.y > yM ) return 0;     //  Impossible greater than ymax
            //  Now checks if p is an end of the segment
            double d1 = fabs(p.x - xA) + fabs(p.y - yA);
            if ( d1 < 2*my_Default_Precision) return 0;                       //  End of segment return 0
            double d2 = fabs(p.x - xB) + fabs(p.y - yB);
            if ( d2 < 2*my_Default_Precision) return 0;                       //  End of segment return 0
            return(1);      //  OK
        }

        inline double sqrDistancePoint(const Point &p) const
        {
            //  If dot product AB.AC si < 0 return distance(A, C)
            //  Indeed, if this is the case C is in the half plane which is outer segment
            if ( (xB - xA)*(p.x - xA) + (yB-yA)*(p.y-yA) < 0 )
                return( (p.x - xA)*(p.x- xA) + (p.y-yA)*(p.y-yA));
            //  If dot product BA.BC si > 0 return distance(A, C)
            if ( (xA - xB)*(p.x - xB) + (yA-yB)*(p.y-yB) < 0 )
                return( (p.x - xB)*(p.x- xB) + (p.y-yB)*(p.y-yB));
            //  Distance point to line is a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b
            double d1 = (a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b);
            return d1;
        }
        //  Check if segment this  cross segment s1
        //  In any case, set ResIntersect to the point of intersection of the two lines
        //  If this point belongs to the 2 segments, return true
        //  If one segment is included in another, set the bool IncludedSegment but return false
        inline int isCrossing(const Segment *s1, Point *ResIntersect, bool *IncludedSegment=nullptr)
        {
            Point Ai;
            if ( IncludedSegment != nullptr ) *IncludedSegment = false;     //  Set default answer
            if ( !hasIntersect(s1, Ai) )
            {
                //  Check if same segment (same line...)
                if ( IncludedSegment != nullptr && fabs(s1->a - a) < my_Default_Precision && fabs(s1->b - b) < my_Default_Precision && fabs(s1->c - c) < my_Default_Precision )
                    *IncludedSegment = true;
                return 0;     // No intersection,
            }
            if ( ResIntersect != NULL) *ResIntersect = Ai;
            return InSegment(Ai) && s1->InSegment(Ai);
        }
        //  Check if segment this  cross segment s1
        //  In any case, set ResIntersect to the point of intersection of the two lines
        //  If this point belongs to the 2 segments, return true
        //  If one segment is included in another, set the bool IncludedSegment but return false
        inline int isCrossingNoEnd(const Segment *s1, Point *ResIntersect, bool *IncludedSegment=nullptr) const
        {
            Point Ai;
            if ( IncludedSegment != nullptr ) *IncludedSegment = false;     //  Set default answer
            if ( !hasIntersect(s1, Ai) )
            {
                //  Check if same segment (same line...), that is A point of this is on s1
                if ( IncludedSegment != nullptr && fabs(s1->a*xA + s1->b*yA  + s1->c) < my_Default_Precision )
                    *IncludedSegment = true;
                return 0;     // No intersection,
            }
            if ( ResIntersect != NULL) *ResIntersect = Ai;
            return InSegmentNoEnd(Ai) && s1->InSegmentNoEnd(Ai);
        }

    /*
		Overloaded ostream << operator to check for print object of class point to STDOUT
	*/
	friend ostream& operator<<(ostream& output,const Segment& Seg)
	{
		output<<"Segment("<<Seg.xA<<","<<Seg.yA<<") to " << Seg.xB << "," << Seg.yB << ") )";
		return output;
	}

};

enum PolyType {
	EmptyPoly = 0,
	VectorPoly = 1,
	RotatedPoly = 2,
	TranslatedPoly = 3,
	TransformedPoly = 4,
	CachedPoly = 5
};

class Polygon
{
    public:
        Polygon();
        Polygon(vector<Point> &points);
        virtual ~Polygon();
        void addVertice(Point &p);
        double area();
        void setPrecision(double epsilon);
        double distance(Point &p);

        int nVertices;
        inline int isClockWise()
        {
           //  1 : clockwise, -1 counter clockwise, 0 not computed yet
           return(poly_clockwise);
        }
        inline int isClosed()                       //  Return true if polygon is closed, that is first and last oint are the same
        {
            return poly_closed;
        }
        int isInPoly(const Point *p, int *SegmentIdx = NULL);
        int Poly_in_Poly(Polygon *Big);

        Rectangle GetBoundingBox();
        void ReCalcBoundingBox();
        double ReCalcArea();
        void delVertex(int index);
        void Reverse();                             //  set vertices order in order to have counter clockwise ordering
        void setStartPointMinY();                   //  Set vertices order such as first point has min Y
        Polygon *enlarge(double Diff);
        Point getCentroid();
        inline Point *GetVertex(int idx) { return(&_Vertices[idx]); }
        inline vector<Point> GetVertices() { return _Vertices;}
        inline Segment GetSegment(int idx) { return(_Segments[idx]);}
        inline int CheckClosed() { if (_Vertices[0] == _Vertices[nVertices-1]) poly_closed = 1; else poly_closed = 0; return poly_closed; }
        Polygon *ConvexHull();
        Polygon *Translate(const Point *pT, int idx_vertex);
        Polygon *Translate(const Point *pT);
        Polygon *Rotate(double angle);
        void ComputeAngles();
        Polygon *Transform(double angle, const Point *pT);
        bool Intersect(Polygon *Poly2);
        inline double getRotation() {return angle_rotation;}
        inline Point getTranslation() { return pt_Translation;}
        inline double getSlope(int idx) { return _Angles[idx]; }
        void align_vertices();
        void BreakLongerEdges(double max_length, double Diff);
        inline double getAngle(int idx)
        {
            double a;
            if ( idx == 0 )
                a = _Angles[0] - _Angles[nVertices-2];
            else
                a = _Angles[idx] - _Angles[idx-1];
            if ( a < 0 ) a += 2*M_PI;           //  Return value between 0 and 2*pi
            return a;
        }
        void Simplify(double max_error);
        void CalcSegments();
        int TypePoly;

        /*
		Overloaded ostream << operator to check for print object of class point to STDOUT
        */
        friend ostream& operator<<(ostream& output,const Polygon& Poly)
        {
            output<<"Polygon with " << Poly.nVertices << " vertices : " ;
            output << Poly._Vertices[0];
            for ( int i = 1; i < Poly.nVertices; i++ )
                output << ", " << Poly._Vertices[i];
            output << '\n';
            return output;
        }


    protected:
        std::vector<Point> _Vertices;
        std::vector<Segment> _Segments;
        std::vector<double> _Angles;
        Rectangle BoundingBox;
        double poly_area;
        Point poly_centroid;
        int poly_centroid_ok;
        int poly_clockwise;
        int poly_closed;
        int poly_changed;
        double precision;
        double angle_rotation;
        Point pt_Translation;

    private:
};

extern uint64_t nbEmptyPoly;
extern uint64_t nbVectorPoly;
extern uint64_t nbRotatedPoly;
extern uint64_t nbTranslatedPoly;
extern uint64_t nbTransformedPoly;
extern uint64_t nbCachedPoly;
extern uint64_t nbCreatedEmptyPoly;
extern uint64_t nbCreatedVectorPoly;
extern uint64_t nbCreatedRotatedPoly;
extern uint64_t nbCreatedTranslatedPoly;
extern uint64_t nbCreatedTransformedPoly;
extern uint64_t nbCreatedCachedPoly;


#endif // GEOMETRY_H
