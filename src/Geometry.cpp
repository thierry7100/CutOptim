#include "Geometry.h"

#define DEBUG_CONSOLE   0

uint64_t nbEmptyPoly = 0;
uint64_t nbVectorPoly = 0;
uint64_t nbRotatedPoly = 0;
uint64_t nbTranslatedPoly = 0;
uint64_t nbTransformedPoly = 0;
uint64_t nbCachedPoly = 0;
uint64_t nbCreatedEmptyPoly = 0;
uint64_t nbCreatedVectorPoly = 0;
uint64_t nbCreatedRotatedPoly = 0;
uint64_t nbCreatedTranslatedPoly = 0;
uint64_t nbCreatedTransformedPoly = 0;
uint64_t nbCreatedCachedPoly = 0;

Polygon::Polygon()
{
    //constructor, build an empty vector
    //  Bounding box is also empty
    BoundingBox.A = Point(0,0);
    BoundingBox.B = Point(0,0);
    nVertices = 0;
    poly_centroid_ok = 0;           //  Not yet computed;
    poly_area = 0.0;
    poly_clockwise = 0;
    poly_closed = 0;
    poly_changed = 1;
    angle_rotation = 0;
    precision = my_Default_Precision;
    TypePoly = EmptyPoly;
    nbEmptyPoly++;
    nbCreatedEmptyPoly++;
}

Polygon::Polygon(vector<Point> &points)
{
int n = points.size();

    BoundingBox.A = Point(0,0);
    BoundingBox.B = Point(0,0);
    nVertices = 0;
    poly_centroid_ok = 0;           //  Not yet computed;
    poly_area = 0.0;
    poly_clockwise = 0;
    poly_closed = 0;
    poly_changed = 1;
    precision = my_Default_Precision;
    for (int i = 0; i < n; i++ )
    {
        addVertice(points[i]);      //  Add this point, but with the desired precision
    }
    //  Add last point which is first
    addVertice(points[0]);
    TypePoly = VectorPoly;
    nbVectorPoly++;
    nbCreatedVectorPoly++;
}

Polygon::~Polygon()
{
    if ( TypePoly == EmptyPoly )
    {
        nbEmptyPoly--;
    }
    else if (TypePoly == RotatedPoly )
    {
        nbEmptyPoly--;
        nbRotatedPoly--;
    }
    else if (TypePoly == TranslatedPoly )
    {
        nbEmptyPoly--;
        nbTranslatedPoly--;
    }
    else if (TypePoly == TransformedPoly )
    {
        nbEmptyPoly--;
        nbTransformedPoly--;
    }
    else if (TypePoly == CachedPoly )
    {
        nbEmptyPoly--;
        nbCachedPoly--;
    }
    else
    {
        nbVectorPoly--;
    }
}

void Polygon::setPrecision(double epsilon)
{
    precision = epsilon;
}

void Polygon::addVertice(Point &p)
{
    Point lP = Point(round(p.x/precision)*precision, round(p.y/precision)*precision);
    if ( nVertices > 0 && lP == _Vertices[nVertices-1] ) return;     //  Do NOT Add point if this is the same as the previous one.
    _Vertices.push_back(lP);
    nVertices = _Vertices.size();
    if ( nVertices == 1 )
    {
        BoundingBox.A = lP;
        BoundingBox.B = lP;
    }
    else
    {
        if ( lP.x < BoundingBox.A.x ) BoundingBox.A.x = lP.x;       //  To avoid rounding errors
        if ( lP.y < BoundingBox.A.y ) BoundingBox.A.y = lP.y;       //  To avoid rounding errors
        if ( lP.x > BoundingBox.B.x ) BoundingBox.B.x = lP.x;       //  To avoid rounding errors
        if ( lP.y > BoundingBox.B.y ) BoundingBox.B.y = lP.y;       //  To avoid rounding errors
    }
    if ( nVertices > 2)
        poly_closed = _Vertices[0] == _Vertices[nVertices-1];
    else
        poly_closed = 0;
    poly_changed = 1;
    poly_area = area();
}

double Polygon::area()
{
double A = 0.0;

    if ( nVertices < 3) return 0.0;
    if ( poly_closed == 0 ) return(0.0);
    if ( poly_changed == 0 && poly_area > 0 ) return(poly_area);
    for (int i = 0; i < nVertices - 1; i++)
    {
        A += _Vertices[i].x*_Vertices[i+1].y - _Vertices[i+1].x*_Vertices[i].y;
    }
    if ( A < 0 )
    {
        poly_clockwise = 1;
        A /= -2.0;
    }
    else
    {
        poly_clockwise = -1;
        A /= 2.0;
    }
    poly_area = A;
    poly_changed = 0;
	return poly_area;
}

double Polygon::distance(Point &p)
{
Point lP = Point(round(p.x/precision)*precision, round(p.y/precision)*precision);
double d1;
double dmin = 1e20;

    for ( int i = 0; i < nVertices-1; i++)
    {
        //  Compute distance from point lP to segment vertice[i], vertice[i+1]
        //  If dot product AB.AP is < 0 return distance(A, P)
        //  Indeed, if this is the case P is outer segment further than A
        if ( (_Vertices[i+1].x - _Vertices[i].x)*(lP.x - _Vertices[i].x) + (_Vertices[i+1].y-_Vertices[i].y)*(lP.y-_Vertices[i].y) < 0 )
            d1 = ( (lP.x - _Vertices[i].x)*(lP.x- _Vertices[i].x) + (lP.y-_Vertices[i].y)*(lP.y-_Vertices[i].y));
        //  If dot product BA.BP si > 0 return distance(B, P)
        else if ( (_Vertices[i].x - _Vertices[i+1].x)*(lP.x - _Vertices[i+1].x) + (_Vertices[i].y-_Vertices[i+1].y)*(lP.y-_Vertices[i+1].y) < 0 )
            d1 = ( (lP.x - _Vertices[i+1].x)*(lP.x- _Vertices[i+1].x) + (lP.y-_Vertices[i+1].y)*(lP.y-_Vertices[i+1].y));
        else
        //  Distance point to line is a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b
        //  line coefficients : a = yA - yB  b = xB - xA c = A.x * B.y - B.x * A.y
        {
            double a = _Vertices[i].y - _Vertices[i+1].y;
            double b = _Vertices[i+1].x - _Vertices[i].x;
            double c = _Vertices[i].x * _Vertices[i+1].y- _Vertices[i+1].x * _Vertices[i].y ;
            d1 = (a* lP.x + b * lP.y + c)*(a * lP.x + b * lP.y + c)/(a*a + b*b);
        }
        dmin = fmin(d1, dmin);
    }
    return(dmin);
}

//  Change the order of vertices in order to be counter clockwise

void Polygon::Reverse(void)
{
    if ( poly_clockwise < 0) return;    //  Oups OK before entering
    if ( poly_clockwise == 0 ) return;  //  Oups not known, unable to change
    for ( int i = 0, j = nVertices-1; i < j ; i++, j--)
    {
        Point p = _Vertices[i];
        _Vertices[i] = _Vertices[j];
        _Vertices[j] = p;
    }
    poly_clockwise = -1;
}

//  Set polygon vertices such as the minimum Y is reached at offset 0 (and nVertice-1)

void Polygon::setStartPointMinY()
{
    double ymin = 1e30;
    int idxmin = 0;
    for ( int i = 0; i < nVertices-1; i++)
    {
        if ( _Vertices[i].y < ymin )
        {
            ymin = _Vertices[i].y;
            idxmin = i;
        }
    }
    if ( idxmin == 0 ) return;      //  Already ordered, return;
    //  Create a new vector
    std::vector<Point>  tmp_vertices;
    for ( int i = 0; i < nVertices - 1; i++)
    {
        tmp_vertices.push_back(_Vertices[(i + idxmin)%(nVertices-1)]);
    }
    tmp_vertices[nVertices-1] = tmp_vertices[0];
    for ( int i = 0; i < nVertices; i++)
    {
        _Vertices[i] = tmp_vertices[i];
    }
}

//  Create a polygon larger than this, which is at a distance min Diff
//  Return a pointer on the created polygon

Polygon *Polygon::enlarge(double Diff)
{
int First = 1;
Point OldA2 = Point(0, 0);
Point OldB2 = Point(0, 0);
Segment A2B2 = Segment(OldA2, OldB2);
Line OldLine(1, 1, 0);
Line FirstLine(1, 1, 0);
int newVertices = 0;
double a, b, c;

    //  Then create larger polygon.
    Polygon *LargePoly = new Polygon();

    for ( int i = 0; i < nVertices - 1; i++)
    {
        //  For each segment create a segment which is parallel to the original segment, but translated at distance Diff
        a = _Vertices[i].y - _Vertices[i+1].y;
        b = _Vertices[i+1].x - _Vertices[i].x;
        c = _Vertices[i].x * _Vertices[i+1].y - _Vertices[i+1].x * _Vertices[i].y;
        //  If equation of first segment is ax + by + c = 0, the translated segment has an equation like ax + by + C = 0 where C = c + Diff*sqrt(a*a + b*b)
        //  Here we will first build a segment which begins in A' and end in B' where A' is the point on the line perpendicular to the segment and which starts in A
        double C = c + Diff*sqrt(a*a+b*b);
        Line A1B1 = Line(a, b, C);
        //  Now the lines which are perpendicular in A or B to the segment
        Line AB_perpendicular_A = Line(-b, a, b*_Vertices[i].x - a*_Vertices[i].y);
        Line AB_perpendicular_B = Line(-b, a, b*_Vertices[i+1].x - a*_Vertices[i+1].y);
        // Compute intersection of AB_perpendicular_A and A1B1 (point A1)
        Point A1 = A1B1.Intersect(AB_perpendicular_A);
        Point B1 = A1B1.Intersect(AB_perpendicular_B);
        //  Now enlarge segment to reach A2 which is Diff further from A1 and B2 which is also Diff further from B1
        double distA1B1 = A1.Distance(B1);
        double Coeff = Diff/distA1B1;
        //  Now compute A2 and B2
        //  A2B1 = Coeff*A1B1
        Point A2 = A1 + Coeff*(A1-B1);
        Point B2 = B1 + Coeff*(B1-A1);
        A2B2 = Segment(A2, B2);
        //  Compute intersection between segment A2B2 and Old Segment A2B2 if not processing the first vertice
        if ( !First )
        {
            Point A3 = A1B1.Intersect(OldLine);
            //  If A3 is between A2 and B2, keep segment A3 - B2
            if (A2B2.InSegment(A3) )
            {
                //Change end of previous segment to A3, so change vertice of largePoly
                LargePoly->_Vertices[newVertices-1] = A3;
                //  And add new point B2
                LargePoly->addVertice(B2);
                newVertices = LargePoly->_Vertices.size();
#if DEBUG_CONSOLE > 1
                printf("Changing last point to A3 (%.3f,%.3f) and add one point (%d) B2 : (%.3f,%.3f)\n", A3.x, A3.y, newVertices, B2.x, B2.y);
#endif // DEBUG_CONSOLE
            }
            else
            {
                // Add vertice A2 and B2
                LargePoly->addVertice(A2);
                LargePoly->addVertice(B2);
                newVertices = LargePoly->_Vertices.size();
#if DEBUG_CONSOLE > 1
                printf("Add 2 points (%d and %d) A2 point (%.3f,%.3f) and B2 (%.3f,%.3f)\n", newVertices-1, newVertices, A2.x, A2.y, B2.x, B2.y);
#endif // DEBUG_CONSOLE
            }
        }
        else
        {
            LargePoly->addVertice(A2);
            LargePoly->addVertice(B2);
            newVertices = LargePoly->_Vertices.size();
            FirstLine = A1B1;
            First = 0;
#if DEBUG_CONSOLE > 1
            printf("Add First %d points : (%.3f, %.3f)   (%.3f,%.3f)\n", newVertices, A2.x, A2.y, B2.x, B2.y);
#endif // DEBUG_CONSOLE
        }
        OldA2 = A2;
        OldB2 = B2;
        OldLine = A1B1;
    }
    //  Now change first segment, compute intersection with OldLine
    Point A3;
    if ( OldLine.hasIntersect(&FirstLine, A3 ) )
    {
        //  If A3 is between A2 and B2, keep segment A3 - B2
        if (A2B2.InSegment(A3) )
        {
            //Change end of previous segment to A3, so change vertice of largePoly
            LargePoly->_Vertices[newVertices-1] = A3;
            //  And set First point of Poly to be the same
            LargePoly->_Vertices[0] = A3;
#if DEBUG_CONSOLE > 1
            printf("Changing first/last point (%d) A3 : (%.3f,%.3f)\n", newVertices, A3.x, A3.y);
    #endif // DEBUG_CONSOLE
        }
        else
        {
            // Add vertice between B2 and first point
            LargePoly->addVertice(LargePoly->_Vertices[0]);
            newVertices = LargePoly->_Vertices.size();
    #if DEBUG_CONSOLE > 1
            printf("Add 1 point (%d) (%3.f,%.3f) , closing polygon\n", newVertices, LargePoly->_Vertices[0].x, LargePoly->_Vertices[0].y);
    #endif // DEBUG_CONSOLE
        }
    }
    else
    {       //  No intersection, add point
            // Add vertice between B2 and first point
        LargePoly->addVertice(LargePoly->_Vertices[0]);
        newVertices = LargePoly->_Vertices.size();
#if DEBUG_CONSOLE > 1
        printf("NoIntersect, Add 1 point (%d) (%3.f,%.3f) , closing polygon\n", newVertices, LargePoly->_Vertices[0].x, LargePoly->_Vertices[0].y);
#endif // DEBUG_CONSOLE
    }
    LargePoly->CheckClosed();
    LargePoly->ReCalcBoundingBox();

    //  Now check that all edges of the enlarged polygon dont intersect
    //  If edge I (_Vertice[i] --> _Vertice[i+1]) intersect edge j (_Vertice[j] --> _Vertice[j+1] ) in point Ai with j > i
    //  1) Change _Vertice[i+1] to Ai
    //  2) Delete All vertices from i+2 to j
    //  The resulting polygon will have an edge from _Vertice[i] --> Ai and Ai --> Vertice[j] (with old indexes i and j).
#if DEBUG_CONSOLE > 0
    printf("Large Poly before post processing with %d (%zu) vertices\n", newVertices, LargePoly->_Vertices.size());
    for (int v = 0; v < newVertices; v++)
    {
        printf("Vertex %d : (%.3f,%.3f)\n", v, LargePoly->_Vertices[v].x, LargePoly->_Vertices[v].y);
    }
#endif
    //  After this step, remove vertices which are inside the polygon.
    //  The edges should never intersect
    int modified = true;
    Point Intersection;
    while ( modified )
    {
        modified = false;
        for (int i = 0; i < newVertices - 2; i++)
        {
            Segment AiBi = Segment(LargePoly->_Vertices[i], LargePoly->_Vertices[i+1]);
            for ( int j = i+2; j < newVertices - 1; j++)
            {
                if ( j == newVertices-2 && i == 0 ) break;      //  Do not check this case, last segment will always intersect first one !
                Segment AjBj = Segment(LargePoly->_Vertices[j],LargePoly->_Vertices[j+1]);
                if ( AiBi.isCrossing(&AjBj, &Intersection))
                {
                //  The two segments are crossing
                //  First change Bi to Intersection
#if DEBUG_CONSOLE > 0
                    printf("Intersection detected between segment %d (%.3f,%.3f -- %.3f,%.3f) and segment %d (%.3f,%.3f -- %.3f,%.3f)\n",
                                i, LargePoly->_Vertices[i].x, LargePoly->_Vertices[i].y, LargePoly->_Vertices[i+1].x, LargePoly->_Vertices[i+1].y,
                                j, LargePoly->_Vertices[j].x, LargePoly->_Vertices[j].y, LargePoly->_Vertices[j+1].x, LargePoly->_Vertices[j+1].y);
                    printf("Changing end of segment to intersection point (%.3f,%.3f)\n", Intersection.x, Intersection.y);
#endif
                //  Build 2 polygons, first one contains segments between i+1 (modified to intersection) and j also modified as intersection
                //  The second one contains all other points
                    Polygon *Poly1 = new Polygon;
                    Polygon *Poly2 = new Polygon;
                    for ( int k = 0; k < newVertices; k++)
                    {
                        if ( k < i || k > j )
                        {
                            Poly1->addVertice(LargePoly->_Vertices[k]);
                        }
                        else if ( k == i )
                        {
                            Poly1->addVertice(LargePoly->_Vertices[k]);
                            Poly1->addVertice(Intersection);
                        }
                        else if ( k == i+1 )
                        {
                            Poly2->addVertice(Intersection);
                            Poly2->addVertice(LargePoly->_Vertices[k]);
                        }
                        else if ( k < j )
                        {
                            Poly2->addVertice(LargePoly->_Vertices[k]);
                        }
                        else
                        {
                            Poly2->addVertice(LargePoly->_Vertices[k]);
                            Poly2->addVertice(Intersection);
                        }
                    }
                    //  Keep the largest one
                    double area1 = fabs(Poly1->area());
                    double area2 = fabs(Poly2->area());
                    if ( area1 > area2 )
                    {
#if DEBUG_CONSOLE > 0
                        printf("Area1 = %.3f, Area2 = %.3f, keep polygon1", area1, area2);
#endif
                        delete Poly2;
                        delete LargePoly;
                        LargePoly = Poly1;
                    }
                    else
                    {
#if DEBUG_CONSOLE > 0
                        printf("Area1 = %.3f, Area2 = %.3f, keep polygon2", area1, area2);
#endif
                        delete Poly1;
                        delete LargePoly;
                        LargePoly = Poly2;
                    }
                    newVertices = LargePoly->_Vertices.size();
                    modified = true;
                    break;
                }
            }
            if ( modified ) break;
        }

    }
    for (int i = 0; i < newVertices; i++)
    {
        if ( isInPoly(&LargePoly->_Vertices[i]) > 0 || distance(LargePoly->_Vertices[i]) < Diff - precision )
        {
            // Inside poly or too close, delete vertex
            LargePoly->delVertex(i);
            newVertices = LargePoly->nVertices;
        }
    }

    LargePoly->CheckClosed();
    LargePoly->ReCalcArea();
#if DEBUG_CONSOLE > 0
    Rectangle BBox = LargePoly->GetBoundingBox();
    printf("After postprocessing, %d vertices, area = %.3f, Bounding box (%.3f,%.3f) - (%.3f,%.3f), closed = %d\n",
        LargePoly->nVertices, LargePoly->area(), BBox.A.x, BBox.A.y, BBox.B.x, BBox.B.y, LargePoly->isClosed());
#endif
    return(LargePoly);
}

//  Force the computation of polygon area.
//  Useful if _Vertices was modified
//  Return the polygon area

double Polygon::ReCalcArea()
{
    poly_changed = 1;
    return(area());
}

//  Compute the bounding box of the polygon

void Polygon::ReCalcBoundingBox()
{
    BoundingBox.A = _Vertices[0];
    BoundingBox.B = _Vertices[0];
    for ( int i = 1; i < nVertices - 1; i++)
    {
        if ( _Vertices[i].x < BoundingBox.A.x ) BoundingBox.A.x = _Vertices[i].x ;
        if ( _Vertices[i].y < BoundingBox.A.y ) BoundingBox.A.y = _Vertices[i].y ;
        if ( _Vertices[i].x > BoundingBox.B.x ) BoundingBox.B.x = _Vertices[i].x ;
        if ( _Vertices[i].y > BoundingBox.B.y ) BoundingBox.B.y = _Vertices[i].y ;
    }
}

//  Return a rectangle which is the bounding box of the polygon

Rectangle Polygon::GetBoundingBox()
{
Rectangle Ret;

    Ret.A.x = BoundingBox.A.x - precision/3;      //    Avoid rounding errors, enlarge slightly bounding box
    Ret.A.y = BoundingBox.A.y - precision/3;      //    Avoid rounding errors, enlarge slightly bounding box
    Ret.B.x = BoundingBox.B.x + precision/3;      //    Avoid rounding errors, enlarge slightly bounding box
    Ret.B.y = BoundingBox.B.y + precision/3;      //    Avoid rounding errors, enlarge slightly bounding box
    return Ret;
}

//  Return TRUE if x is betwenn d1 and d2
//  Indeed,
static inline bool isBetween(double x, double d1, double d2)
{
    if ( d1 <= d2 )
    {
        return (x >= d1 - my_Default_Precision && x <= d2 + my_Default_Precision);
    }
    return (x >= d2 - my_Default_Precision && x <= d1 + my_Default_Precision);
}

//  Return True if the point is inside the polygon and false otherwise
//  If the point is a vertex of the polygon, return a negative value (- index vertex - 1)
//  If the point is on an edge of the polygon, return 0 but set the segment number
//  Use a ray tracing algorithm for this, launch a ray and if this ray cross the polygon edges an odd number of times the point is inside.
//  The bounding box of the polygon should be OK before entering this function.
//  If the _Vertices points were modified without using addVertice, the call of ReCalcBoundingBox is mandatory.

int Polygon::isInPoly(const Point *p, int *SegmentIdx)
{
double lastX, lastY;
int inside = false;

    if ( SegmentIdx != nullptr ) *SegmentIdx = -1;      //  If provided, initialize to - 1 (not on an edge
    //  Round coordinates to avoid errors
    Point lP = Point(round(p->x/precision)*precision, round(p->y/precision)*precision);

    if ( lP.x < BoundingBox.A.x - precision/3 ) return false;       //  Outside bounding box, impossible, sub precision to avoid rounding errors
    if ( lP.y < BoundingBox.A.y - precision/3) return false;       //  Outside bounding box, impossible
    if ( lP.x > BoundingBox.B.x + precision/3) return false;       //  Outside bounding box, impossible
    if ( lP.y > BoundingBox.B.y + precision/3) return false;       //  Outside bounding box, impossible
    //  If the point is a vertex, return false.
    //  Also checks if the point is on one edge, in this case return false
    lastX = _Vertices[0].x;
    lastY = _Vertices[0].y;
    if ( fabs(lastX-lP.x) + fabs(lastY-lP.y) < precision ) return -1;    //  On vertex 0
    for (int i = 1; i < nVertices; i++)
    {
        if ( fabs(_Vertices[i].x-lP.x) < 50*precision && fabs(_Vertices[i].y-lP.y) < 50*precision) return -(i+1);      //  this is a vertex, return negative value
        //  Check if lP is on edge between LastX,LastY and _Vertice[i]
        if ( isBetween(lP.x, lastX, _Vertices[i].x) && isBetween(lP.y, lastY, _Vertices[i].y) )
        {
            //  Could be on edge, check if points are aligned : compute determinant
            double deltaX =_Vertices[i].x - lastX;
            double deltaY =_Vertices[i].y - lastY;
            double det = (lP.x - lastX)*deltaY - (lP.y - lastY)*deltaX;
            //  With infinite precision, the determinant is null when the points are aligned
            //  Here, we don't have infinite precision (rounding)
            //  It could be quite complex to determine the alignment
            //  Say if fabs(det) is lower than precision*fabs(deltaX + deltaY)
            double max_error_det = 2*(fabs(deltaX) + fabs(deltaY));
#define OUTPUT_DETERMINANT 0

#if OUTPUT_DETERMINANT > 0
            if ( fabs(det) < 10 )
                cout << "det =" << det << ", deltaX=" << deltaX << ", deltaY=" << deltaY << ", computed precision=" << max_error_det*precision<< '\n';
#endif // OUTPUT_DETERMINANT
            if ( fabs(det) < max_error_det * precision )
            {
                if ( SegmentIdx != nullptr ) *SegmentIdx = i - 1; //    Between i - 1 and i.
                return false;      //  On edge, not inside !
            }
        }
        lastX = _Vertices[i].x;
        lastY = _Vertices[i].y;
    }

    lastX = _Vertices[0].x;
    lastY = _Vertices[0].y;
    for (int i = 1; i < nVertices; i++)
    {
        double curX = _Vertices[i].x;
        double curY = _Vertices[i].y;
        //  Note the <= comparison with last Y, to avoid double counting when line intersect edge at the vertex.
        if (!(lP.y <= lastY && lP.y < curY))           //   Not below segment so intersect is possible
        {
            if (!(lP.y >= lastY && lP.y > curY))       //  not Above current segment so intersect is possible
            {
                if ( !(lP.x > lastX && lP.x > curX))  //  Not after current segment, so intersect is possible
                {
                    if (curY == lP.y)                   //  Intersect at vertex ?
                    {
                        double DiffY;
                        //  In this case, change inside if H line really cross the edge
                        //  If DiffY is not null (not horizontal) and DiffY has not the same sign as curY-LastY, this is not a true intersection.
                        if ( i < nVertices - 1)
                            DiffY =  _Vertices[i+1].y - curY;
                        else
                            DiffY = _Vertices[0].y - curY;
                        //  If DiffY is null, do the same thing with next point
                        if ( DiffY == 0 )
                        {
                            if ( i < nVertices - 2)
                                DiffY =  _Vertices[i+2].y - curY;
                            else if ( i == nVertices - 2)
                                DiffY = _Vertices[0].y - curY;
                            else
                                DiffY = _Vertices[1].y - curY;
                        }
                        //  In this case intersect only if curY-LastY has the same sign as DiffY
                        double sign = -1;
                        sign = (curY - lastY) * DiffY;
                        if ( sign > 0 && lP.x <= curX )
                            inside = !inside;
                    }
                    else if ( curY != lastY)         //  Only if segment is not horizontal
                    {
                        double intersect = lastX + (lP.y - lastY) * (curX - lastX) / (curY - lastY);
                        if ( lP.x <= intersect )
                            inside = ! inside;
                    }
                    else
                    {
                        inside = !inside;
                    }
                }
            }
        }
        lastX = curX;
        lastY = curY;
    }
    return inside;
}

//  Delete vertex by index in the polygon
//  do NOT check validity of removing this vertex.

void Polygon::delVertex(int index)
{
    if ( nVertices > 3 )
    {
        _Vertices.erase(_Vertices.begin()+index);
        nVertices = _Vertices.size();
        poly_changed = 1;
        poly_centroid_ok = 0;
    }
}

//  Compute the convex hull of the polygon and return it as a new polygon
extern vector<Point> makeConvexHull(const vector<Point> &points) ;

Polygon * Polygon::ConvexHull()
{
Polygon *Hull = new Polygon;

    Hull->_Vertices = makeConvexHull(_Vertices);
    Hull->nVertices = Hull->_Vertices.size();
    Hull->CheckClosed();
    if ( Hull->poly_closed == 0 )
    {
        //  Not closed, close it
        Hull->addVertice(Hull->_Vertices[0]);
    }
    Hull->ReCalcBoundingBox();
    Hull->ReCalcArea();
    return Hull;
}

//  Return the centroid of the Polygon
//  The formula is described here : https://en.wikipedia.org/wiki/Centroid
//  Use a cache to compute this centroid once.

Point Polygon::getCentroid()
{
double p_area;
double cx = 0, cy = 0;

    if ( poly_changed ) poly_centroid_ok = 0;
    if ( poly_centroid_ok == 0 )    //  Not yet computed, do it
    {
        p_area = area();
        for ( int i = 0; i < nVertices - 1; i++)
        {
            double d = _Vertices[i].x*_Vertices[i+1].y - _Vertices[i+1].x*_Vertices[i].y;
            cx += (_Vertices[i].x + _Vertices[i+1].x)*d;
            cy += (_Vertices[i].y + _Vertices[i+1].y)*d;
        }
        cx /= 6*p_area;
        cy /= 6*p_area;
        poly_centroid = Point(cx, cy);
        poly_centroid_ok = 1;
    }
    return poly_centroid;
}

/**
*  Translate the polygon such as centroid will be moved to pT
*  Centroid must be ok before call
*  Return a new Polygon
*/

Polygon *Polygon::Translate(const Point *pT)
{
Polygon *Translated = new Polygon;

    nbTranslatedPoly++;
    nbCreatedTranslatedPoly++;
    Translated->TypePoly = TranslatedPoly;
    if ( poly_centroid_ok == 0 ) getCentroid();
    Translated->pt_Translation = *pT - poly_centroid;
    Translated->_Vertices = _Vertices;          //      Copy vertice vector
    Translated->nVertices = Translated->_Vertices.size();
    for ( int i = 0; i < nVertices; i++)
    {
        Translated->_Vertices[i].x += Translated->pt_Translation.x;
        Translated->_Vertices[i].y += Translated->pt_Translation.y;
    }
    Translated->_Angles = _Angles;
    Translated->poly_area = poly_area;
    Translated->poly_centroid = poly_centroid + Translated->pt_Translation;
    Translated->angle_rotation = angle_rotation;
    Translated->poly_centroid_ok = 1;
    Translated->poly_changed = 0;
    Translated->poly_closed = poly_closed;
    Translated->poly_clockwise = poly_clockwise;
    Translated->BoundingBox.A = BoundingBox.A + Translated->pt_Translation;
    Translated->BoundingBox.B = BoundingBox.B + Translated->pt_Translation;
    Translated->precision = precision;
    return Translated;
}

//  Translate the polygon such as vertex I is placed on point pT
//  Return a new Polygon

Polygon *Polygon::Translate(const Point *pT, int idx_vertex)
{
Polygon *Translated = new Polygon;

    nbTranslatedPoly++;
    nbCreatedTranslatedPoly++;
    Translated->TypePoly = TranslatedPoly;
    Translated->pt_Translation = *pT - _Vertices[idx_vertex];
    Translated->_Vertices = _Vertices;          //      Copy vertice vector
    Translated->nVertices = Translated->_Vertices.size();
    for ( int i = 0; i < nVertices; i++)
    {
        double newx = Translated->_Vertices[i].x + Translated->pt_Translation.x;
        newx = round(newx/precision)*precision;
        Translated->_Vertices[i].x = newx;
        double newy = Translated->_Vertices[i].y + Translated->pt_Translation.y;
        newy = round(newy/precision)*precision;
        Translated->_Vertices[i].y = newy;
    }
    Translated->_Angles = _Angles;
    Translated->poly_area = poly_area;
    Translated->poly_centroid = poly_centroid + Translated->pt_Translation;
    Translated->angle_rotation = angle_rotation;
    Translated->poly_centroid_ok = 1;
    Translated->poly_changed = 0;
    Translated->poly_closed = poly_closed;
    Translated->poly_clockwise = poly_clockwise;
    Translated->BoundingBox.A = BoundingBox.A + Translated->pt_Translation;
    Translated->BoundingBox.B = BoundingBox.B + Translated->pt_Translation;
    Translated->precision = precision;
    return Translated;
}

//  Rotate a polygon by the angle
//  Rotate from the centroid of the polygon, so it will not be modified
//  Return a new polygon

Polygon *Polygon::Rotate(double angle)
{
Polygon *Rotated = new Polygon;
double cos_a = cos(angle);
double sin_a = sin(angle);
double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

    nbRotatedPoly++;
    nbCreatedRotatedPoly++;
    Rotated->TypePoly = RotatedPoly;
    if ( poly_centroid_ok == 0 ) getCentroid();
    Rotated->precision = precision;
    Rotated->pt_Translation = pt_Translation;
    Rotated->_Vertices = _Vertices;          //      Copy vertice vector
    Rotated->nVertices = Rotated->_Vertices.size();

    for ( int i = 0; i < nVertices; i++)
    {
        double newx = (_Vertices[i].x - poly_centroid.x) * cos_a - (_Vertices[i].y - poly_centroid.y)* sin_a + poly_centroid.x;
        newx = round(newx/precision)*precision;
        Rotated->_Vertices[i].x = newx;
        double newy = (_Vertices[i].x - poly_centroid.x) * sin_a + (_Vertices[i].y - poly_centroid.y)* cos_a + poly_centroid.y;
        newy = round(newy/precision)*precision;
        Rotated->_Vertices[i].y = newy;
        if ( i == 0 )
        {
            xmin = newx;
            xmax = newx;
            ymin = newy;
            ymax = newy;
        }
        else
        {
            if ( newx < xmin ) xmin = newx;
            if ( newy < ymin ) ymin = newy;
            if ( newx > xmax ) xmax = newx;
            if ( newy > ymax ) ymax = newy;
        }
    }
    Rotated->BoundingBox.A = Point(xmin, ymin);
    Rotated->BoundingBox.B = Point(xmax, ymax);
    Rotated->poly_area = poly_area;
    Rotated->poly_centroid = poly_centroid;
    Rotated->angle_rotation = angle;
    Rotated->poly_centroid_ok = 1;
    Rotated->poly_changed = 0;
    Rotated->poly_closed = poly_closed;
    Rotated->poly_clockwise = poly_clockwise;
    return Rotated;
}

//  Transform source polygon : rotate it around its centroid by angle then translate it by pT
//  Centroid must be ok before call
//  Return a new Polygon

Polygon *Polygon::Transform(double angle, const Point *pT)
{
Polygon *Transformed = new Polygon;
double cos_a = cos(angle);
double sin_a = sin(angle);
double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

    nbTransformedPoly++;
    nbCreatedTransformedPoly++;
    Transformed->TypePoly = TransformedPoly;
    if ( poly_centroid_ok == 0 ) getCentroid();
    Transformed->angle_rotation = angle;
    Transformed->pt_Translation = *pT;
    Transformed->precision = precision;
    Transformed->_Vertices = _Vertices;          //      Copy vertice vector
    Transformed->nVertices = Transformed->_Vertices.size();
    Transformed->poly_centroid = poly_centroid + Transformed->pt_Translation;
    for ( int i = 0; i < nVertices; i++)
    {
        double newx = (_Vertices[i].x - poly_centroid.x) * cos_a - (_Vertices[i].y - poly_centroid.y)* sin_a + Transformed->poly_centroid.x;
        newx = round(newx/precision)*precision;
        Transformed->_Vertices[i].x = newx;
        double newy = (_Vertices[i].x - poly_centroid.x) * sin_a + (_Vertices[i].y - poly_centroid.y)* cos_a + Transformed->poly_centroid.y;
        newy = round(newy/precision)*precision;
        Transformed->_Vertices[i].y = newy;
        if ( i == 0 )
        {
            xmin = newx;
            xmax = newx;
            ymin = newy;
            ymax = newy;
        }
        else
        {
            if ( newx < xmin ) xmin = newx;
            if ( newy < ymin ) ymin = newy;
            if ( newx > xmax ) xmax = newx;
            if ( newy > ymax ) ymax = newy;
        }
    }
    Transformed->BoundingBox.A = Point(xmin, ymin);
    Transformed->BoundingBox.B = Point(xmax, ymax);
    Transformed->poly_area = poly_area;
    Transformed->poly_centroid_ok = 1;
    Transformed->poly_changed = 0;
    Transformed->poly_closed = poly_closed;
    Transformed->poly_clockwise = poly_clockwise;
    return Transformed;
}

/**
*  Compute angles of each polygon edge
*  The result is stored in a vector
*/

void Polygon::ComputeAngles()
{
    _Angles.resize(nVertices, 0.0);
    for ( int i = 0; i < nVertices-1; i++)
    {
        double a = atan2(_Vertices[i+1].y - _Vertices[i].y, _Vertices[i+1].x - _Vertices[i].x );
        //  Beware of rounding errors
//        _Angles[i] = round(a/angle_precision)*angle_precision;
        _Angles[i] = a;
    }
    //  Last value is not defined as Pt n-1 = Pt 0
    _Angles[nVertices-1] = 0.0;
}

/**
*   Return true if one segment of current polygon intersect one segment of Poly2
*   Do do this, checks for each segment of Poly2 if it crosses one segment of this
*   If one segment of Poly2 is included in one segment of this (or the reverse), if segments share same direction return true
*   Indeed two polygons couldn't share a segment in the same direction without crossing or being included.
*/

bool Polygon::Intersect(Polygon *Poly2)
{
Segment *pSeg2;
bool ShareSegment = false;      //  If true, one segment is included in another

    for ( int i = 0; i < Poly2->nVertices - 1; i++)
    {
        pSeg2 = &Poly2->_Segments[i];
        if ( pSeg2->xM <= BoundingBox.A.x ) continue;           //  Segment couldn't cross, right from polygon 1
        if ( pSeg2->xm >= BoundingBox.B.x ) continue;           //  Segment couldn't cross, left from polygon 1
        if ( pSeg2->yM <= BoundingBox.A.y ) continue;           //  Segment couldn't cross, below from polygon 1
        if ( pSeg2->ym >= BoundingBox.B.y ) continue;           //  Segment couldn't cross, above from polygon 1
        //  Possible, so check with all edges of polygon 1
        for ( int j = 0; j < nVertices - 1; j++)
        {
            if ( _Segments[j].isCrossingNoEnd(pSeg2, NULL, &ShareSegment) ) //  Cross return true
                return true;
            //  If lines are parallel checks for inclusion, at least partial
            if ( ShareSegment )
            {
                Point ptA1 = pSeg2->getPointA();
                Point ptB1 = pSeg2->getPointB();
                Point ptA2 = _Segments[j].getPointA();
                Point ptB2 = _Segments[j].getPointB();
                if ( _Segments[j].InSegmentNoEnd(ptA1) )
                {
                    //cout << "Shared point A of " << *pSeg2 << " in " << _Segments[j] << " with same rotation side\n";
                    return true;
                }
                if ( _Segments[j].InSegmentNoEnd(ptB1) )
                {
                    //cout << "Shared point B of " << *pSeg2 << " in " << _Segments[j] << " with same rotation side\n";
                    return true;
                }
                if ( pSeg2->InSegmentNoEnd(ptA2) )
                {
                    //cout << "Shared point A of " << _Segments[j] << " in " <<  pSeg2 << " with same rotation side\n";
                    return true;
                }
                if ( pSeg2->InSegmentNoEnd(ptB2) )
                {
                    //cout << "Shared point B of " << _Segments[j] << " in " <<  pSeg2 << " with same rotation side\n";
                    return true;
                }
            }
        }
    }
    return false;
}

/**
*  Delete all vertices which are not necessary, i.e. when removing then will lead to an error less than max_error
*/

void Polygon::Simplify(double max_error)
{
Point A = _Vertices[nVertices-1];           //  Start at end of polygon, to keep indexes OK
Point B;
Point C;
int hasSimplify = 1;

    while ( hasSimplify )
    {
        hasSimplify = 0;
        if ( nVertices == 3 ) break;    //  Not able to remove one more vertex
        //  Now walk through all path elements to delete vertices if they are aligned
        for ( int i = nVertices-2; i > 0; i--)
        {
            B = _Vertices[i];
            C = _Vertices[i-1];
            Segment AC = Segment(A, C);
            if ( AC.sqrDistancePoint(B) < max_error)
            {
                //  Point B is not necessary, remove
                delVertex(i);
                hasSimplify = 1;
                break;              //  Will start a new loop
            }
            A = B;
        }
        //  Now the case of first vertex (and last one)
        A = _Vertices[1];
        B = _Vertices[0];
        C = _Vertices[nVertices-2];
        Segment AC = Segment(A, C);
        if ( AC.sqrDistancePoint(B) < max_error)
        {
            //  Point B, with index 0 is not neceassry, remove it
            delVertex(0);               // Beware, this is also index nVertices - 1
            _Vertices[nVertices-1] = _Vertices[0];  //  New index 0, was 1 before removing point
            hasSimplify = 1;
        }
    }
    ReCalcArea();
    ReCalcBoundingBox();
}

//  Check if polygon is included in the Big one
//  Return true if all vertices of this polygon are inside Big

int Polygon::Poly_in_Poly(Polygon *Big)
{
    //  First check bounding box, this will save time in most cases.
    if ( BoundingBox.A.x < Big->BoundingBox.A.x ) return false;     //  Impossible, bounding box not included
    if ( BoundingBox.A.y < Big->BoundingBox.A.y ) return false;     //  Impossible, bounding box not included
    if ( BoundingBox.B.x > Big->BoundingBox.B.x ) return false;     //  Impossible, bounding box not included
    if ( BoundingBox.B.y > Big->BoundingBox.B.y ) return false;     //  Impossible, bounding box not included
    //  True algoritm, check all vertices
    for ( int i = 0; i < nVertices; i++)
    {
        if ( ! (Big->isInPoly(&_Vertices[i]) > 0) )
            return false;           //  Vertex i is not in Big, not included !
    }
    return true;
}

void Polygon::CalcSegments()
{
    _Segments.resize(nVertices-1);      //  One segment between 2 vertices...
    for (int i = 0; i < nVertices-1; i++ )
    {
        _Segments[i] = Segment(_Vertices[i], _Vertices[i+1]);
    }
}

//  Examine each edges of the Polygon.
//  If an edge is too long, break it in several segments. Each segment will be roughly inline, but not exactly in line to avoid errors linked to parallel segments

void Polygon::BreakLongerEdges(double max_length, double Diff)
{
std::vector<Point>::iterator it;
double Sign = 1;
int Modified = 0;

    do {
        Modified = 0;
        //  Do it reverse to avoid problem with indexes when inserting new points
        for ( int i = nVertices - 1; i > 0; i-- )
        {
            double l = _Vertices[i].Distance(_Vertices[i-1]);       //  Segment length
            if ( l > max_length )
            {
                //  Too long, break it
                int nSegment = (int)(l / max_length) + 1;       //  Number of segments, will add nSegment - 1 vertices
                //  Compute Delta off line
                //  For each segment create a segment which is parallel to the original segment, but translated at distance Diff
                //  Original segment
                double a = _Vertices[i-1].y - _Vertices[i].y;
                double b = _Vertices[i].x - _Vertices[i-1].x;
                double c = _Vertices[i-1].x * _Vertices[i].y - _Vertices[i].x * _Vertices[i-1].y;
                //  If equation of first segment is ax + by + c = 0, the translated segment has an equation like ax + by + C = 0 where C = c + Diff*sqrt(a*a + b*b)
                //  Here we will first build a segment which begins in A' and end in B' where A' is the point on the line perpendicular to the segment and which starts in A
                double C = c + Diff*sqrt(a*a+b*b);
                Line A1B1 = Line(a, b, C);
                //  Now the lines which are perpendicular in A or B to the segment
                Line AB_perpendicular_A = Line(-b, a, b*_Vertices[i-1].x - a*_Vertices[i-1].y);
                // Compute intersection of AB_perpendicular_A and A1B1 (point A1)
                Point A1 = A1B1.Intersect(AB_perpendicular_A);
                Point Delta = Point(A1.x - _Vertices[i-1].x, A1.y - _Vertices[i-1].y );
                //  Insert only one point at a time to avoid indices mismatch
                Point As = Point(_Vertices[i-1].x + (_Vertices[i].x-_Vertices[i-1].x) / double(nSegment), _Vertices[i-1].y + (_Vertices[i].y-_Vertices[i-1].y) / double(nSegment));
                As.x += Sign * Delta.x;
                As.y += Sign * Delta.y;
                As.x = round(As.x/precision)*precision;
                As.y = round(As.y/precision)*precision;
                Sign *= -1;
                it = _Vertices.begin();
                _Vertices.insert(it + i, As);
                BoundingBox.A.x = fmin(BoundingBox.A.x, As.x);
                BoundingBox.A.y = fmin(BoundingBox.A.y, As.y);
                BoundingBox.B.x = fmax(BoundingBox.B.x, As.x);
                BoundingBox.B.y = fmax(BoundingBox.B.y, As.y);
                Modified = 1;
                poly_changed = 1;
                nVertices = _Vertices.size();
                break;
            }
        }
    } while ( Modified );
    ReCalcArea();
    ReCalcBoundingBox();
}
