#define NANOSVG_IMPLEMENTATION

#include "SvgDoc.h"

#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <string>
#include <math.h>
#include "Geometry.h"
#include "ConvexHull.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

#define SIMPLIFY_POLYGON 1

#define FIXED_BIT       16

#define DEBUG_BEZIER    0


enum ReasonPlacementFailure {
	FAIL_OUTSIDE_SHEET = 1,
	FAIL_VERTEX_IN_POLYGON = 2,
	FAIL_EDGE_INTERSECT = 3
};

//static int DebugPassage = 0;

//  Compute the max distance between the bezier curve and a line segment
//  Actually return the square of the distance, easier to compute
//  The distance is the max between the control points and the segment joining the 2

static double Distance_bezier_segment(double x, double y, NSVGPathElt *Elt)
{
    //  Compute parameters of the line segment between Start and End
    double b = Elt->EndX - x;
    double a = y - Elt->EndY;
    double c = x * Elt->EndY - Elt->EndX*y;
    //  Compute distance between first control point and the line segment Start - End
    // Distance between point and line is (a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b)
    double d1 = a * Elt->BzCtrl1X + b * Elt->BzCtrl1Y + c;
    d1 *= d1;
    d1 /= a*a + b*b;
    //  Compute distance between second control point and the line segment Start - End
    // Distance between point and line is (a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b)
    double d2 = a * Elt->BzCtrl2X + b * Elt->BzCtrl2Y + c;
    d2 *= d2;
    d2 /= a*a + b*b;
    return(fmax(d1, d2));
}

//  Break each element of type bezier cubic path into segments

static void SplitBezierCurves(NSVGpath *path, double flat_factor)
{
NSVGPathElt *Elt = path->PathElts;
float x, y;
float xmin, xmax, ymin, ymax;

    NSVGpath *dPath = (NSVGpath *) malloc(sizeof(NSVGpath));
    path->PolyLine = dPath;
    dPath->StartX = path->StartX;          //  Copy elements
    dPath->StartY = path->StartY;
    dPath->nElts = path->nElts;
    dPath->next = 0;
    dPath->PolyLine = dPath;
    dPath->ClosePolygon = NULL;
    dPath->hasBezier = 0;
    dPath->closed = path->closed;
    //  Copy all elements of path
    NSVGPathElt *dElt, *oldElt = 0;
    while ( Elt )
    {
        dElt = (NSVGPathElt *) malloc(sizeof(NSVGPathElt));
        *dElt = *Elt;
        dElt->next = NULL;
        if ( oldElt )
            oldElt->next = dElt;
        else
            path->PolyLine->PathElts = dElt;
        oldElt = dElt;
        Elt = Elt->next;
    }
    //  Now walk through all path elements to transform bezier curves into line segments
    x = path->StartX;
    y = path->StartY;
    xmin = x;
    xmax = x;
    ymin = y;
    ymax = y;
    Elt = path->PolyLine->PathElts;
    while ( Elt != NULL )
    {
        if ( Elt->Type == NSVG_PATH_BEZIER_C )
        {
            //  Break bezier curves into line segments
            //  For this, split bezier curve into two parts, until parts are nearly line segments
            double d_bezier = Distance_bezier_segment(x, y, Elt);
            if ( d_bezier <=  flat_factor )
            {
                Elt->Type = NSVG_PATH_LINE;     //  nearly a line segment, transform type
            }
            else
            {
                //  split bezier into two parts, then resume processing at the same point.
                NSVGPathElt *Elt2 = (NSVGPathElt *) malloc(sizeof(NSVGPathElt));    //  Create a new element
                //  First create the linked list
                Elt2->next = Elt->next;
                Elt->next = Elt2;
                Elt2->Type = NSVG_PATH_BEZIER_C;
                //  The end point of the second half is the previous end point.
                Elt2->EndX = Elt->EndX;
                Elt2->EndY = Elt->EndY;
                //  Now compute new bezier control points.
                //  First is halfway between start and first control point
                float m1x, m1y;
                m1x = (x + Elt->BzCtrl1X) / 2;
                m1y = (y + Elt->BzCtrl1Y) / 2;
                //  Temp point, halfway between the 2 two control points
                float m2x, m2y;
                m2x = (Elt->BzCtrl2X + Elt->BzCtrl1X) / 2;
                m2y = (Elt->BzCtrl2Y + Elt->BzCtrl1Y) / 2;
                // Then a point halfway between the second control point and the end point
                float m3x, m3y;
                m3x = (Elt->BzCtrl2X + Elt->EndX) / 2;
                m3y = (Elt->BzCtrl2Y + Elt->EndY) / 2;
                // Then compute a new point in the middle of m1, m2
                float m4x, m4y;
                m4x = (m1x + m2x)/2;
                m4y = (m1y + m2y)/2;
                // Then compute a new point in the middle of m2, m3
                float m5x, m5y;
                m5x = (m2x + m3x)/2;
                m5y = (m2y + m3y)/2;
                // And at last compute a new point in the middle of m4, m5
                float m6x, m6y;
                m6x = (m4x + m5x)/2;
                m6y = (m4y + m5y)/2;
                //  The new control points for the first half are m1 and m4, and it goes up to m6
                Elt->EndX = m6x;
                Elt->EndY = m6y;
                Elt->BzCtrl1X = m1x;
                Elt->BzCtrl1Y = m1y;
                Elt->BzCtrl2X = m4x;
                Elt->BzCtrl2Y = m4y;
                //  The new control points for the second half are m5 and m3
                Elt2->BzCtrl1X = m5x;
                Elt2->BzCtrl1Y = m5y;
                Elt2->BzCtrl2X = m3x;
                Elt2->BzCtrl2Y = m3y;
                path->PolyLine->nElts++;
                continue;           //  Restart procrssing of the segment, which is smaller now.
            }

        }
        //  Change start point for next, which is end of this one
        x = Elt->EndX;
        y = Elt->EndY;
        xmin = fmin(xmin, x);
        xmax = fmax(xmax, x);
        ymin = fmin(ymin, y);
        ymax = fmax(ymax, y);
        //  Then next element
        Elt = Elt->next;
    }
    path->PolyLine->bounds[0] = xmin;
    path->PolyLine->bounds[1] = ymin;
    path->PolyLine->bounds[2] = xmax;
    path->PolyLine->bounds[3] = ymax;
}

static void Path2Polygon(Polygon *Poly, NSVGpath *path)
{
struct NSVGPathElt *Elt = path->PathElts;
Point p;

    p = Point(path->StartX, path->StartY);
    Poly->addVertice(p);
    for (int n = 0; n < path->nElts; n++)
    {
        p = Point(Elt->EndX, Elt->EndY);
        Poly->addVertice(p);
        Elt = Elt->next;
    }
}

//  Walk through the polyline and delete vertices if the error is smaller than max_error

void SvgDoc::SimplifyPath(NSVGpath *path, double max_error)
{
NSVGPathElt *Elt, *Elt2, *oldElt = NULL;
float x, y;
float xmin, xmax, ymin, ymax;

    //  Now walk through all path elements to delete vertices if they are aligned
    x = path->StartX;
    y = path->StartY;
    xmin = x;
    xmax = x;
    ymin = y;
    ymax = y;
    Elt = path->PathElts;
    while ( Elt != NULL )
    {
        Elt2 = Elt->next;
        if ( Elt2 != NULL )
        {
            //  Build segment x,y --> Elt2->EndX, EndY
            //  Compute parameters of the line segment between Start and End
            double b = Elt2->EndX - x;
            double a = y - Elt2->EndY;
            double c = x * Elt2->EndY - Elt2->EndX*y;
            //  Compute distance between Elt->EndX, EndY and the segment
            // Distance between point and line is (a * p.x + b * p.y + c)*(a * p.x + b * p.y + c)/(a*a + b*b)
            double d1 = a * Elt->EndX + b * Elt->EndY + c;
            d1 *= d1;
            d1 /= a*a + b*b;
            if ( d1 < max_error )
            {
                //  remove Elt
                if ( oldElt == NULL)
                    path->PathElts = Elt2;
                else
                    oldElt->next = Elt2;
                path->nElts--;
                Elt = Elt2;
                continue;
            }
            else
            {
                x = Elt->EndX;
                y = Elt->EndY;
                xmin = fmin(xmin, x);
                xmax = fmax(xmax, x);
                ymin = fmin(ymin, y);
                ymax = fmax(ymax, y);
                //  Then next element
                oldElt = Elt;
                Elt = Elt2;
                continue;
            }
        }
        break;
    }
    path->PolyLine->bounds[0] = xmin;
    path->PolyLine->bounds[1] = ymin;
    path->PolyLine->bounds[2] = xmax;
    path->PolyLine->bounds[3] = ymax;
}


//  Transform each path into a polyline, if not already the case
//  Then simplify path to save some running time

void SvgDoc::TransformPaths(double flat_factor, bool KeepNested)
{
    int nShape, nPath;

    nShape = 0;
    if ( SvgData == NULL ) return;
    if ( debug_level > 0 )
    {
        OutDebug << "Entering TransformPaths\n" << std::flush;
    }
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        nPath = 0;
        if ( debug_level > 0 )
        {
            OutDebug << "Shape " << nShape << " id:" << shape->id << " Stroke[" << (int) shape->stroke.type << " " << setw(6) << setfill('0') << std::hex << (shape->stroke.color&0xFFFFFF) << std::dec << "] width:" << shape->strokeWidth << "\n" << std::flush;
        }
        //  Transform each path to a polygon
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next, nPath++)
		{
            snprintf(path->id, 99, "%s_%d", shape->id, nPath);
            path->id[99] = 0;
            //  Walk through the path and delete vertices which have trange coordinates (nan)
            //  This has been observed in some cases.
            //  Indeed vertices are very close in these cases, so it is safe to remove vertex with nan coordinates.
            //  However, this can't be the first one.
            int StrangeControlPoints = 0;
            int StrangeEndPoints = 0;
            NSVGPathElt *Elt = path->PathElts;
            NSVGPathElt *PrevElt = NULL;
            for (int i = 0; i < path->nElts; i++)
            {
                int ShouldDelete = 0;
                if ( isnan(Elt->EndX) || isnan(Elt->EndY) )
                {
                    StrangeEndPoints++;
                    ShouldDelete = 1;
                }
                if ( isnan(Elt->BzCtrl1X) || isnan(Elt->BzCtrl1Y) || isnan(Elt->BzCtrl2X) || isnan(Elt->BzCtrl2Y) )
                {
                    StrangeControlPoints++;
                    ShouldDelete = 1;
                }
                if ( ShouldDelete )
                {
                    if ( PrevElt == NULL )
                    {
                        if ( debug_level > 0 )OutDebug << "nan encoutered as first element in path, aborting\n";
                        cerr << "nan encoutered as first element in path, aborting\n";
                        exit(1);
                    }
                    //  Remove vertex
                    path->nElts--;
                    PrevElt->next = Elt->next;
                    Elt = Elt->next;
                    continue;
                }
                PrevElt = Elt;
                Elt = Elt->next;
            }

            if ( debug_level > 2)
            {
                OutDebug << "Before Transform path, Shape " << nShape << " Path " << nPath << " with " << path->nElts << " HasBezier:" << path->hasBezier << "\n";
                OutDebug << "Path has " << StrangeEndPoints << " StrangeEndPoints and "<< StrangeControlPoints << " StrangeControlPoints (removed)\n";
                OutDebug << "Starting point=" << path->StartX << "," << path->StartY << "\n";
                OutDebug << "  Bounds (" << path->bounds[0] << "," << path->bounds[1] << ") --> (" << path->bounds[2] << "," << path->bounds[3] << "\n";
                NSVGPathElt *Elt = path->PathElts;
                for (int i = 0; i < path->nElts; i++)
                {
                    OutDebug << "   Elt " << i << " type " << Elt->Type << " EndPoint=(" << Elt->EndX << "," << Elt->EndY;
                    OutDebug << ")  Control points (" << Elt->BzCtrl1X << "," << Elt->BzCtrl1Y << ") and (" << Elt->BzCtrl2X << "," <<  Elt->BzCtrl2Y << ")\n";
                    Elt = Elt->next;
                }
                OutDebug << std::flush;
            }
			if ( path->closed == 0 )
			{
                if ( debug_level > 0 )
                {
                    OutDebug << "Path is NOT closed, skip !\n";
                }
                //  Path is not closed, do NOT process
                continue;
			}
            //  If path has bezier curves, break them in small segments, if not copy data to polyline path
            SplitBezierCurves(path, flat_factor*flat_factor);
            if ( debug_level > 2 )
            {
                OutDebug << "After flatening bezier, Shape " << nShape << " Path " << nPath << " with " << path->nElts << " HasBezier:" << path->hasBezier << "\n";
                OutDebug << "  Start (" << path->PolyLine->StartX << "," << path->PolyLine->StartY << ") Bounds (" << path->bounds[0] << "," << path->bounds[1] << ") --> (" << path->bounds[2] << "," << path->bounds[3] << "\n";
                NSVGPathElt *Elt = path->PolyLine->PathElts;
                for (int i = 0; i < path->PolyLine->nElts; i++)
                {
                    OutDebug << "   Elt " << i << " type " << Elt->Type << " EndPoint=(" << Elt->EndX << "," << Elt->EndY;
                    OutDebug << ")  Control points (" << Elt->BzCtrl1X << "," << Elt->BzCtrl1Y << ") and (" << Elt->BzCtrl2X << "," <<  Elt->BzCtrl2Y << ")\n";
                    Elt = Elt->next;
                }
            }

            SimplifyPath(path->PolyLine, flat_factor*flat_factor);

            if ( debug_level > 2 )
            {
                OutDebug << "After Simplification, Shape " << nShape << " Path " << nPath << " with " << path->nElts << " HasBezier:" << path->hasBezier << "\n";
                OutDebug << "  Start (" << path->PolyLine->StartX << "," << path->PolyLine->StartY << ") Bounds (" << path->bounds[0] << "," << path->bounds[1] << ") --> (" << path->bounds[2] << "," << path->bounds[3] << "\n";
                NSVGPathElt *Elt = path->PolyLine->PathElts;
                for (int i = 0; i < path->PolyLine->nElts; i++)
                {
                    OutDebug << "   Elt " << i << " type " << Elt->Type << " EndPoint=(" << Elt->EndX << "," << Elt->EndY;
                    OutDebug << ")  Control points (" << Elt->BzCtrl1X << "," << Elt->BzCtrl1Y << ") and (" << Elt->BzCtrl2X << "," <<  Elt->BzCtrl2Y << ")\n";
                    Elt = Elt->next;
                }
            }
            //  Create polygon structure
            path->ClosePolygon = new Polygon();
            Path2Polygon(path->ClosePolygon, path->PolyLine);
            if ( debug_level > 1 )
            {
                OutDebug << "Adding polygon " << nPath << " new name " << path->id << " with " << path->ClosePolygon->nVertices << " vertices, area:" << path->ClosePolygon->area() << "mm²\n";
            }
            if ( path->ClosePolygon->isClockWise() )
            {
                if ( debug_level  > 1 )
                {
                    OutDebug << "Polygon is clockwise, make it counter clockwise\n";
                }
                path->ClosePolygon->Reverse();
            }
		}
		if ( debug_level > 0)
            OutDebug << std::flush;
	}
	if ( KeepNested )
	{
        nShape = 0;
        //  Second pass, check if a path is included in an other one
        if ( debug_level > 0 )
        {
            OutDebug << "Transform paths pass 2 : checking if paths are included in other ones\n";
        }
        for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
        {
            nPath = 0;
            NSVGpath *OldPath = NULL;
            NSVGpath *next_path = NULL;
            //  For each path of the shape, check if it is included in a larger o
            for (NSVGpath *path = shape->paths; path != NULL; path = next_path, nPath++)
            {
                next_path = path->next;     //  Memorize value, because will be modified if path is removed
                if ( RemoveIfIncluded(path, shape, OldPath, nShape, nPath) == 0 )
                {
                //  Not removed, new OldPath now. IF removed, DO NOT change OldPath, as the current entry has been moved
                    OldPath = path;
                }
            }
        }
	}
    if ( debug_level > 0 )
    {
        OutDebug << "End of Transform paths\n" << std::flush;
    }

}

//  If a path is included in an other one, remove this path for the primary list and move it to the child list of the large path
//  Do NOT take into account nested inclusion, there is a single level.

int SvgDoc::RemoveIfIncluded(NSVGpath *Inpath, NSVGshape *Inshape, NSVGpath *OldPath, int nShapeIn, int nPathIn)
{
    int nShape = 0;
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        int nPath = 0;
        //  For each path of the shape, check if the path passed as parameter (Inpath) is included in this one
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next, nPath++)
		{
            if ( path->ClosePolygon != NULL && Inpath->ClosePolygon != NULL )
            {
                if ( nShapeIn == 0 && nPathIn == 1 && nShape== 0 && nPath == 2 )
                {
                    OutDebug << "Inpath->ClosePolygon->Poly_in_Poly(path->ClosePolygon) returns " << Inpath->ClosePolygon->Poly_in_Poly(path->ClosePolygon) << " \n";
                }
                if ( Inpath->ClosePolygon->Poly_in_Poly(path->ClosePolygon) )
                {
                    if ( debug_level > 0 )
                    {
                        OutDebug << "Shape_" << nShapeIn << "/Path_" << nPathIn << "(" << Inpath->id << ") is included in " << path->id << "\n";
                    }
                    //  First, remove from primary list
                    if ( OldPath )
                    {
                        OldPath->next = Inpath->next;       //  Not first in list, change next element of previous one
                    }
                    else
                    {
                        Inshape->paths = Inpath->next;      //  First in list for this shape, change head of list
                    }
                    // Now move from primary list to child list, insert head of list
                    Inpath->next = path->child;
                    path->child = Inpath;
                    //  If Inpath has children, move them to child list
                    NSVGpath *pCur = Inpath->child;
                    while (  pCur != NULL )
                    {
                        NSVGpath *pNext = pCur->next;
                        pCur->next = path->child;
                        path->child = pCur;
                        pCur = pNext;
                    }
                    if ( debug_level > 1 )
                    {
                        OutDebug << path->id << " has " << NbChildren(path) << " children\n";
                    }
                    return 1;
                }
            }
            else if (  path->ClosePolygon == NULL )
            {
                OutDebug << " Shape " << nShapeIn << " Shape->ClosePolygon is NULL\n";
            }
            else
            {
                OutDebug << " InShape->ClosePolygon is NULL\n";
            }
        }
	}
	return 0;
}

void SvgDoc::EnlargePaths(double Diff)
{
    int nShape, nPath;

    if ( debug_level > 0 )
    {
        OutDebug << "---------------- Entering EnlargePaths, elapsed time: " <<  1000.0*(clock() - StartClock)/CLOCKS_PER_SEC << "ms\n";
    }
    nShape = 0;
    if ( SvgData == NULL ) return;
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        nPath = 0;
        //  Transform each path to a polygon
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next, nPath++)
		{
            if ( path->ClosePolygon == NULL )
            {
                OutDebug << path->id << " has no closed polygon, skipping large polygon generation\n";
                continue;
            }
            if ( debug_level > 1 )
            {
                OutDebug << path->id << " before elarging polygon with " << path->ClosePolygon->nVertices << " vertices and area " << path->ClosePolygon->area() << "mm²\n" << std::flush;
            }
            path->LargePolygon = path->ClosePolygon->enlarge(Diff);
            if ( debug_level > 1 )
            {
                OutDebug << path->id << " before simplify large polygon with " << path->LargePolygon->nVertices << " vertices and area " << path->LargePolygon->area() << "mm²\n" << std::flush;
            }
            path->LargePolygon->Simplify(Diff*Diff/5);
            if ( debug_level > 1 )
            {
                OutDebug << path->id << " after simplify large polygon with " << path->LargePolygon->nVertices << " vertices and area " << path->LargePolygon->area() << "mm²\n" << std::flush;
            }
		}
	}
    if ( debug_level > 0 )
    {
        OutDebug << "---------------- Exit EnlargePaths, elapsed time: " <<  1000.0*(clock() - StartClock)/CLOCKS_PER_SEC << "ms\n\n\n" << std::flush;
    }
}

//  For each path with a large polygon, modify this large polygon in order to keep max edge length lower tahn max_length
//  This will add vertices to polygons, so it will add some positions which could lead to better opimization.
//  But beware, more vertices also mean more execution time !

void SvgDoc::BreakLongerEdges(double max_length, double Diff)
{
    int nShape, nPath;

    if ( debug_level > 0 )
    {
        OutDebug << "---------------- Entering BreakLongerEdges, elapsed time: " <<  1000.0*(clock() - StartClock)/CLOCKS_PER_SEC << "ms\n";
    }
    nShape = 0;
    if ( SvgData == NULL ) return;
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        nPath = 0;
        //  Transform each large polygon to add vertices when edges are too long
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next, nPath++)
		{
            if ( path->LargePolygon == NULL )
            {
                OutDebug << " Shape " << shape->id << " Index " << nShape << " Path " << nPath << " has no large polygon\n";
                continue;
            }
            path->LargePolygon->BreakLongerEdges(max_length, Diff);
            if ( debug_level > 1 )
            {
                OutDebug << " Shape " << shape->id << " Index " << nShape << " Path " << nPath << " Large polygon with " << path->LargePolygon->nVertices << " vertices and area " << path->LargePolygon->area() << "mm²\n";
            }
		}
	}
    if ( debug_level > 0 )
    {
        OutDebug << "---------------- Exit BreakLongerEdges, elapsed time: " <<  1000.0*(clock() - StartClock)/CLOCKS_PER_SEC << "ms\n\n\n" << std::flush;
    }
}

//  Read SVG file and create Image object

SvgDoc::SvgDoc(string &FileName)
{

    StartClock = clock();
    Name = FileName;
	// Load SVG
	SvgData = nsvgParseFromFile(FileName.c_str(), "mm", 25.4);
	if (SvgData )
	{
        SheetSizeX = round(SvgData->width * 100) / 100.0;
        SheetSizeY = round(SvgData->height * 100) / 100.0;
	}
    nbTranslation = 0;
    nbRotation = 0;
    nbCheckAngles = 0;
    nbPointInPoly = 0;
    nbIntersectPoly = 0;
    nbPlacementImpossible = 0;
    CacheMiss = 0;
    CacheHit_OK = 0;
    CacheHit_KO = 0;

}

//  Return the number of children for a given path

int SvgDoc::NbChildren(NSVGpath *path)
{
int nb = 0;

    for ( NSVGpath *child = path->child; child != NULL; child = child->next)
    {
        nb++;
    }
    return nb;
}

SvgDoc::~SvgDoc()
{
    if ( OutDebug.is_open() )
        OutDebug.close();
}

void SvgDoc::WritePath(ostream &Out, NSVGpath *path, NSVGshape *shape)
{
    Out << "    <path" << endl;
    //  Line style from SVG file
    Out << std::hex;
    Out << "       style=\"fill:none;stroke:#" << setw(6) << setfill('0') << std::hex << (shape->stroke.color&0xFFFFFF);
    Out << std::dec;
    Out << ";stroke-width:" << shape->strokeWidth << "\"" << endl;
    // Then path itself
    Out << "       d=\"M " << path->StartX << "," << path->StartY;
    NSVGPathElt *Elt = path->PathElts;
    for (int i = 0; i < path->nElts; i++)
    {
        if ( Elt->Type == NSVG_PATH_BEZIER_C )
            Out << " C " << Elt->BzCtrl1X << "," << Elt->BzCtrl1Y << "," << Elt->BzCtrl2X << "," << Elt->BzCtrl2Y << "," << Elt->EndX << "," << Elt->EndY;
        else
            Out << " L " << Elt->EndX << "," << Elt->EndY;
        Elt = Elt->next;
    }
    if ( path->closed == 1 )
    {
        Out << " z";
    }
    Out << "\"" << endl;
    Out << "       id=\"" << path->id << "\"" << endl;
    Out << "       inkscape:connector-curvature=\"0\" />" << endl;
}

void SvgDoc::WritePlacedPath(ostream &Out, NSVGpath *ref_path, NSVGpath *out_path, NSVGshape *shape, int hasgroup)
{
    if (hasgroup)
        Out << "   ";
    Out << "    <path" << endl;
    //  Line style from SVG file
    Out << std::hex;
    if (hasgroup)
        Out << "   ";
    Out << "       style=\"fill:none;stroke:#" << setw(6) << setfill('0') << std::hex << (shape->stroke.color&0xFFFFFF);
    Out << std::dec;
    Out << ";stroke-width:" << shape->strokeWidth << "\"" << endl;
    // Get transformation values (rotation and translation) from reference path
    double angle = ref_path->PlacedPolygon->getRotation();
    Point PtTrans = ref_path->PlacedPolygon->getTranslation();
    Point Centroid = ref_path->PlacedPolygon->getCentroid() - PtTrans;      //  Should get centroid of unplaced polygon
    //  Then write path
    Point Start = Point(out_path->StartX, out_path->StartY);
    Start.Transform(angle, Centroid, PtTrans);
    if ( hasgroup )
        Out << "   ";
    Out << "       d=\"M " << Start.x << "," << Start.y;
    NSVGPathElt *Elt = out_path->PathElts;
    for (int i = 0; i < out_path->nElts; i++)
    {
        if ( Elt->Type == NSVG_PATH_BEZIER_C )
        {
            Point EltEnd = Point(Elt->EndX, Elt->EndY);
            EltEnd.Transform(angle, Centroid, PtTrans);
            Point EltBzCtrl1 = Point(Elt->BzCtrl1X, Elt->BzCtrl1Y);
            EltBzCtrl1.Transform(angle, Centroid, PtTrans);
            Point EltBzCtrl2 = Point(Elt->BzCtrl2X, Elt->BzCtrl2Y);
            EltBzCtrl2.Transform(angle, Centroid, PtTrans);
            Out << " C " << EltBzCtrl1.x << "," << EltBzCtrl1.y << "," << EltBzCtrl2.x << "," << EltBzCtrl2.y << "," << EltEnd.x << "," << EltEnd.y;
        }
        else
        {
            Point EltEnd = Point(Elt->EndX, Elt->EndY);
            EltEnd.Transform(angle, Centroid, PtTrans);
            Out << " L " << EltEnd.x << "," << EltEnd.y;
        }
        Elt = Elt->next;
    }
    if ( out_path->closed == 1 )
    {
        Out << " z";
    }
    Out << "\"" << endl;
    if ( hasgroup )
        Out << "   ";
    Out << "       id=\"Placed_" << out_path->id << "\"" << endl;
    if ( hasgroup )
        Out << "   ";
    Out << "       inkscape:connector-curvature=\"0\" />" << endl;
}
void SvgDoc::WriteOrginalLayer(ostream& Out)
{
    //  Then the data, begin with original layer
    Out << "  <g\n";
    Out << "     inkscape:label=\"Orginal_Layer\"" << endl;
    Out << "     inkscape:groupmode=\"layer\"" << endl;
    Out << "     id=\"Orginal_Layer1\">" << endl;

    int nShape = 0;
    Out << std::fixed;
    Out << std::setprecision(6);
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            WritePath(Out, path, shape);       //  Write path itself.
            for ( NSVGpath *path_child = path->child; path_child != NULL; path_child = path_child->next)
            {
                WritePath(Out, path_child, shape);
            }

        }
    }
    Out << "  </g>" << endl;
}

//  Write the polygon layer to output stream
void SvgDoc::WritePolygonLayer(ostream& Out)
{

    //  Then the data, begin with original layer
    Out << "  <g\n";
    Out << "     inkscape:label=\"Polygon_Layer\"" << endl;
    Out << "     inkscape:groupmode=\"layer\"" << endl;
    Out << "     id=\"Polygon_Layer1\">" << endl;

    int nShape = 0;
    Out << std::fixed;
    Out << std::setprecision(6);
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            Polygon *poly_path = path->ClosePolygon;
            if ( poly_path == NULL ) continue;
            Out << "    <path" << endl;
            //  Line style : cyan
            Out << "       style=\"fill:none;stroke:#00FFFF;stroke-width:" << shape->strokeWidth << "\"" << endl;
            // Then path itself
            Out << "       d=\"M " << poly_path->GetVertex(0)->x << "," << poly_path->GetVertex(0)->y;
			for (int i = 1; i < poly_path->nVertices; i++)
			{
                Out << " L " << poly_path->GetVertex(i)->x << "," << poly_path->GetVertex(i)->y;
			}
			if ( path->closed == 1 )
			{
                Out << " z";
			}
			Out << "\"" << endl;
            Out << "       id=\"poly_" << path->id << "\"" << endl;
            Out << "       inkscape:connector-curvature=\"0\" />" << endl;
        }
    }
    Out << "  </g>" << endl;
}

//  Write large polygon layer to output stream

void SvgDoc::WriteLargePolygonLayer(ostream& Out)
{
    //  Then the data, begin with original layer
    Out << "  <g\n";
    Out << "     inkscape:label=\"Large_Polygon_Layer\"" << endl;
    Out << "     inkscape:groupmode=\"layer\"" << endl;
    Out << "     id=\"LargePolygon_Layer1\">" << endl;

    int nShape = 0;
    Out << std::fixed;
    Out << std::setprecision(6);
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            Polygon *poly_path = path->LargePolygon;
            if ( poly_path == NULL ) continue;
            Out << "    <path" << endl;
            //  Line style : Blue
            Out << "       style=\"fill:none;stroke:#0000FF;stroke-width:" << shape->strokeWidth << "\"" << endl;
            // Then path itself
            Out << "       d=\"M " << poly_path->GetVertex(0)->x << "," << poly_path->GetVertex(0)->y;
			for (int i = 1; i < poly_path->nVertices; i++)
			{
                Out << " L " << poly_path->GetVertex(i)->x << "," << poly_path->GetVertex(i)->y;
			}
			if ( path->closed == 1 )
			{
                Out << " z";
			}
			Out << "\"" << endl;
            Out << "       id=\"Largepoly_" << shape->id << "_" << iPath << "\"" << endl;
            Out << "       inkscape:connector-curvature=\"0\" />" << endl;
        }
    }
    Out << "  </g>" << endl;
}

//  Write Hull to output file

void SvgDoc::WriteHullLargePolygonLayer(ostream& Out)
{
    //  Then the data, begin with original layer
    Out << "  <g\n";
    Out << "     inkscape:label=\"Hull_Placed_Layer\"" << endl;
    Out << "     inkscape:groupmode=\"layer\"" << endl;
    Out << "     id=\"Hull_Placed_Layer1\">" << endl;

    int nShape = 0;
    Out << std::fixed;
    Out << std::setprecision(6);
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            Polygon *poly_path = path->LargeConvexHull;
            if ( poly_path == NULL ) continue;
            Out << "    <path" << endl;
            //  Line style : Green
            Out << "       style=\"fill:none;stroke:#00FF00;stroke-width:" << shape->strokeWidth << "\"" << endl;
            // Then path itself
            Out << "       d=\"M " << poly_path->GetVertex(0)->x << "," << poly_path->GetVertex(0)->y;
			for (int i = 1; i < poly_path->nVertices; i++)
			{
                Out << " L " << poly_path->GetVertex(i)->x << "," << poly_path->GetVertex(i)->y;
			}
			if ( path->closed == 1 )
			{
                Out << " z";
			}
			Out << "\"" << endl;
            Out << "       id=\"Hull_" << shape->id << "_" << iPath << "\"" << endl;
            Out << "       inkscape:connector-curvature=\"0\" />" << endl;
        }
    }
    Out << "  </g>" << endl;

}

//  Write Placed POLYGON layer to output file

void SvgDoc::WritePlacedPolygonLayer(ostream& Out)
{
    //  Then the data, begin with original layer
    Out << "  <g\n";
    Out << "     inkscape:label=\"Placed_Polygon_Layer\"" << endl;
    Out << "     inkscape:groupmode=\"layer\"" << endl;
    Out << "     id=\"Placed_Polygon_Layer1\">" << endl;

    int nShape = 0;
    Out << std::fixed;
    Out << std::setprecision(6);
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            Polygon *poly_path = path->PlacedPolygon;
            if ( poly_path == NULL ) continue;
            Out << "    <path" << endl;
            //  Line style : Magenta
            Out << "       style=\"fill:none;stroke:#FF00FF;stroke-width:" << shape->strokeWidth << "\"" << endl;
            // Then path itself
            Out << "       d=\"M " << poly_path->GetVertex(0)->x << "," << poly_path->GetVertex(0)->y;
			for (int i = 1; i < poly_path->nVertices; i++)
			{
                Out << " L " << poly_path->GetVertex(i)->x << "," << poly_path->GetVertex(i)->y;
			}
			if ( path->closed == 1 )
			{
                Out << " z";
			}
			Out << "\"" << endl;
            Out << "       id=\"PP_" << shape->id << "_" << iPath << "\"" << endl;
            Out << "       inkscape:connector-curvature=\"0\" />" << endl;
        }
    }
    Out << "  </g>" << endl;

}

//  Write the result to output stream
void SvgDoc::WritePlacedLayer(ostream& Out)
{

    Out << "  <g\n";
    Out << "     inkscape:label=\"Placed_Layer\"" << endl;
    Out << "     inkscape:groupmode=\"layer\"" << endl;
    Out << "     id=\"Placed_Layer1\">" << endl;


    int nShape = 0;
    Out << std::fixed;
    Out << std::setprecision(6);
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            int hasgroup = 0;
            if ( path->PlacedPolygon == NULL ) continue;        //  Not placed, skip this one
            if ( path->child != NULL )          //  There is at least one child, build a group
            {
                Out << "    <g\n";
                Out << "       id=\"group_" << shape->id << "\">\n";
                hasgroup = 1;
            }
            WritePlacedPath(Out, path, path, shape, hasgroup);
            for ( NSVGpath *path_child = path->child; path_child != NULL; path_child = path_child->next)
            {
                WritePlacedPath(Out, path, path_child, shape, 1);
            }
            if ( hasgroup )
            {
                Out << "    </g>\n";
            }
        }
    }
    Out << "  </g>" << endl;
}

//  Write SVG header

void SvgDoc::WriteHeader(ostream &Out)
{
    Out << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    Out << "<!-- Created with CutOptim (http://www.fablab-lannion.org/) -->\n";
    Out << "\n";
    Out << "<svg\n";
    Out << "   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n";
    Out << "   xmlns:cc=\"http://creativecommons.org/ns#\"\n";
    Out << "   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n";
    Out << "   xmlns:svg=\"http://www.w3.org/2000/svg\"\n";
    Out << "   xmlns=\"http://www.w3.org/2000/svg\"\n";
    Out << "   xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"\n";
    Out << "   xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"\n";
    Out << "   width=\"" << round(SheetSizeX) << "mm\"" << endl;
    Out << "   height=\"" << round(SheetSizeY) << "mm\"\n";
    Out << "      viewBox=\"0 0 " << round(SheetSizeX) << " " << round(SheetSizeY) << "\"" << endl;
    Out << "      version=\"1.1\"\n";
    Out << "      id=\"svg8\"\n";
    Out << "      sodipodi:docname=\"" << Name << "\">" << endl;
    Out << "     <defs\n";
    Out << "        id=\"defs2\" />\n";
    Out << "     <sodipodi:namedview\n";
    Out << "        id=\"base\"\n";
    Out << "        pagecolor=\"#ffffff\"\n";
    Out << "        bordercolor=\"#666666\"\n";
    Out << "        borderopacity=\"1.0\"\n";
    Out << "        inkscape:pageopacity=\"0.0\"\n";
    Out << "        inkscape:pageshadow=\"2\"\n";
    Out << "        inkscape:document-units=\"mm\"\n";
    Out << "        inkscape:current-layer=\"Placed_Layer\"\n";
    Out << "        showgrid=\"false\"\n";
    Out << "        inkscape:window-maximized=\"1\" />\n";
    Out << "     <metadata\n";
    Out << "        id=\"metadata5\">\n";
    Out << "       <rdf:RDF>\n";
    Out << "         <cc:Work\n";
    Out << "            rdf:about=\"\">\n";
    Out << "           <dc:format>image/svg+xml</dc:format>\n";
    Out << "           <dc:type\n";
    Out << "              rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n";
    Out << "           <dc:title></dc:title>\n";
    Out << "         </cc:Work>\n";
    Out << "       </rdf:RDF>\n";
    Out << "     </metadata>\n";
}

void SvgDoc::WriteFile(ostream &Out, int output_layers)
{
    WriteHeader(Out);
    if ( output_layers & 1 )        //  Original layer
        WriteOrginalLayer(Out);
    if ( output_layers & 2 )        //  Close Polygon layer
        WritePolygonLayer(Out);
    if ( output_layers & 4 )        //  Large Polygon layer
        WriteLargePolygonLayer(Out);
    if ( output_layers & 8 )        //  Hull Placed layer
        WriteHullLargePolygonLayer(Out);
    if ( output_layers & 16 )        //  Placed Polygon layer
        WritePlacedPolygonLayer(Out);
    WritePlacedLayer(Out);
    Out << "</svg>" << endl;
}

int SvgDoc::WriteDoc(string &FileName, int FlagFile, int Output_Layers)
{

ofstream Out;

    if ( SvgData == NULL ) return(0);
    if ( FlagFile )
    {
        Out.open(FileName.c_str());
        if ( !Out.is_open()  ) return(0);
        WriteFile(Out, Output_Layers);
        Out.close();
    }
    else
        WriteFile(cout, Output_Layers);
    return 1;
}

//  Comparison by area

bool compare_area (const NSVGpath *first, const NSVGpath *second)
{
    return(first->LargePolygon->area() > second->LargePolygon->area());
}



//  Build a single list of all paths, sorted with large area first
void SvgDoc::BuilSingleListPath()
{
    int nShape = 0;
	for (NSVGshape *shape = SvgData->shapes; shape != NULL; shape = shape->next, nShape++)
	{
        //  Then for each path in the shape
        int iPath = 0;
        for (NSVGpath *path = shape->paths; path != NULL; path = path->next, iPath++)
        {
            path->isFixed = 0;              //  Not yet fixed !
            path->CachePoly = NULL;         //  No cache associated yet
            if ( path->LargePolygon == NULL ) continue;     //  Do NOT Add to list if no large polygon
            listPath.push_back(path);
        }
    }
    if ( debug_level > 0 )
    {
        OutDebug << "--------  BuilSingleListPath : There are " << listPath.size() << " paths in the list" << endl;
    }
    listPath.sort(compare_area);
    if ( debug_level > 0 )
    {
        int idx = 0;
        for (std::list<NSVGpath *>::iterator it=listPath.begin(); it != listPath.end(); ++it)
        {
            OutDebug << "  Path " << idx << " id: " << (*it)->id << " area " <<  (*it)->LargePolygon->area() << "mm2 and " << NbChildren(*it) << " children\n";
            idx++;
        }
        OutDebug << "-----------   End of BuilSingleListPath  ---------------------------\n\n\n";
    }
}
//  Sort polygons by area, larger first

void SvgDoc::SortbyArea()
{
    listPath.sort(compare_area);
}

//  For each path (polygon) compute the convex hull of this polygon

void SvgDoc::ComputeConvexHulls()
{
    for (std::list<NSVGpath *>::iterator it=listPath.begin(); it != listPath.end(); ++it)
    {
        NSVGpath *path = *it;
        path->LargeConvexHull = path->LargePolygon->ConvexHull();
    }
}

bool CheckOKAngles(double alpha, double beta, double a)
{
    //  First case, alpha is negative
    if ( alpha < 0 )
    {
        if ( beta < alpha )
        {
            //  Beta is also < 0
            //  In this case, OK if a >= alpha or a <= beta
            return (a >= alpha - angle_precision || a <= beta + angle_precision);
        }
        //  For all other cases, OK if a is betwen alpha and beta
        return (a >= alpha - angle_precision && a <= beta + angle_precision);
    }
    else
    {
    //  Alpha is positive.
        if ( beta >= alpha )
        {
            //  In this case beta is also positive, OK between alpha and beta
            return ( a >= alpha - angle_precision && a <= beta + angle_precision);
        }
        //  All other cases, NOK is between beta and alpha, so OK if >= alpha or <= beta
        return ( a >= alpha - angle_precision || a <= beta + angle_precision );
    }
}
//  Check if it is possible to place vertex idx_m of polygon mP at vertex idx_f of polygon fP
//  The edges of the moving polygon should NOT be inside the polygon fP

static bool CheckAngles(Polygon *mP, int idx_m, Polygon *fP, int idx_f)
{
double sm1, sm2;        //  Slopes for moving (new) polygon at requested vertex (just before and after)
double sf1, sf2;        //  Slopes for fixed polygon at requested vertex (just before and after)

    if ( idx_m == 0)
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sm1 = mP->getSlope(mP->nVertices-2) + M_PI;
        if ( sm1 > M_PI ) sm1 -= 2*M_PI;
        sm2 = mP->getSlope(0);
    }
    else
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sm1 = mP->getSlope(idx_m-1) + M_PI;
        if ( sm1 > M_PI ) sm1 -= 2*M_PI;
        sm2 = mP->getSlope(idx_m);
    }
    if ( idx_f == 0)
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sf1 = fP->getSlope(fP->nVertices-2) + M_PI;
        if ( sf1 > M_PI ) sf1 -= 2*M_PI;
        sf2 = fP->getSlope(0);
    }
    else
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sf1 = fP->getSlope(idx_f-1) + M_PI;
        if ( sf1 > M_PI ) sf1 -= 2*M_PI;
        sf2 = fP->getSlope(idx_f);
    }
    //  check if sm1 is OK, given angles sf1 and sf2
    int res1 = CheckOKAngles(sf1, sf2, sm1);
    if ( res1 == 0 ) return false;
    //  Idem with sm2
    int res2 = CheckOKAngles(sf1, sf2, sm2);
    if (res2 == 0 ) return false;
    //  Also check sf1 upon sm1 and sm2
    int res3 = CheckOKAngles(sm1, sm2, sf1);
    if ( res3 == 0 ) return false;
    //  and finally sf2 upon sm1 and sm2
    int res4 = CheckOKAngles(sm1, sm2, sf2);
    if (res4 == 0 ) return false;
    //  All is OK, return true
    return true;
}

//  Compute the angle index iRot when vertex iVertexFloat of float polygon FloatPoly if placed on vertex iVertexFixed of fixed polygon FixedPoly
//  Check if the other edge angle is OK with this rotation.
//  If position is OK, return 1, if not return 0
//  The parameter angle is upated by this function.

bool ComputeFreeRotAngle(Polygon *FloatPoly, int iVertexFloat, Polygon *FixedPoly, int iVertexFixed, int iRot, double &angle)
{
double sf1, sf2;
double sm1, sm2;

    //  First compute angles of edges ogf Fixed polygon
    if ( iVertexFixed == 0)
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sf1 = FixedPoly->getSlope(FixedPoly->nVertices-2) + M_PI;
        if ( sf1 > M_PI ) sf1 -= 2*M_PI;
        sf2 = FixedPoly->getSlope(0);
    }
    else
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sf1 = FixedPoly->getSlope(iVertexFixed-1) + M_PI;
        if ( sf1 > M_PI ) sf1 -= 2*M_PI;
        sf2 = FixedPoly->getSlope(iVertexFixed);
    }
    //  Then angles for the floating polygon
    if ( iVertexFloat == 0)
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sm1 = FloatPoly->getSlope(FloatPoly->nVertices-2) + M_PI;
        if ( sm1 > M_PI ) sm1 -= 2*M_PI;
        sm2 = FloatPoly->getSlope(0);
    }
    else
    {
        //  Arrival segment, angle is pi + slope, and if angle is greater than pi, substract 2*pi to keep value between -pi and +pi
        sm1 = FloatPoly->getSlope(iVertexFloat-1) + M_PI;
        if ( sm1 > M_PI ) sm1 -= 2*M_PI;
        sm2 = FloatPoly->getSlope(iVertexFloat);
    }
    switch ( iRot )
    {
    case 0:
    //  if iRot is 0, align sm1 on sf1
        angle = sf1 - sm1;
        break;
    case 1:
    //  if iRot is 1, align sm2 on sf1
        angle = sf1 - sm2;
        break;
    case 2:
    //  if iRot is 2, align sm1 on sf2
        angle = sf2 - sm1;
        break;
    case 3:
    //  if iRot is 3, align sm2 on sf2
        angle = sf2 - sm2;
        break;
    default:
        printf("Invalid choice !!\n");
        FloatPoly = 0;
        FloatPoly->nVertices = 0;           //  Segment violation
    }
    sm1 += angle;
    sm2 += angle;
    //  check if sm1 is OK, given angles sm1, smé, sf1 and sf2
    int res1 = CheckOKAngles(sf1, sf2, sm1);
    if ( res1 == 0 ) return false;
    //  Idem with sm2
    int res2 = CheckOKAngles(sf1, sf2, sm2);
    if (res2 == 0 ) return false;
    //  Also check sf1 upon sm1 and sm2
    int res3 = CheckOKAngles(sm1, sm2, sf1);
    if ( res3 == 0 ) return false;
    //  and finally sf2 upon sm1 and sm2
    int res4 = CheckOKAngles(sm1, sm2, sf2);
    if (res4 == 0 ) return false;
    //  All is OK, return true
    return true;
}

//  Optimize placement procedure when rotation is chosen among fixed angles
//  parameters
//  StepAngle : roration step in radian, when 0, not roration is allowed
//  Optimizing level : optimize placement of a group ogf optimizing level polygons. Beware, value of Optimizing level above 3 will lead to very slow processing
//  FirstPos : hint to lace first (largest polygon).
//  Flag_File : if true, output in a file, so information output on consoleis allowed.
//
//  Return : cost (area of convex hull) of all placed polygons.

double SvgDoc::Optimize(double StepAngle, int OptimizingLevel, int FirstPos, int Flag_file)
{
std::list<NSVGpath *> FloatingPaths = listPath;          //  All path are floating
std::list<NSVGpath *> FixedPath;
std::list<NSVGpath *> CurFloating;
Polygon *FixedPolyHull;
double BestCost = 1e100;
clock_t StartOptimize;
clock_t StartLevel;
clock_t LastLevel;
int NotPlaced = 0;
int num_rot = int((2*M_PI-0.001) / StepAngle) + 1;

    LastLevel = StartOptimize = clock();
    if ( debug_level > 0 )
    {
        OutDebug << "Time: " << fixed <<  ((double) 1000.0*(StartOptimize - StartClock))/CLOCKS_PER_SEC << "ms Enter Optimize with Optimizing level set to " << OptimizingLevel << endl;
    }
    //  First remove largest polygon from floating list
    NSVGpath *First = FloatingPaths.front();
    FloatingPaths.pop_front();
    //  Place first polygon, specific algorithm
    if ( debug_level > 0 )
    {
        OutDebug << "-----   Placing first polygon (" << First->id << ") -----------\n";
    }
    First->PlacedPolygon = PlaceFirst(First->LargePolygon, FirstPos, StepAngle);
    //  Check placement OK ! Should stay in the sheet.
    if ( First->PlacedPolygon == NULL )
    {
        cerr << "Unable to place first polygon, too large for sheet\n";
        if ( debug_level > 0 )
        {
            OutDebug << "Optimize : Unable to place first polygon, too large for sheet\n";
        }
        if ( Flag_file )
        {
            cout << "Unable to complete optimization, Unable to place first polygon, too large for sheet\n";
        }
        return 0;
    }
    //  Compute slopes of all edges of first polygon
    //  Should be done after polygon is placed because rotation change angles
    First->PlacedPolygon->ComputeAngles();
    First->PlacedPolygon->CalcSegments();
    First->isFixed = true;                                  //  OK Now !
    double FixedArea = First->PlacedPolygon->area();

    FixedPolyHull = First->LargeConvexHull = First->PlacedPolygon->ConvexHull();
    FixedPath.push_back(First);
    int idxPoly = 1;
    while ( !FloatingPaths.empty() )
    {
        StartLevel = clock();
        if ( debug_level > 0 )
        {
            OutDebug << "--------------------   Placing Polygon " << idxPoly << ": " << FloatingPaths.front()->id << "  ------------------------------\n";
            OutDebug << "Elapsed time " << ((double)1000.0*(StartLevel - StartOptimize))/CLOCKS_PER_SEC;
            OutDebug << "ms, for last polygon " << ((double)1000.0*(StartLevel - LastLevel))/CLOCKS_PER_SEC << "ms\n";
        }
        LastLevel = StartLevel;
        //  Prepare a list with at most OptimizingLevel elements which is a copy of the OptimizingLevel first elements
        int index = 0;
        for (std::list<NSVGpath *>::iterator it=FloatingPaths.begin(); it != FloatingPaths.end() && index < OptimizingLevel; ++it, index++)
        {
            CurFloating.push_back(*it);         //  Put OptimizingLevel elements in working list
        }
        //  Limit optimizing level when number of ploygon is low
        if ( OptimizingLevel > index ) OptimizingLevel = index ;

        NSVGpath *NewFixed = OptimizeLevel(FixedPolyHull, FixedPath, CurFloating, num_rot, StepAngle, FixedArea, BestCost, 1e100, idxPoly, idxPoly-1, OptimizingLevel);
        if ( NewFixed != NULL )                 //  If NULL, couldn't be placed, continue with next one.
        {
            NewFixed->isFixed = true;
            FixedPath.push_back(NewFixed);    //    Put it in fixed list
            //  Could (and should) free cache memory for this one
            delete NewFixed->CachePoly;
            NewFixed->CachePoly = NULL;
            FixedArea += NewFixed->LargePolygon->area();
            FixedPolyHull = NewFixed->LargeConvexHull;
        }
        else
        {
            NotPlaced++;                    //  COunt not placed elements
            cerr << "Unable to place polygon " << idxPoly << "\n";
            if ( debug_level > 0)
            {
                OutDebug << "######   Unable to place polygon " << idxPoly << ", aborting\n";
            }
        }

        FloatingPaths.pop_front();          //  Remove first in list
        CurFloating.clear();                //  Erase list before creating new one
        idxPoly++;
        if ( Flag_file )
        {
            StartLevel = clock();
            cout << "Polygon " << idxPoly << "  Not placed " << NotPlaced << "  Elapsed time " << ((double)1000.0*(StartLevel - LastLevel))/CLOCKS_PER_SEC << "\n";
#if CHECK_MEM > 0
            cout << " mem Poly_empty=" << nbEmptyPoly << "/" << nbCreatedEmptyPoly <<  " mem Poly_Vector=" << nbVectorPoly << "/" << nbCreatedVectorPoly;
            cout <<  " Rot_Poly=" << nbRotatedPoly << "/" << nbCreatedRotatedPoly;
            cout <<  " Trans_Poly=" << nbTranslatedPoly << "/" << nbCreatedTranslatedPoly;
//            cout <<  " Transformed_Poly=" << nbTransformedPoly << "/" << nbCreatedTransformedPoly;
            cout <<  " Cached_Poly=" << nbCachedPoly << "/" << nbCreatedCachedPoly;
            cout << "\n";
#endif
        }
    }
    if ( debug_level > 0 )
    {
        OutDebug << "Optimization finished, elapsed time " << ((double) (clock() - StartOptimize)) / CLOCKS_PER_SEC << "s\n";
        OutDebug << "nbRotation: " << nbRotation << "\n";
        OutDebug << "nbTranslation: " << nbTranslation << "\n";
        OutDebug << "nbPointInPoly: " << nbPointInPoly << "\n";
        OutDebug << "nbIntersectPoly: " << nbIntersectPoly << "\n";
        OutDebug << "nbPlacementImpossible: " << nbPlacementImpossible << "\n";
        OutDebug << "nbCheckAngles=" << nbCheckAngles << "\n";
        OutDebug << "CacheMiss=" << CacheMiss << "\n";
        OutDebug << "CacheHit OK=" << CacheHit_OK<< "\n";
        OutDebug << "CacheHit KO=" << CacheHit_KO<< "\n";
    }
    if ( Flag_file )
    {
        cout << "Optimization finished, cost " << BestCost << " elapsed time " << ((double) (clock() - StartOptimize)) / CLOCKS_PER_SEC << "\n";
        cout << "Performance stats :\n";
        cout << "nbTranslation=" << nbTranslation << "\n";
        cout << "nbRotation=" << nbRotation << "\n";
        cout << "nbCheckAngles=" << nbCheckAngles << "\n";
        cout << "nbPointInPoly=" << nbPointInPoly << "\n";
        cout << "nbIntersectPoly=" << nbIntersectPoly << "\n";
        cout << "nbPlacementImpossible=" << nbPlacementImpossible << "\n";
        cout << "CacheMiss=" << CacheMiss << "\n";
        cout << "CacheHit OK=" << CacheHit_OK<< "\n";
        cout << "CacheHit KO=" << CacheHit_KO<< "\n";
        if ( NotPlaced == 0 )
        {
            cout << "All elements have been successfully placed \n";
        }
        else
        {
            cout << NotPlaced << " elements couldn't have been placed \n";
        }
    }
    return BestCost;
}

//  Optimize the placement of the elements in list CurFloating.
//  The Fixed paths list and Hull are also passed as parameters
//  Return the best placement (least cost) for the FIRST element in CurFloating, or NULL if it could not be placed.
//  For each element in CurFloating, for all possible angles, try to place a floating polygon on a vertex of a fixed one (or semi fixed indeed)
//  This is a brute force algorithm !!
//  Try to keep running time reasonnable by trying only possible placements.

NSVGpath *SvgDoc::OptimizeLevel(Polygon *FixedPolyHull, std::list<NSVGpath *> FixedList, std::list<NSVGpath *> CurFloating, int num_rot, double StepAngle, double FixedArea,
                                double &Cost, double InBestCost, int FloatPolygon, int nbFixedPoly, int LevelOptimize)
{
Polygon *rot_p = NULL;
Polygon *trans_rot_p = NULL;
vector<Point> FixedHull;
vector<Point> tmpFixedHull;
double BestCost = 1e100;
int Best_i = -1, Best_j = -1;
int Best_Poly = -1;
double BestAngle;
Point BestTranslation;
Point OverAllCentroid;
Rectangle OverAll;
int Reason;
int FixedPolygon = 0;
Point RefPoint;

    if ( debug_level > 1 )
    {
        OutDebug << "Entering OptimizeLevel(" << LevelOptimize <<  "): FixedHullArea " <<  FixedPolyHull->area() << " Fixed area " << FixedArea;
        OutDebug << " InBestCost=" << InBestCost << "\n";
        OutDebug << "Elements in list:" << CurFloating.size() << " first is " << CurFloating.front()->id << "\n";
    }
    FixedHull = FixedPolyHull->GetVertices();
    OverAllCentroid = FixedPolyHull->getCentroid();
    OverAll = FixedPolyHull->GetBoundingBox();

    NSVGpath *WorkingFloatingEntry = CurFloating.front();
    //  And remove element from temporary floating list
    CurFloating.pop_front();
    BestCost = 1e100;          //  Very large, each try should be better
    BestAngle = -1;            //   Impossible value

    if ( debug_level > 1 )
    {
        OutDebug << "Polygon to be placed : " << WorkingFloatingEntry->id << "  has " << WorkingFloatingEntry->LargePolygon->nVertices << " vertices.\n";
    }
    double tmpPolygonArea = FixedArea + WorkingFloatingEntry->LargePolygon->area();
    int nPreviousVertices = getFixedNbVertex(FixedList);            //  Number of vertices of all polygons before this one
    //  Build a cache object which will speed up computation if not yet present
    if ( WorkingFloatingEntry->CachePoly == NULL )
    {
        WorkingFloatingEntry->CachePoly = new CachePosition(WorkingFloatingEntry->LargePolygon->nVertices, num_rot, nPreviousVertices);
    }
#ifdef UNDEF
    if ( FloatPolygon == 2)
    {
        OutDebug << "Cache object for Polygon 2=" << WorkingFloatingEntry->CachePoly << " cache value for 0/3/2 = " << WorkingFloatingEntry->CachePoly->isOKPlaced(3, 0, 2, &trans_rot_p)<< "\n";
    }
#endif
    //  Try each possible angle
    int iRot = 0;
    for ( double angle = 0; angle < 2*M_PI - 0.001; angle += StepAngle, iRot++)
    {
        if ( debug_level > 1 )
            OutDebug << "L" << LevelOptimize << ": Rotation set to " << round(angle*180/M_PI) << "\n";
        if ( rot_p ) delete rot_p;
        rot_p = WorkingFloatingEntry->LargePolygon->Rotate(angle);
        nbRotation++;
        rot_p->ComputeAngles();     //  Compute angles to be more efficient later
        rot_p->CalcSegments();

        //  Then try all vertices of all placed polygons
        FixedPolygon = 0;
        int iFixedVertex = 0;
        for (std::list<NSVGpath *>::iterator it=FixedList.begin(); it != FixedList.end(); ++it, FixedPolygon++)
        {
            //  Get next placed polygon
            NSVGpath *path = *it;
            Polygon  *curFixed = path->PlacedPolygon;
            // For each vertex of this placed polygon
            for ( int i = 0; i < curFixed->nVertices-1; iFixedVertex++, i++)
            {
                //  Try all vertices of the floating polygon, and select the best
                for (int j = 0; j < rot_p->nVertices-1; j++)
                {
#ifdef UNDEF
                    cout << "Passage=" << DebugPassage << " iRot=" << iRot << " Fixed=" << FixedPolygon << " i=" << i << " j=" << j << " rot_p->nVertices=" << rot_p->nVertices << "\n";
                    if ( iRot == 5 && i == 22 && j == 3 )
                    {
                        DebugPassage++;
                        cout << " Test rot_p (" << rot_p << ") rot_p->nVertices=" << rot_p->nVertices << "\n";
                    }
                    if ( FloatPolygon == 2 && j == 3 && i == 2 && FixedPolygon == 0 )
                    {
                        if ( debug_level > 2 )
                        {
                            OutDebug << "Test Check Angles, align vertex " << j << " on vertex " << i << " : " << *curFixed->GetVertex(i) << "\n";
                        }
                    }
#endif
                    //  Use cache if available
                    trans_rot_p = NULL;
                    int Cached = WorkingFloatingEntry->CachePoly->isOKPlaced(j, iRot, iFixedVertex, &trans_rot_p);
#ifdef UNDEF
                    if (FloatPolygon == 2 && iRot == 0 && iFixedVertex == 2 && j == 3 )
                    {
                        OutDebug << "Vertex 3 of polygon 2 on Vertex 2 of Polygon 0, Cached =" << Cached << "  trans_rot_p = " << trans_rot_p << "\n";
                    }
#endif
                    if ( Cached < 0 )       //  Cache says impossible !
                    {
                        CacheHit_KO++;
                        if ( debug_level > 2 )
                            OutDebug << "Cache says placement polygon " << FloatPolygon << " vertex " << j << " on vertex " << i << " of Polygon " << FixedPolygon << "(" << iFixedVertex << ")" << "impossible " << *curFixed->GetVertex(i) << "\n";
                        continue;           //  Goto next vertex...
                    }
                    else if ( Cached == 0 )
                    {
                        CacheMiss++;        // Not yet in cache, check angles first
                        nbCheckAngles++;
                        if ( rot_p->nVertices == 0 )
                        {
                            cout << " rot_p (" << rot_p << ") invalid with 0 vertices, level =" << LevelOptimize << " Float poly = " << FloatPolygon << "\n";
                        }
                        if (CheckAngles(rot_p, j, curFixed, i) == 0)
                        {
                            //  Impossible, invalid cache entry
                            WorkingFloatingEntry->CachePoly->addOKPlaced(j, iRot, iFixedVertex, NULL);
                            if ( debug_level > 2 )
                                OutDebug << "Placement polygon " << FloatPolygon << " vertex " << j << " on vertex " << i << " of Polygon " << FixedPolygon << " impossible " << *curFixed->GetVertex(i) << ", reason angle mismatch\n";
                            continue;    //  Not possible goto next without further processing
                        }
                        //  Compute polygon translated and rotated...
                        trans_rot_p = rot_p->Translate(curFixed->GetVertex(i), j);
                        nbTranslation++;
                        trans_rot_p->CalcSegments();
                    }
                    else
                    {
                        CacheHit_OK++;
                    }

                    RefPoint = Point(WorkingFloatingEntry->StartX, WorkingFloatingEntry->StartY);
                    RefPoint.Transform(angle, WorkingFloatingEntry->LargePolygon->getCentroid(), trans_rot_p->getTranslation());
                    if ( (Reason = PlacementNotPossible(trans_rot_p, FixedList, &OverAll, &RefPoint, Cached)) )
                    {
                        if ( debug_level > 2 )
                            OutDebug << "Placement polygon " << FloatPolygon << " vertex " << j << " on vertex " << i << " of Polygon " << FixedPolygon << " impossible " << *curFixed->GetVertex(i) << ", reason " << Reason << "\n";
                        //  If impossible due to conflict with fixed polygon, invalidate cache entry
                        if ( path->isFixed && (Reason & FIXED_BIT) != 0 )
                        {
                            WorkingFloatingEntry->CachePoly->addOKPlaced(j, iRot, iFixedVertex, NULL);
                            delete trans_rot_p;
                            trans_rot_p = NULL;
                        }
                        //  Free memory associated with polygon
                        if ( Cached == 0 && trans_rot_p != NULL )
                        {
                            delete trans_rot_p;
                            trans_rot_p = NULL;
                        }

                        nbPlacementImpossible++;
                        if ( rot_p->nVertices == 0 )
                        {
                            cout << " rot_p (" << rot_p << ") invalid with 0 vertices, level =" << LevelOptimize << " Float poly = " << FloatPolygon << "\n";
                        }
                        continue;
                    }
                    if ( path->isFixed && Cached == 0)        //  If current vertex belongs to a fixed path, cache value
                    {
                        WorkingFloatingEntry->CachePoly->addOKPlaced(j, iRot, iFixedVertex, trans_rot_p);
                        Cached = 1;
#ifdef UNDEF
                        if ( LevelOptimize == 1 && FloatPolygon == 2 && iRot == 0 && iFixedVertex == 2 && j == 3 )
                        {
                            OutDebug << "iRot=0, Vertex 3 of Polygon 2 placed on vertex 2 of polygon 0, add cache entry "  << trans_rot_p << "\n";
                            Polygon *TestPoly = NULL;
                            int TestCache = WorkingFloatingEntry->CachePoly->isOKPlaced(j, iRot, iFixedVertex, &TestPoly);
                            OutDebug << "Check Cached=" << TestCache << " Value Cached=" << TestPoly << "\n";
                        }
#endif
                    }
                    //  Compute cost function
                    //  We will use only the Hull area.
                    //  As the Hull will only grow with subsequent calls, if a temp value is greater than the best one, no need to process further polygons
                    vector<Point> tmpHull = CombineHull(FixedHull, trans_rot_p->GetVertices());
                    double tmpHullArea = HullArea(tmpHull);
                    if ( debug_level > 1)
                        OutDebug << "  Vertex " << j << " on " << i << " of polygon " << FixedPolygon << " Hull area " << tmpHullArea << "--> DiffArea=" << tmpHullArea - tmpPolygonArea << "\n";
                    //  In this case, call again this function but with a list without the first element
                    //  First copy the list, the first element has already been removed
                    std::list<NSVGpath *> tmpFloatList = CurFloating;
                    //  Then add this element, as placed to the fixed list.
                    std::list<NSVGpath *> tmpFixedList = FixedList;
                    tmpFixedList.push_back(WorkingFloatingEntry);
                    //  Move Polygon at its place, but before delete existing polygon to save memory
                    if ( WorkingFloatingEntry->PlacedPolygon ) delete WorkingFloatingEntry->PlacedPolygon;
                    if ( WorkingFloatingEntry->LargeConvexHull ) delete WorkingFloatingEntry->LargeConvexHull;
                    Point temp_pT = trans_rot_p->getTranslation();
                    WorkingFloatingEntry->PlacedPolygon = WorkingFloatingEntry->LargePolygon->Transform(angle, &temp_pT);
                    WorkingFloatingEntry->PlacedPolygon->ComputeAngles();
                    WorkingFloatingEntry->PlacedPolygon->CalcSegments();
                    //  Update variables for Fixed polygons lists
                    double tmpFixedArea = FixedArea + WorkingFloatingEntry->LargePolygon->area();
                    tmpFixedHull = CombineHull(FixedHull, WorkingFloatingEntry->PlacedPolygon->GetVertices());
                    WorkingFloatingEntry->LargeConvexHull = new Polygon(tmpFixedHull);
                    WorkingFloatingEntry->LargeConvexHull->ReCalcArea();
                    if ( debug_level > 1 && !tmpFloatList.empty() && tmpHullArea >= InBestCost )
                    {
                        OutDebug << " TmpHullarea (" << tmpHullArea << ") larger than input best values (" << InBestCost << ") no need to process other polygons in list\n";
                    }
                    if ( debug_level > 1 && !tmpFloatList.empty() && tmpHullArea >= BestCost )
                    {
                        OutDebug << " TmpHullarea (" << tmpHullArea << ") larger than current best values (" << BestCost << ") no need to process other polygons in list\n";
                    }
                    if ( !tmpFloatList.empty() && (tmpHullArea <= InBestCost && tmpHullArea <= BestCost) ) //  No need to go deeper, if nothing in list or cost already too high
                    {
                        if ( debug_level > 1 )
                        {
                            OutDebug << "  Applying current translation : " << temp_pT << " to polygon " << FloatPolygon << "\n";
                            OutDebug << "  new hull with " << FixedHull.size() << " vertices, area " << HullArea(FixedHull) << " Fixed polygon area " << tmpFixedArea <<"\n";
                            OutDebug << "  LargeConvexHull area :" << WorkingFloatingEntry->LargeConvexHull->area() << "\n";
                            OutDebug << "  Then call Optimize level again...\n";
                        }
                        OverAll = Hull2BoundingBox(FixedHull);
                        //  Then call Optimize level if the floating list is non empty
                        NSVGpath *Res = OptimizeLevel(WorkingFloatingEntry->LargeConvexHull, tmpFixedList, tmpFloatList, num_rot, StepAngle, tmpFixedArea, tmpHullArea, BestCost, FloatPolygon+1, nbFixedPoly, LevelOptimize - 1);
                        if ( Res == NULL )
                        {
                            //  Placement impossible, make sure that tmpHullArea will be greater than BestCost
                            tmpHullArea = BestCost + 10;
                        }
                        //  After return, the tmpHullArea is updated
                        if ( debug_level > 1 )
                        {
                            OutDebug << "-----   Returning at level " << LevelOptimize << ",  Cost is " << tmpHullArea << "\n";
                        }
                    }

                    if ( tmpHullArea < BestCost  )
                    {   //  Yes record new best value and corresponding position
                        BestCost = tmpHullArea;
                        BestAngle = angle;
                        BestTranslation = trans_rot_p->getTranslation();
                        Best_i = i;
                        Best_j = j;
                        Best_Poly = FixedPolygon;
                        if ( debug_level > 1 )
                        {
                            OutDebug << "New best found placing polygon " << FloatPolygon << "\n";
                            OutDebug << "  Vertex " << j << " on vertex " << i << " of polygon " <<  FixedPolygon << "\n";
                            OutDebug << " Rotation=" << round(angle*180/M_PI) << " Translation " << BestTranslation << " from " << *(rot_p->GetVertex(j)) << " to " <<  *(curFixed->GetVertex(i)) << "\n";
                            OutDebug << "  Hull area:" << tmpHullArea << " Polygon area:" << tmpPolygonArea << "\n";

                        }
                    }
                    if ( Cached == 0 )
                    {
                        delete trans_rot_p;     //  No need to keep, free memory
                        if ( rot_p->nVertices == 0 )
                        {
                            cout << " rot_p (" << rot_p << ") invalid with 0 vertices, level =" << LevelOptimize << " Float poly = " << FloatPolygon << "\n";
                        }
                    }
                }
             }
        }
    }
    //  Do some cleaning...
    if ( rot_p->nVertices == 0 )
    {
        cout << " rot_p (" << rot_p << ") invalid with 0 vertices, level =" << LevelOptimize << " Float poly = " << FloatPolygon << "\n";
    }
    if ( rot_p != NULL ) delete rot_p;
    rot_p = NULL;
    //  Place polygon at BestAngle and Best Translation
    if ( BestAngle < 0 )
    {
        Cost = 1.0e50;        //      Return a large cost, all possible solutions will be better than this one
        return NULL;
    }
    if ( BestCost > InBestCost )    //  Not better at this level
    {
        Cost = BestCost;
        if ( debug_level > 1 )
        {
            OutDebug << "Level " << LevelOptimize << " No better move than " <<  InBestCost << " found, return at upper level\n";
        }
        return WorkingFloatingEntry;
    }
    if ( WorkingFloatingEntry->PlacedPolygon ) delete WorkingFloatingEntry->PlacedPolygon;
    if ( WorkingFloatingEntry->LargeConvexHull ) delete WorkingFloatingEntry->LargeConvexHull;
    WorkingFloatingEntry->PlacedPolygon = WorkingFloatingEntry->LargePolygon->Transform(BestAngle, &BestTranslation);
    WorkingFloatingEntry->PlacedPolygon->ComputeAngles();
    WorkingFloatingEntry->PlacedPolygon->CalcSegments();
    FixedArea += WorkingFloatingEntry->PlacedPolygon->area();
    FixedHull = CombineHull(FixedHull, WorkingFloatingEntry->PlacedPolygon->GetVertices());
    WorkingFloatingEntry->LargeConvexHull = new Polygon(FixedHull);
    WorkingFloatingEntry->LargeConvexHull->ReCalcArea();
    WorkingFloatingEntry->LargeConvexHull->getCentroid();
    Cost = BestCost;
    if ( debug_level > 0 )
    {
        OutDebug << "Level " << LevelOptimize << ": Applying best move : " << BestTranslation << "/" << round(BestAngle*180/M_PI) << "° to polygon " << FloatPolygon << " will place vertex " << Best_j << " on vertex " << Best_i << " of polygon " << Best_Poly << "\n";
        OutDebug << "  new hull with " << FixedHull.size() << " vertices, area " << HullArea(FixedHull) << " Fixed polygon area " << FixedArea <<"\n";
        OutDebug << "  LargeConvexHull area :" << WorkingFloatingEntry->LargeConvexHull->area() << " Cost= " << Cost << "\n";
        OutDebug << "  CurFloating.Size=" << CurFloating.size() << "\n";
        int idx_list = 1;
        for (std::list<NSVGpath *>::iterator it=CurFloating.begin(); it != CurFloating.end() ; ++it, idx_list++)
        {
            NSVGpath  *lPath = *it;
            OutDebug << " Element " << idx_list << ": Translation=" << lPath->PlacedPolygon->getTranslation() << " rot=" << round(lPath->PlacedPolygon->getRotation()*180/M_PI) ;
            OutDebug << " Hull area=" << lPath->LargeConvexHull->area() << "\n";
        }
    }
    return WorkingFloatingEntry;
}

Point SvgDoc::ComputeRefPoint(int FirstPos, int xSize, int ySize)
{
Point A;

    if ( FirstPos < CenterLeft)     //  Top
    {
        A.y = 0;
    }
    else if ( FirstPos < BottomLeft )
    {
        A.y = (SheetSizeY - ySize)/2;       //  Sheet center
    }
    else
    {
        A.y = SheetSizeY - ySize;           //  Bottom
    }
    if ( (FirstPos &3) == 0 )     //  Left
    {
        A.x = 0;
    }
    else if ( (FirstPos &3) == 1 )     //  Center ?
    {
        A.x = (SheetSizeX - xSize)/2;       //  Sheet center
    }
    else
    {
        A.x = SheetSizeX - xSize;           //  Right
    }
    return A;
}

Polygon *SvgDoc::PlaceFirst(Polygon *p, int FirstPos, double StepAngle)
{
Rectangle BBox;
Polygon *Rot = NULL;
Polygon *cur;

    if ( debug_level > 0 )
    {
        OutDebug << "Entering PlaceFirst...\n";
    }
    //  If possible place polygon at page center
    for (double angle = 0; angle < 2*M_PI - 0.001; angle += StepAngle)
    {
        if ( angle > 0 )
        {
            if (Rot != NULL) delete Rot;
            Rot = p->Rotate(angle);
            cur = Rot;
        }
        else
        {
            cur = p;
        }
        BBox = cur->GetBoundingBox();
        double xSize = BBox.B.x - BBox.A.x;
        double ySize = BBox.B.y - BBox.A.y;
        if ( xSize < SheetSizeX && ySize < SheetSizeY )
        {
            //  Possible, compute new position. new Bounding box point A will be in
            Point newA = ComputeRefPoint(FirstPos, xSize, ySize);
            Point Translation = newA - BBox.A;
            Translation.align_grid();
            Point newPt0 = *(cur->GetVertex(0)) + Translation;
            newPt0.align_grid();
            Polygon *Translated = cur->Translate(&newPt0, 0);        //  Move point 0 to this point
            if ( Rot != NULL ) delete Rot;
            if ( debug_level > 0 )
            {
                OutDebug << "First polygon placed with rotation " <<  round(angle*180/M_PI) << "° and translation " << newPt0 << "\n";
            }
            return Translated;
        }
    }
    return NULL;
}

//  Return 0 if placement of polygon CurPoly is feasible with the polygons which are in the list FixedPath
//  If not return the reason of failure
int SvgDoc::PlacementNotPossible(Polygon *CurPoly, std::list<NSVGpath *> FixedPath, Rectangle *OverAll, const Point *RefPoint, int Cached)
{
int PlacedPolygon = 0;
int FixedBit = 0;
int FixedVertex = 0;

    if ( Cached == 0 )      //  If cached, no need to recompute !
    {
        //  If bounding box is outside the sheet, return false
        Rectangle BBox = CurPoly->GetBoundingBox();
        if ( BBox.A.x < 0 ) return FAIL_OUTSIDE_SHEET | FIXED_BIT;
        if ( BBox.A.y < 0 ) return FAIL_OUTSIDE_SHEET | FIXED_BIT;
        if ( BBox.B.x > SheetSizeX ) return FAIL_OUTSIDE_SHEET | FIXED_BIT;
        if ( BBox.B.y > SheetSizeY ) return FAIL_OUTSIDE_SHEET | FIXED_BIT;
        //  If BBox outside OverAll return true
        if ( BBox.A.x >= OverAll->B.x ) return 0;
        if ( BBox.A.y >= OverAll->B.y ) return 0;
        if ( BBox.B.x <= OverAll->A.x ) return 0;
        if ( BBox.B.y <= OverAll->A.y ) return 0;
    }

    PlacedPolygon = 0;
    //  Now checks if some vertex of CurPoly lies into a polygon from the list
    FixedVertex = 0;
    for (std::list<NSVGpath *>::iterator it=FixedPath.begin(); it != FixedPath.end(); ++it, PlacedPolygon++)
    {
        NSVGpath *path = *it;
        if ( path->isFixed )
            FixedBit = FIXED_BIT;
        else
            FixedBit = 0;
        FixedVertex += path->LargePolygon->nVertices;
        if ( Cached != 0 && path->isFixed && Cached > FixedVertex) continue;       //  No need to reprocess fixed polygons already computed
        //  Check if start point of real path is inside polygon. This test is mandatory to avoid placing two polygons at the same place !
        if ( path->PlacedPolygon->isInPoly(RefPoint))
            return FAIL_VERTEX_IN_POLYGON | FixedBit;
        //  Then check for all vertices of polygon
        for ( int i = 0; i < CurPoly->nVertices; i++ )
        {
            nbPointInPoly++;
            if ( path->PlacedPolygon->isInPoly(CurPoly->GetVertex(i)))
                return FAIL_VERTEX_IN_POLYGON | FixedBit;
        }
    }

    PlacedPolygon = 0;
    FixedVertex = 0;
    //  Last, also checks if polygon overlap with segment crossing, without point inclusion.
    //  Check if one segment of CUrPoly intersect with one segment of one of the fixed polys
    for (std::list<NSVGpath *>::iterator it=FixedPath.begin(); it != FixedPath.end(); ++it, PlacedPolygon++)
    {
        NSVGpath *path = *it;
        if ( path->isFixed )
            FixedBit = FIXED_BIT;
        else
            FixedBit = 0;
        FixedVertex += path->LargePolygon->nVertices;
        if ( Cached != 0 && path->isFixed && Cached > FixedVertex) continue;       //  No need to reprocess fixed polygons already computed
        Polygon *fPoly = path->PlacedPolygon;
        nbIntersectPoly++;
        if ( CurPoly->Intersect(fPoly) )
            return FAIL_EDGE_INTERSECT | FixedBit;
    }
    return 0;
}

//  Return the number of vertex of all fixed polygons

int SvgDoc::getFixedNbVertex(std::list<NSVGpath *> FixedList)
{
int num = 0;

    //  For each element in list, add the number of vertices
    for (std::list<NSVGpath *>::iterator it=FixedList.begin(); it != FixedList.end(); ++it)
    {
        NSVGpath *path = *it;
        num += path->LargePolygon->nVertices;
    }
    return num;
}

//  Optimize placement procedure when rotation is chosen among fixed angles
//  parameters
//  Optimizing level : optimize placement of a group ogf optimizing level polygons. Beware, value of Optimizing level above 3 will lead to very slow processing
//  FirstPos : hint to lace first (largest polygon).
//  Flag_File : if true, output in a file, so information output on consoleis allowed.
//
//  Return : cost (area of convex hull) of all placed polygons.

double SvgDoc::OptimizeFreeRot(int OptimizingLevel, int FirstPos, int Flag_file)
{
std::list<NSVGpath *> FloatingPaths = listPath;          //  All path are floating
std::list<NSVGpath *> FixedPath;
std::list<NSVGpath *> CurFloating;
Polygon *FixedPolyHull;
double BestCost = 1e100;
clock_t StartOptimize;
clock_t StartLevel;
clock_t LastLevel;
int NotPlaced = 0;

    LastLevel = StartOptimize = clock();
    if ( debug_level > 0 )
    {
        OutDebug << "Time: " << fixed <<  ((double) 1000.0*(StartOptimize - StartClock))/CLOCKS_PER_SEC << "ms Enter Optimize with Optimizing level set to " << OptimizingLevel << endl;
    }
    //  First remove largest polygon from floating list
    NSVGpath *First = FloatingPaths.front();
    FloatingPaths.pop_front();
    //  Place first polygon, specific algorithm
    if ( debug_level > 0 )
    {
        OutDebug << "-----   Placing first polygon (" << First->id << ") -----------\n";
    }
    //  In this mode (free rotation), try first with rotation angles of 90°
    First->PlacedPolygon = PlaceFirst(First->LargePolygon, FirstPos, M_PI/2.0);
    if ( First->PlacedPolygon == NULL )     //  If not OK, try with 9°
        First->PlacedPolygon = PlaceFirst(First->LargePolygon, FirstPos, M_PI/20.0);
    //  Check placement OK ! Should stay in the sheet.
    if ( First->PlacedPolygon == NULL )
    {
        cerr << "Unable to place first polygon, too large for sheet\n";
        if ( debug_level > 0 )
        {
            OutDebug << "Optimize : Unable to place first polygon, too large for sheet\n";
        }
        if ( Flag_file )
        {
            cout << "Unable to complete optimization, Unable to place first polygon, too large for sheet\n";
        }
        return 0;
    }
    //  Compute slopes of all edges of first polygon
    //  Should be done after polygon is placed because rotation change angles
    First->PlacedPolygon->ComputeAngles();
    First->PlacedPolygon->CalcSegments();
    First->isFixed = true;                                  //  OK Now !
    double FixedArea = First->PlacedPolygon->area();

    FixedPolyHull = First->LargeConvexHull = First->PlacedPolygon->ConvexHull();
    FixedPath.push_back(First);
    int idxPoly = 1;
    while ( !FloatingPaths.empty() )
    {
        StartLevel = clock();
        if ( debug_level > 0 )
        {
            OutDebug << "--------------------   Placing Polygon (Free rotation)" << idxPoly << ": " << FloatingPaths.front()->id << "  ------------------------------\n";
            OutDebug << "Elapsed time " << ((double)1000.0*(StartLevel - StartOptimize))/CLOCKS_PER_SEC;
            OutDebug << "ms, for last polygon " << ((double)1000.0*(StartLevel - LastLevel))/CLOCKS_PER_SEC << "ms\n";
        }
        LastLevel = StartLevel;
        //  Prepare a list with at most OptimizingLevel elements which is a copy of the OptimizingLevel first elements
        int index = 0;
        for (std::list<NSVGpath *>::iterator it=FloatingPaths.begin(); it != FloatingPaths.end() && index < OptimizingLevel; ++it, index++)
        {
            CurFloating.push_back(*it);         //  Put OptimizingLevel elements in working list
        }
        //  Limit optimizing level when number of ploygon is low
        if ( OptimizingLevel > index ) OptimizingLevel = index ;

        NSVGpath *NewFixed = OptimizeLevelFreeRot(FixedPolyHull, FixedPath, CurFloating, FixedArea, BestCost, 1e100, idxPoly, idxPoly-1, OptimizingLevel);
        if ( NewFixed != NULL )                 //  If NULL, couldn't be placed, continue with next one.
        {
            NewFixed->isFixed = true;
            FixedPath.push_back(NewFixed);    //    Put it in fixed list
            //  Could (and should) free cache memory for this one
            delete NewFixed->CachePoly;
            NewFixed->CachePoly = NULL;
            FixedArea += NewFixed->LargePolygon->area();
            FixedPolyHull = NewFixed->LargeConvexHull;
        }
        else
        {
            NotPlaced++;                    //  COunt not placed elements
            cerr << "Unable to place polygon " << idxPoly << "\n";
            if ( debug_level > 0)
            {
                OutDebug << "######   Unable to place polygon " << idxPoly << ", aborting\n";
            }
        }

        if ( Flag_file )
        {
            StartLevel = clock();
            cout << "FreeRotOptim : Polygon " << idxPoly << "(" << FloatingPaths.front()->id << ")" << "  Not placed " << NotPlaced << "  Elapsed time " << ((double)1000.0*(StartLevel - LastLevel))/CLOCKS_PER_SEC << "\n";
#if CHECK_MEM > 0
            cout << " mem Poly_empty=" << nbEmptyPoly << "/" << nbCreatedEmptyPoly <<  " mem Poly_Vector=" << nbVectorPoly << "/" << nbCreatedVectorPoly;
            cout <<  " Rot_Poly=" << nbRotatedPoly << "/" << nbCreatedRotatedPoly;
            cout <<  " Trans_Poly=" << nbTranslatedPoly << "/" << nbCreatedTranslatedPoly;
//            cout <<  " Transformed_Poly=" << nbTransformedPoly << "/" << nbCreatedTransformedPoly;
            cout <<  " Cached_Poly=" << nbCachedPoly << "/" << nbCreatedCachedPoly;
            cout << "\n";
#endif
        }
        FloatingPaths.pop_front();          //  Remove first in list
        CurFloating.clear();                //  Erase list before creating new one
        idxPoly++;
    }
    if ( debug_level > 0 )
    {
        OutDebug << "Optimization finished, elapsed time " << ((double) (clock() - StartOptimize)) / CLOCKS_PER_SEC << "s\n";
        OutDebug << "nbRotation: " << nbRotation << "\n";
        OutDebug << "nbTranslation: " << nbTranslation << "\n";
        OutDebug << "nbPointInPoly: " << nbPointInPoly << "\n";
        OutDebug << "nbIntersectPoly: " << nbIntersectPoly << "\n";
        OutDebug << "nbPlacementImpossible: " << nbPlacementImpossible << "\n";
        OutDebug << "nbCheckAngles=" << nbCheckAngles << "\n";
        OutDebug << "CacheMiss=" << CacheMiss << "\n";
        OutDebug << "CacheHit OK=" << CacheHit_OK<< "\n";
        OutDebug << "CacheHit KO=" << CacheHit_KO<< "\n";
    }
    if ( Flag_file )
    {
        cout << "Optimization finished, cost " << BestCost << " elapsed time " << ((double) (clock() - StartOptimize)) / CLOCKS_PER_SEC << "\n";
        cout << "Performance stats :\n";
        cout << "nbTranslation=" << nbTranslation << "\n";
        cout << "nbRotation=" << nbRotation << "\n";
        cout << "nbCheckAngles=" << nbCheckAngles << "\n";
        cout << "nbPointInPoly=" << nbPointInPoly << "\n";
        cout << "nbIntersectPoly=" << nbIntersectPoly << "\n";
        cout << "nbPlacementImpossible=" << nbPlacementImpossible << "\n";
        cout << "CacheMiss=" << CacheMiss << "\n";
        cout << "CacheHit OK=" << CacheHit_OK<< "\n";
        cout << "CacheHit KO=" << CacheHit_KO<< "\n";
        if ( NotPlaced == 0 )
        {
            cout << "All elements have been successfully placed \n";
        }
        else
        {
            cout << NotPlaced << " elements couldn't have been placed \n";
        }
    }
    return BestCost;
}

//  Optimize the placement of the elements in list CurFloating.
//  The Fixed paths list and Hull are also passed as parameters
//  Return the best placement (least cost) for the FIRST element in CurFloating, or NULL if it could not be placed.
//  For each element in CurFloating, try to place a floating polygon on a vertex of a fixed one (or semi fixed indeed)
//  This is a brute force algorithm, with an heuristic, only the angles which will lead to colinear edges are used (at most 2 angles for each test)
//  Try to keep running time reasonnable by trying only possible placements.

NSVGpath *SvgDoc::OptimizeLevelFreeRot(Polygon *FixedPolyHull, std::list<NSVGpath *> FixedList, std::list<NSVGpath *> CurFloating, double FixedArea,
                                double &Cost, double InBestCost, int FloatPolygon, int nbFixedPoly, int LevelOptimize)
{
Polygon *rot_p = NULL;
Polygon *trans_rot_p = NULL;
vector<Point> FixedHull;
vector<Point> tmpFixedHull;
double BestCost = 1e100;
int Best_i = -1, Best_j = -1;
int Best_Poly = -1;
double BestAngle;
Point BestTranslation;
Point OverAllCentroid;
Rectangle OverAll;
int Reason;
int FixedPolygon = 0;
Point RefPoint;

    if ( debug_level > 1 )
    {
        OutDebug << "Entering OptimizeLevel(" << LevelOptimize <<  "): FixedHullArea " <<  FixedPolyHull->area() << " Fixed area " << FixedArea;
        OutDebug << " InBestCost=" << InBestCost << "\n";
        OutDebug << "Elements in list:" << CurFloating.size() << " first is " << CurFloating.front()->id << "\n";
    }
    FixedHull = FixedPolyHull->GetVertices();
    OverAllCentroid = FixedPolyHull->getCentroid();
    OverAll = FixedPolyHull->GetBoundingBox();

    NSVGpath *WorkingFloatingEntry = CurFloating.front();
    //  And remove element from temporary floating list
    CurFloating.pop_front();
    BestCost = 1e100;          //  Very large, each try should be better
    BestAngle = -500;            //   Impossible value
    if ( debug_level > 1 )
    {
        OutDebug << "Polygon to be placed : " << WorkingFloatingEntry->id << "  has " << WorkingFloatingEntry->LargePolygon->nVertices << " vertices.\n";
    }
    double tmpPolygonArea = FixedArea + WorkingFloatingEntry->LargePolygon->area();
    int nPreviousVertices = getFixedNbVertex(FixedList);            //  Number of vertices of all polygons before this one
    //  Build a cache object which will speed up computation if not yet present
    if ( WorkingFloatingEntry->CachePoly == NULL )
    {
        //  Use 4 for num_rot because 4 possibilities will be tested, even if 2 at most are possible.
        WorkingFloatingEntry->CachePoly = new CachePosition(WorkingFloatingEntry->LargePolygon->nVertices, 4, nPreviousVertices);
    }
    //  Precomputing of segments and angle to save time
    WorkingFloatingEntry->LargePolygon->ComputeAngles();     //  Compute angles to be more efficient later
    WorkingFloatingEntry->LargePolygon->CalcSegments();

    //  Try all vertices of all placed polygons
    FixedPolygon = 0;
    int iFixedVertex = 0;
    for (std::list<NSVGpath *>::iterator it=FixedList.begin(); it != FixedList.end(); ++it, FixedPolygon++)
    {
        //  Get next placed polygon
        NSVGpath *path = *it;
        Polygon  *curFixed = path->PlacedPolygon;
        if ( debug_level > 2 )
        {
            OutDebug << "\n\n-------    Trying matching polygon " << FloatPolygon << " with Polygon " << FixedPolygon << "\n";
        }
        // For each vertex of this placed polygon
        for ( int iVFixed = 0; iVFixed < curFixed->nVertices-1; iFixedVertex++, iVFixed++)
        {
            double angle = 0;
            if ( debug_level > 2 )
            {
                OutDebug << "\nTrying on vertex " << iVFixed << " of fixed polygon\n";
            }
            //  Try all vertices of the floating polygon, and select the best
            for (int iVFloat = 0; iVFloat < WorkingFloatingEntry->LargePolygon->nVertices-1; iVFloat++)
            {
                if ( debug_level > 3 )
                {
                    OutDebug << "\nTrying with vertex " << iVFloat << " on vertex " << iVFixed << "\n";
                }
                //  For each case, try 4 angles. Angles are computed such as edge of flaot polygon are aligned to edge of fixed polygon with matching vertex.
                //  As there is 2 edges for each vertex, there are 4 possible angles, even if at most 2 will be possible
                for ( int iRot = 0; iRot < 4; iRot++)
                {
                    //  Use cache if available

                    trans_rot_p = NULL;
                    int Cached = WorkingFloatingEntry->CachePoly->isOKPlacedFreeRot(iVFloat, iRot, iFixedVertex, &trans_rot_p);
//                    Cached = 0;
                    if ( Cached < 0 )       //  Cache says impossible !
                    {
                        CacheHit_KO++;
                        if ( debug_level > 2 )
                            OutDebug << "Cache says placement polygon " << FloatPolygon << " vertex " << iVFloat << " on vertex " << iVFixed << " of Polygon " << FixedPolygon << "(" << iFixedVertex << ")" << "impossible " << *curFixed->GetVertex(iVFixed) << "\n";
                        continue;           //  Goto next vertex...
                    }
                    else if ( Cached == 0 )
                    {
                        CacheMiss++;        // Not yet in cache, check angles first
                        nbCheckAngles++;
                        angle = 0;
                        //  Compute angle associated with floating polygon vertex j placed on vertex i of fixed polygon curFixed polygon (index i_rot)
                        //  If not possible, invalidate cache
                        int Res = ComputeFreeRotAngle(WorkingFloatingEntry->LargePolygon, iVFloat, curFixed, iVFixed, iRot, angle);
                        if ( debug_level > 3 )
                        {
                            OutDebug << "iRot =" << iRot << " angle=" << round(angle * 180 / M_PI);
                            OutDebug << "  ComputeFreeRotAngle: Res=" << Res << "\n";
                        }
                        if ( Res == 0)
                        {
                            //  Impossible, invalid cache entry
                            WorkingFloatingEntry->CachePoly->addOKPlacedFreeRot(iVFloat, iRot, iFixedVertex, NULL);
                            if ( debug_level > 2 )
                                OutDebug << "ComputeFreeRotAngle: Placement polygon " << FloatPolygon << " vertex " << iVFloat << " with rotation " << angle << "(" << iRot << ")" << " on vertex " << iVFixed << " of Polygon " << FixedPolygon << " impossible " << *curFixed->GetVertex(iVFixed) << ", reason angle mismatch\n";
                            continue;    //  Not possible goto next without further processing
                        }
                        //  Compute polygon translated and rotated...
                        if ( rot_p ) delete rot_p;
                        rot_p = WorkingFloatingEntry->LargePolygon->Rotate(angle);
                        nbRotation++;
                        rot_p->ComputeAngles();     //  Compute angles to be more efficient later
                        rot_p->CalcSegments();
                        trans_rot_p = rot_p->Translate(curFixed->GetVertex(iVFixed), iVFloat);
                        nbTranslation++;
                        trans_rot_p->CalcSegments();
                    }
                    else
                    {
                        CacheHit_OK++;
                        angle = trans_rot_p->getRotation();
                    }

                    RefPoint = Point(WorkingFloatingEntry->StartX, WorkingFloatingEntry->StartY);
                    RefPoint.Transform(angle, WorkingFloatingEntry->LargePolygon->getCentroid(), trans_rot_p->getTranslation());
                    if ( (Reason = PlacementNotPossible(trans_rot_p, FixedList, &OverAll, &RefPoint, Cached)) )
                    {
                        if ( debug_level > 2 )
                            OutDebug << "Placement polygon " << FloatPolygon << " vertex " << iVFloat << " on vertex " << iVFixed << " of Polygon " << FixedPolygon << " impossible " << *curFixed->GetVertex(iVFixed) << ", reason " << Reason << "\n";
                        //  If impossible due to conflict with fixed polygon, invalidate cache entry
                        if ( path->isFixed && (Reason & FIXED_BIT) != 0 )
                        {
                            WorkingFloatingEntry->CachePoly->addOKPlacedFreeRot(iVFloat, iRot, iFixedVertex, NULL);
                            delete trans_rot_p;
                            trans_rot_p = NULL;
                        }
                        //  Free memory associated with polygon
                        if ( Cached == 0 && trans_rot_p != NULL )
                        {
                            delete trans_rot_p;
                            trans_rot_p = NULL;
                        }

                        nbPlacementImpossible++;
                        continue;
                    }
                    if ( path->isFixed && Cached == 0)        //  If current vertex belongs to a fixed path, cache value
                    {
                        WorkingFloatingEntry->CachePoly->addOKPlacedFreeRot(iVFloat, iRot, iFixedVertex, trans_rot_p);
                        Cached = 1;
                    }
                    //  Compute cost function
                    //  We will use only the Hull area.
                    //  As the Hull will only grow with subsequent calls, if a temp value is greater than the best one, no need to process further polygons
                    vector<Point> tmpHull = CombineHull(FixedHull, trans_rot_p->GetVertices());
                    double tmpHullArea = HullArea(tmpHull);
                    if ( debug_level > 1)
                        OutDebug << "  Vertex " << iVFloat << " on " << iVFixed << " of polygon " << FixedPolygon << " Hull area " << tmpHullArea << "--> DiffArea=" << tmpHullArea - tmpPolygonArea << "\n";
                    //  In this case, call again this function but with a list without the first element
                    //  First copy the list, the first element has already been removed
                    std::list<NSVGpath *> tmpFloatList = CurFloating;
                    //  Then add this element, as placed to the fixed list.
                    std::list<NSVGpath *> tmpFixedList = FixedList;
                    tmpFixedList.push_back(WorkingFloatingEntry);
                    //  Move Polygon at its place, but before delete existing polygon to save memory
                    if ( WorkingFloatingEntry->PlacedPolygon ) delete WorkingFloatingEntry->PlacedPolygon;
                    if ( WorkingFloatingEntry->LargeConvexHull ) delete WorkingFloatingEntry->LargeConvexHull;
                    Point temp_pT = trans_rot_p->getTranslation();
                    WorkingFloatingEntry->PlacedPolygon = WorkingFloatingEntry->LargePolygon->Transform(angle, &temp_pT);
                    WorkingFloatingEntry->PlacedPolygon->ComputeAngles();
                    WorkingFloatingEntry->PlacedPolygon->CalcSegments();
                    //  Update variables for Fixed polygons lists
                    double tmpFixedArea = FixedArea + WorkingFloatingEntry->LargePolygon->area();
                    tmpFixedHull = CombineHull(FixedHull, WorkingFloatingEntry->PlacedPolygon->GetVertices());
                    WorkingFloatingEntry->LargeConvexHull = new Polygon(tmpFixedHull);
                    WorkingFloatingEntry->LargeConvexHull->ReCalcArea();
                    if ( debug_level > 1 && !tmpFloatList.empty() && tmpHullArea >= InBestCost )
                    {
                        OutDebug << " TmpHullarea (" << tmpHullArea << ") larger than input best values (" << InBestCost << ") no need to process other polygons in list\n";
                    }
                    if ( debug_level > 1 && !tmpFloatList.empty() && tmpHullArea >= BestCost )
                    {
                        OutDebug << " TmpHullarea (" << tmpHullArea << ") larger than current best values (" << BestCost << ") no need to process other polygons in list\n";
                    }
                    if ( !tmpFloatList.empty() && (tmpHullArea <= InBestCost && tmpHullArea <= BestCost) ) //  No need to go deeper, if nothing in list or cost already too high
                    {
                        if ( debug_level > 1 )
                        {
                            OutDebug << "  Applying current translation : " << temp_pT << " to polygon " << FloatPolygon << "\n";
                            OutDebug << "  new hull with " << FixedHull.size() << " vertices, area " << HullArea(FixedHull) << " Fixed polygon area " << tmpFixedArea <<"\n";
                            OutDebug << "  LargeConvexHull area :" << WorkingFloatingEntry->LargeConvexHull->area() << "\n";
                            OutDebug << "  Then call Optimize level again...\n";
                        }
                        OverAll = Hull2BoundingBox(FixedHull);
                        //  Then call Optimize level if the floating list is non empty
                        NSVGpath *Res = OptimizeLevelFreeRot(WorkingFloatingEntry->LargeConvexHull, tmpFixedList, tmpFloatList, tmpFixedArea, tmpHullArea, BestCost, FloatPolygon+1, nbFixedPoly, LevelOptimize - 1);
                        if ( Res == NULL )
                        {
                            //  Not possible, make sure that tmpHullArea will be greater than BestCost
                            tmpHullArea = BestCost + 10;
                        }
                        //  After return, the tmpHullArea is updated
                        if ( debug_level > 1 )
                        {
                            OutDebug << "-----   Returning at level " << LevelOptimize << ",  Cost is " << tmpHullArea << "\n";
                        }
                    }
                    if ( tmpHullArea < BestCost  )
                    {   //  Yes record new best value and corresponding position
                        BestCost = tmpHullArea;
                        BestAngle = angle;
                        BestTranslation = trans_rot_p->getTranslation();
                        Best_i = iVFixed;
                        Best_j = iVFloat;
                        Best_Poly = FixedPolygon;
                        if ( debug_level > 1 )
                        {
                            OutDebug << "New best found placing polygon " << FloatPolygon << "\n";
                            OutDebug << "  Vertex " << iVFloat << " on vertex " << iVFixed << " of polygon " <<  FixedPolygon << "\n";
                            OutDebug << " Rotation=" << round(angle*180/M_PI) << " Translation " << BestTranslation << " from " << *(rot_p->GetVertex(iVFloat)) << " to " <<  *(curFixed->GetVertex(iVFixed)) << "\n";
                            OutDebug << "  Hull area:" << tmpHullArea << " Polygon area:" << tmpPolygonArea << "\n";

                        }
                    }
                    if ( Cached == 0 )
                    {
                        delete trans_rot_p;     //  No need to keep, free memory
                        if ( rot_p->nVertices == 0 )
                        {
                            cout << " rot_p (" << rot_p << ") invalid with 0 vertices, level =" << LevelOptimize << " Float poly = " << FloatPolygon << "\n";
                        }
                    }
               }
            }
         }
    }
    //  Do some cleaning...
    if ( rot_p != NULL ) delete rot_p;
    rot_p = NULL;
    //  Place polygon at BestAngle and Best Translation
    if ( BestAngle < -100 )
    {
        Cost = 1.0e50;        //      Return a large cost, all possible solutions will be better than this one
        return NULL;
    }
    if ( BestCost > InBestCost )    //  Not better at this level
    {
        Cost = BestCost;
        if ( debug_level > 1 )
        {
            OutDebug << "Level " << LevelOptimize << " No better move than " <<  InBestCost << " found, return at upper level\n";
        }
        return WorkingFloatingEntry;
    }
    if ( WorkingFloatingEntry->PlacedPolygon ) delete WorkingFloatingEntry->PlacedPolygon;
    if ( WorkingFloatingEntry->LargeConvexHull ) delete WorkingFloatingEntry->LargeConvexHull;
    WorkingFloatingEntry->PlacedPolygon = WorkingFloatingEntry->LargePolygon->Transform(BestAngle, &BestTranslation);
    WorkingFloatingEntry->PlacedPolygon->ComputeAngles();
    WorkingFloatingEntry->PlacedPolygon->CalcSegments();
    FixedArea += WorkingFloatingEntry->PlacedPolygon->area();
    FixedHull = CombineHull(FixedHull, WorkingFloatingEntry->PlacedPolygon->GetVertices());
    WorkingFloatingEntry->LargeConvexHull = new Polygon(FixedHull);
    WorkingFloatingEntry->LargeConvexHull->ReCalcArea();
    WorkingFloatingEntry->LargeConvexHull->getCentroid();
    Cost = BestCost;
    if ( debug_level > 0 )
    {
        OutDebug << "Level " << LevelOptimize << ": Applying best move : " << BestTranslation << "/" << round(BestAngle*180/M_PI) << "° to polygon " << FloatPolygon << " will place vertex " << Best_j << " on vertex " << Best_i << " of polygon " << Best_Poly << "\n";
        OutDebug << "  new hull with " << FixedHull.size() << " vertices, area " << HullArea(FixedHull) << " Fixed polygon area " << FixedArea <<"\n";
        OutDebug << "  LargeConvexHull area :" << WorkingFloatingEntry->LargeConvexHull->area() << " Cost= " << Cost << "\n";
        OutDebug << "  CurFloating.Size=" << CurFloating.size() << "\n";
        int idx_list = 1;
        for (std::list<NSVGpath *>::iterator it=CurFloating.begin(); it != CurFloating.end() ; ++it, idx_list++)
        {
            NSVGpath  *lPath = *it;
            if ( lPath->PlacedPolygon == NULL )
            {
                OutDebug << "Element " << idx_list << "Not placed\n";
            }
            else
            {
                OutDebug << " Element " << idx_list << ": Translation=" << lPath->PlacedPolygon->getTranslation() << " rot=" << round(lPath->PlacedPolygon->getRotation()*180/M_PI) ;
                OutDebug << " Hull area=" << lPath->LargeConvexHull->area() << "\n";
            }
        }
    }
    return WorkingFloatingEntry;
}

