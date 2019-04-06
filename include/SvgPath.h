#ifndef SVGPATH_H
#define SVGPATH_H

#include "Geometry.h"

class PathElt
{
    public:
        int TypeElt;
        Point EndPoint;
        Point EltBezier1, EltBezier2;
        PathElt(int Type, Point End, Point Bezier1=Point(0,0), Point Bezier2=Point(0,0))
        {
            TypeElt = Type;
            EndPoint = End;
            EltBezier1 = Bezier1;
            EltBezier1 = Bezier2;
        }
};

class SvgPath
{
    public:
        SvgPath();
        virtual ~SvgPath();

        inline int isClosed() { return(ClosedPath); }



    protected:
        int ClosedPath;
        Point StartPath;
        Rectangle BoundingBox;
    private:
};

#endif // SVGPATH_H
