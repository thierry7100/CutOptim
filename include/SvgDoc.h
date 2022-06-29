#ifndef SVGDOC_H
#define SVGDOC_H

#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <ostream>
#include <iostream>
#include <fstream>
#include <ctime>
#include "nanosvg.h"

using namespace std;

class SvgDoc
{
    public:
        SvgDoc(string &FileName);
        virtual ~SvgDoc();

        NSVGimage *SvgData;
        double SheetSizeX, SheetSizeY;
        void TransformPaths(double flat_factor, bool KeepNested);
        void EnlargePaths(double Diff);
        void BreakLongerEdges(double max_length, double Diff);
        void BuilSingleListPath();
        void SortbyArea();
        void ComputeConvexHulls();
        double Optimize(double StepAngle, int OptimizingLevel, int FirstPos, int Flag_File);
        double OptimizeFreeRot(int OptimizingLevel, int FirstPos, int Flag_File);
        void setDebugLevel(int level, string DebugFileName) { debug_level = level; OutDebug.open(DebugFileName); }
        void setUseCache(int Val) { UseCache = Val;}
        void setRectCost(double Val) { RectCostFactor = Val;}

        int WriteDoc(string &FileName, int Flag_File, int Output_Layers);

        std::list<NSVGpath *> listPath;

    protected:
        NSVGpath *OptimizeLevel(Polygon *FixedPolyHull, std::list<NSVGpath *> FixedList, std::list<NSVGpath *> CurFloating, int num_rot, double StepAngle, double FixedArea, double &BestCost, double InBestCost, int FloatPolygon, int nbFixedPoly, int LevelOptimize);
        NSVGpath *OptimizeLevelFreeRot(Polygon *FixedPolyHull, std::list<NSVGpath *> FixedList, std::list<NSVGpath *> CurFloating, double FixedArea, double &BestCost, double InBestCost, int FloatPolygon, int nbFixedPoly, int LevelOptimize);
        void SimplifyPath(NSVGpath *path, double max_error);
        void WriteOrginalLayer(ostream& Out);
        void WritePolygonLayer(ostream& Out);
        void  WriteLargePolygonLayer(ostream& Out);
        void  WriteHullLargePolygonLayer(ostream& Out);
        void  WritePlacedPolygonLayer(ostream& Out);
        void  WritePlacedLayer(ostream& Out);
        void  WriteHeader(ostream& Out);
        void  WritePath(ostream &Out, NSVGpath *path, NSVGshape *shape);
        void  WritePlacedPath(ostream &Out, NSVGpath *ref_path, NSVGpath *out_path, NSVGshape *shape, int hasgroup);
        void  WriteFile(ostream &Out, int output_layers);
        Polygon *PlaceFirst(Polygon *p, int FirstPos, double StepAngle);
        int PlacementNotPossible(Polygon *CurPoly, std::list<NSVGpath *> FixedPath, Rectangle *OverAll, const Point *RefPoint, int Cached);
        int RemoveIfIncluded(NSVGpath *path, NSVGshape *shape, NSVGpath *OldPath, int nShape, int nPath );
        int NbChildren(NSVGpath *path);
        int getFixedNbVertex(std::list<NSVGpath *> FixedList);
        Point ComputeRefPoint(int FirstPos, int xSize, int ySize);

        double RectCostFactor;
        string Name;
        int debug_level;
        int UseCache;

        ofstream OutDebug;
        clock_t StartClock;

        long int nbTranslation;
        long int nbRotation;
        long int nbCheckAngles;
        long int nbPointInPoly;
        long int nbIntersectPoly;
        long int nbPlacementImpossible;

        long int CacheMiss;
        long int CacheHit_OK;
        long int CacheHit_KO;

    private:
};

#endif // SVGDOC_H
