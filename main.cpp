#include <iostream>
#include "cxxopts.hpp"
#include "SvgDoc.h"
#include "Geometry.h"
#include <math.h>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;



cxxopts::ParseResult parse(int argc, char* argv[])
{
  try
  {
    cxxopts::Options options(argv[0], " - example command line options");
    options
      .positional_help("[optional args]")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f, file", "File", cxxopts::value<std::string>()->default_value("TestPoly1.svg"), "SVG Input File")
        ("o,output", "Output file", cxxopts::value<std::string>(), "SVG Output File")
        ("positional", "File to be processed", cxxopts::value<std::vector<std::string>>())
        ("h,help", "Print help")
        ("d,distance", "Min distance between paths to be cut", cxxopts::value<double>(), "1.0")
        ("q,rect_cost", "Overall Rectangle cost", cxxopts::value<double>(), "0.0")
        ("m,max_length", "Max length of one segment, break than longer", cxxopts::value<double>(), "1000.0")
        ("l,optimizing_level", "Optimizing level, process list_size elements together", cxxopts::value<int>(), "1")
        ("debug_level", "Level of debug info in specific debug file", cxxopts::value<int>(), "0")
        ("debug_file", "Generate debug info from inkscape", cxxopts::value<bool>()->default_value("true"))
        ("df_name", "Name of debug file", cxxopts::value<std::string>()->default_value("Debug_CutOptim.txt"))
        ("k,original", "Output Original layer", cxxopts::value<bool>()->default_value("false"))
        ("n,nested", "Keep nested path together", cxxopts::value<bool>()->default_value("true"))
        ("c,use_cache", "USe cache to speed up processing", cxxopts::value<bool>()->default_value("false"))
        ("y,layer_output", "Output internal layers : 1 Input layer, 2 Polygon, 4 Large polygon, 8 Hull layer, 16 Placed Polygon layer, OR these values to output multiple layers", cxxopts::value<int>(), "0")
        ("a,angle", "Rotation step", cxxopts::value<double>(), "90.0")
        ("r,free_rot", "allow free rotation", cxxopts::value<bool>()->default_value("true"))
        ("p,firstpos", "Position of largest object", cxxopts::value<std::string>(), "Position of largest object on the sheet")
    ;
    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help({""}) << std::endl;
      exit(0);
    }

    return result;

  } catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

static int ConvertPos(string Pos)
{
int rPos = 0;

    if ( Pos.length() < 2 ) return(TopLeft);
    if ( Pos[0] == 'T' || Pos[0] == 't')
    {
        rPos = 0;
    }
    else if ( Pos[0] == 'C' || Pos[0] == 'c' )
    {
        rPos = 4;
    }
    else if ( Pos[0] == 'B' || Pos[0] == 'b')
    {
        rPos = 8;
    }
    else
    {
        return TopLeft;         //  Top left if unknown
    }
    if ( Pos[1] == 'L' || Pos[1] == 'l')
    {
        rPos += 0;
    }
    else if ( Pos[1] == 'C' || Pos[1] == 'c' )
    {
        rPos += 1;
    }
    else if ( Pos[1] == 'R' || Pos[1] == 'r' )
    {
        rPos += 2;
    }
    else
    {
        return TopLeft;         //  Top left if unknown
    }
    return rPos;
}
#define TEST_ANGLES 0

extern bool CheckOKAngles(double sf1, double sf2, double sm1);

int main(int argc, char *argv[])
{
string InputFileName;
string OutputFileName;
double MinCutDistance = 1.0;
double StepAngle = 0;
double rect_cost = 0.0;
int Flag_file = 0;
int debug_level = 0;
int Output_Layer = 0;
int KeepNested = 1;
int OptimizingLevel = 1;
int FirstPos = CenterCenter;
int UseCache = 0;
bool Flag_free_rot = 1;
string debug_file_name = "Debug_CutOptim.txt";

double max_segment_length = 1000.0;
#ifdef UNDEF
        //  Used to check param passing from inkscape
        ofstream ParamDbg("/home/thierry/Programmes/CutOptim/DebgParam.txt");
        if ( ParamDbg.is_open() == false )
        {
            cerr << "unable to open /home/thierry/Programmes/CutOptim/DebgParam.txt\n";
            exit(1);
        }
         ParamDbg << "--- Params  ------\n";
        ParamDbg << " argc =" << argc << "\n";
        for ( int i = 0; i < argc; i++)
        {
            ParamDbg << "Argv[" << i << "]=" << argv[i] << "\n";
        }
        ParamDbg.close();
#endif
    auto result = parse(argc, argv);
    InputFileName = result["f"].as<std::string>();
    if (result.count("output"))
    {
        OutputFileName = result["output"].as<std::string>();
        Flag_file = 1;
    }
    if (result.count("distance"))
        MinCutDistance = result["distance"].as<double>();
    if (result.count("angle"))
    {
        StepAngle = result["angle"].as<double>();
        Flag_free_rot = 0;
    }
    if (result.count("max_length"))
        max_segment_length = result["max_length"].as<double>();
    if (result.count("rect_cost"))
        rect_cost = result["rect_cost"].as<double>();
    if (result.count("debug_level"))
        debug_level = result["debug_level"].as<int>();
    if (result.count("layer_output"))
        Output_Layer = result["layer_output"].as<int>();
    if (result.count("optimizing_level"))
        OptimizingLevel = result["optimizing_level"].as<int>();
    if (result.count("original"))
    {
        Output_Layer |= result["original"].as<bool>();
    }
    if ( result.count("firstpos"))
    {
        FirstPos = ConvertPos(result["firstpos"].as<string>() );
    }
    if ( result.count("df_name"))
    {
        debug_file_name = result["df_name"].as<string>();
    }
    if (result.count("nested"))
    {
        KeepNested = result["nested"].as<bool>();
    }

    if (result.count("use_cache"))
    {
        UseCache = result["use_cache"].as<bool>();
    }

    if (result.count("free_rot"))
    {
        Flag_free_rot = result["free_rot"].as<bool>();
    }
    if (result.count("debug_file"))
        debug_level += result["debug_file"].as<bool>();;
    if ( Flag_file )
    {
        cout << "Parse options complete" << endl;
        cout << "Input File = " << InputFileName << std::endl;
        cout << "Output File = " << OutputFileName << std::endl;
        cout << "angle = " << StepAngle << std::endl;
        cout << "MinCutDistance = " << MinCutDistance << std::endl;
        cout << "Debug Level = " << debug_level << std::endl;
        cout << "Output layers = " << Output_Layer << std::endl;
        cout << "Optimizing level = " << OptimizingLevel << std::endl;
        cout << "Max_segment_length = " << max_segment_length << std::endl;
        cout << "Use free rotation = " << Flag_free_rot << std::endl;
        cout << "Rectangle cost factor = " << rect_cost << std::endl;
        cout << "Keep nested paths = " << KeepNested << std::endl;
        cout << "Use Cache = " << UseCache << std::endl;
    }

      //    Read and process input file

    SvgDoc InputSVG(InputFileName);
    InputSVG.setDebugLevel(debug_level, debug_file_name);
    InputSVG.setUseCache(UseCache);
    InputSVG.setRectCost(rect_cost);
    if ( InputSVG.SvgData == NULL )
    {
    cerr << " No input to process aborting \n";
    return 1;
    }
    //    Change paths to polygones
    InputSVG.TransformPaths(MinCutDistance / 10.0, KeepNested);
    //    Make larger polygons
    InputSVG.EnlargePaths(MinCutDistance);
    //  Break longer edges if necessary
    InputSVG.BreakLongerEdges(max_segment_length, MinCutDistance / 10.0);
    //    Sort polygons by area
    InputSVG.BuilSingleListPath();
      //    Compute convex hulls of each path
      //  InputSVG.ComputeConvexHulls();
      //    Optimize
    if ( StepAngle == 0 ) StepAngle = 360;
    if ( Flag_free_rot )
        InputSVG.OptimizeFreeRot(OptimizingLevel, FirstPos, Flag_file);
    else
        InputSVG.Optimize(StepAngle * 2 * M_PI / 360, OptimizingLevel, FirstPos, Flag_file);
    InputSVG.WriteDoc(OutputFileName, Flag_file, Output_Layer);


#ifdef TEST_GEO
      //    Now test distance from point to Segment
      Point P00 = Point(0,0);
      Point P11 = Point(1,1);
      Segment Seg1 = Segment(P00, P11);
      Point P01 = Point(0, 1);
      Point P10 = Point(1, 0);
      printf("Distance from P10 to Seg1 =%.3f\n", sqrt(Seg1.sqrDistancePoint(P10)));
      printf("Distance from P01 to Seg1 =%.3f\n", sqrt(Seg1.sqrDistancePoint(P01)));
      Point P2 = Point(2,0);
      printf("Distance from P2 to Seg1 =%.3f\n", sqrt(Seg1.sqrDistancePoint(P2)));
      Point P3 = Point(-1,-1);
      printf("Distance from P3 to Seg1 =%.3f\n", sqrt(Seg1.sqrDistancePoint(P3)));
      Point P4 = Point(2,2);
      printf("Distance from P4 to Seg1 =%.3f\n", sqrt(Seg1.sqrDistancePoint(P4)));
      Point P5 = Point(-1,0);
      printf("Distance from P5 to Seg1 =%.3f\n", sqrt(Seg1.sqrDistancePoint(P5)));
      Segment Seg2 = Segment(P11, P00);
      printf("Distance from P1 to Seg2 =%.3f\n", sqrt(Seg2.sqrDistancePoint(P10)));
      printf("Distance from P2 to Seg2 =%.3f\n", sqrt(Seg2.sqrDistancePoint(P2)));
      printf("Distance from P3 to Seg2 =%.3f\n", sqrt(Seg2.sqrDistancePoint(P3)));
      printf("Distance from P4 to Seg2 =%.3f\n", sqrt(Seg2.sqrDistancePoint(P4)));
      printf("Distance from P5 to Seg2 =%.3f\n", sqrt(Seg2.sqrDistancePoint(P5)));

      //    Test Polygons
      Polygon Poly1 = Polygon();
      Polygon Poly2 = Polygon();
      //    First a rectangle
      Poly1.addVertice(P00);
      Poly1.addVertice(P01);
      Poly1.addVertice(P11);
      Poly1.addVertice(P10);
      Poly1.addVertice(P00);
      printf("Polygon 1 nVertices = %d clockwise =%d, area = %.3f\n", Poly1.nVertices, Poly1.isClockWise(), Poly1.area());

      Poly2.addVertice(P00);
      Poly2.addVertice(P10);
      Poly2.addVertice(P11);
      Poly2.addVertice(P01);
      Poly2.addVertice(P00);
      printf("Polygon 2 nVertices = %d clockwise =%d, area = %.3f\n", Poly2.nVertices, Poly2.isClockWise(), Poly2.area());
      printf("Distance P01 to Poly2 =%.3f\n", sqrt(Poly2.distance(P01)));
      printf("Distance P10 to Poly2 =%.3f\n", sqrt(Poly2.distance(P10)));
      printf("Distance P2 to Poly2 =%.3f\n", sqrt(Poly2.distance(P2)));
      printf("Distance P3 to Poly2 =%.3f\n", sqrt(Poly2.distance(P3)));
      printf("Distance P4 to Poly2 =%.3f\n", sqrt(Poly2.distance(P4)));
      printf("Distance P5 to Poly2 =%.3f\n", sqrt(Poly2.distance(P5)));

      Polygon Poly8 = Polygon();
      for (int i = 0; i <= 8; i++)
      {
        Point p = Point(cos(i*M_PI/4.0), sin(i*M_PI/4.0));
        Poly8.addVertice(p);
      }
      printf("Polygon 8 nVertices = %d clockwise =%d, area = %.3f\n", Poly8.nVertices, Poly8.isClockWise(), Poly8.area());

      printf("Distance P01 to Poly8 =%.3f\n", sqrt(Poly8.distance(P01)));
      printf("Distance P10 to Poly8 =%.3f\n", sqrt(Poly8.distance(P10)));
      printf("Distance P2 to Poly8 =%.3f\n", sqrt(Poly8.distance(P2)));
      printf("Distance P3 to Poly8 =%.3f\n", sqrt(Poly8.distance(P3)));
      printf("Distance P4 to Poly8 =%.3f\n", sqrt(Poly8.distance(P4)));
      printf("Distance P5 to Poly8 =%.3f\n", sqrt(Poly8.distance(P5)));
#endif
#if TEST_ANGLES > 0
    printf("Cas 1 : Départ en -1,0, arrivée en 0,0,  puis va en 0,1\n");
    double sf1 = atan2(-1, 0);
    double sf2 = atan2(-0.2, -1);

    if ( sf1 < 0 ) sf1 += 2*M_PI;
    if ( sf2 < 0 ) sf2 += 2*M_PI;
    double sm2 = 0;
    printf("sf1 = %.1f, sf2 = %.1f\n", sf1 * 180 / M_PI, sf2 *180.0 / M_PI);

    printf("Essai sm1 sm2\n");

    sm2 = 0;
    for ( double a = 1.0; a < 360; a+= 20.0)
    {
        sm2 = a/180.0 * M_PI;
        printf("Angle %.1f, SM1 %d  SM2 %d\n", a, CheckOKAngles(sf1, sf2, sm2), CheckOKAngles(sf2, sf1, sm2));
    }
#endif
      return(0);
}
