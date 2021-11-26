## Presentation

All the users of laser cutter have probably been confronted with the problem: my drawing does not fit in the wood sheet at my disposal or even in the laser cutter! To try to solve this problem, I started writing an optimization program for placing objects on a sheet. 
For interested readers, this is a variant of a well known problem in operational research (bin packing problem). By doing a little search on the Web, there are also links to applications performing this task, but the free versions are often very limited, they are limited optimizing the cutting of rectangles. 
 
This program reads an input SVG file containing the objects to be placed and outputs a second SVG file containing the placed objects. 
It can be used as is (command line, no graphical interface!) Or as an extension inkscape which then provides the GUI. 

## Environment 

As mentioned above, this program can be used alone or as an inkscape extension. 
The processing can take a relatively long time, it is better to run this program on a modern processor, but if you are in no hurry ...
The memory consumption is reasonable, no need to rush to buy new RAMs!
The program has been tested both on Linux and Windows, so it should work on both systems.


## The Software

First of all as it is a program written in C ++ for a performance issue, so it must be compiled on your machine. 
I wrote this program under Linux / Ubuntu (compiled with gcc), but since there is no system dependency, it should work as is under any other version of Linux. For fans of Windows (there are so many!), I created a Visual Studio project that allows to compile on this platform. For Mac users, sorry I do not have access to a Mac, you will have to make the adaptation  yourself, but the C ++ used is really standard, it should work as soon as you have access to a compiler. For information, I did not change the code between Linux and Windows, that's saying! 

### Linux installation

The code is available here: https://github.com/thierry7100/CutOptim
For the beginner, you clone (or download) the directory, it comes in the form of a .zip archive, which must be extracted.
Then you open a terminal, go to the created directory and launch the commands:
- make release
- make install: this will copy the software into the directory  ~/.local / bin which is in the list of executable directories, which will allow you to use it directly (this may be specific to Ubuntu, to you to put the program elsewhere on another system.
- make install_inkscape: this will copy the program to the inkscape extension directory (~/.config/inkscape / extensions). If you want to make this extension available for all accounts on your machine, copy the file cutoptim.inx + the executable into /usr/share/inkscape / extensions (you must be root).

If you have opted for the inkscape extension, at the next start you will have a Fablab / Laser Cutting Optimizer extension


### Windows installation

The code is available here: https://github.com/thierry7100/CutOptim
For the uninitiated, you clone (or download) the directory, it comes in the form of a .zip archive, which must be extracted.
Then you launch Visual Studio, you can get a free version for a specific purpose, see https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=Community&rel=16
Then, once Visual Studio is started, you open the CutOptim project, then:
- You ask to generate the Release version of the project if it is not the one that appears in the menu bar.
- You choose the platform (x86 or x64) of your choice. By default the file is configured in x64, if you have a 32bit version of Windows, change to x86.
- Then choose Generate / Generate Solution, compilation starts and after a few seconds, your program is available.
- Then, under windows, it is better to use this program like extension inkscape, the command line hardly being used? ! To do this, copy the cutoptim.inx and CutOptm / x64 / Release / CutOptim.exe files to the inkscape extensions directory. This can be found via the Edit / Preferences / System command, but it is usually under C: \ USERS \\ Your user name\\ AppData \ Roaming / inkscape / extensions. Attention, to see this directory, you will have to validate the visualization of hidden files under the explorer of files, if it were not done.
As in Linux, the next startup of inkscape you will have a Fablab / Laser Cutting Optimizer

## Operation

The input and output format of the files is the SVG format, available on many software programs. If you use inkscape, it is the native format, the program has been tested in this context. It should also work with files from other drawing software that generates SVG. If you have a problem tell me (thierry@fablab-lannion.org) <br>

Description of the process:
- The inkscape document (at least its size) can generally be of importance, here it sets the size of the sheet used for cutting. You will need to set a document size compatible with your material (and the cutter, of course!).
- First, the program reads the inkscape document, it only considers paths or simple objects. The text strings not transformed into paths, the pictures ... are simply ignored. I therefore advise to turn everything into a path before launching the program. Inkscape Ctrl + A then Objects in Path (SHIFT + Ctrl + C). Unclosed paths are also ignored, the software is only able to process shapes with a closed outline.
- The paths can be placed anywhere, on the sheet or outside, it does not matter. To a certain extent, they can even be superimposed (see below).
- Then, from these paths, the program creates polygons approaching the paths (with an error of less than 0.1mm on average).
- Then the program "enlarges" these polygons to prevent paths from touching each other in the final result. The size of the enlargement is configurable.
- The program then takes these enlarged polygons and will try to place them in a way that is not optimal but good. Why not optimal because the problem is difficult (complete NP in mathematical terms) and requires a very long time even for simple configurations. The basic idea here is to start from the largest polygon, then to place the sorted polygons by decreasing size such that a vertex of the polygon to be placed is positioned on a vertex of an already placed polygon. This reduces the space of possibilities, even if it remains very large!
- The "best" configuration is obtained when the size of the convex hull is minimal. Another mathematical term! The convex hull is the smallest convex form containing all points of all plotted polygons. Intuitively, this maximizes the free space on the sheet, which is the desired result. Be careful, it is not necessarily the smallest rectangle, the convex hull is not usually a rectangle!
- To place the paths, the software is allowed to rotate the objects, unless you block this possibility. Depending on your needs (non-homogeneous material) you may have to limit rotations to 0 and 180 ° for example, or even to block any rotation (this will be the case for example with printed fabric).

### Program options as inkscape extension
![Options dialog box](CutOptim_optionsinkscape.png) 
The program has many options detailed below:
- Units: Always use mm, the program is not tested for other choices. Sorry for inches users, you will have to test this !
- Min distance between objects: This is the size at which the created polygons will be enlarged. This value must be greater than 0.8mm, the approximation by polygons is not perfect.
- Max length of single segment: As explained above the software will try to find a good configuration by positioning vertices on other vertices. It can be interesting in some cases to "add" vertices to have more possibilities. If an edge is longer than the specified size, it will be broken into multiple segments, with additional vertices. Do not abuse this option, too low a value will slow down the treatment tremendously. Do not go below 100mm in most cases, even if the value 0 is allowed to indicate that you do NOT want to use this possibility.
- Optimizing level: as indicated above the program places the polygons in order of decreasing size (we place the largest pebbles first ...). This sometimes leads to clearly suboptimal situations. By increasing this parameter, the software will optimize the placement of a group of N polygons. This gives better results, but be careful, it considerably increases the treatment time. Do not exceed 2 or 3, if the default value of 1 does not give good results. If you draw has 300 vertices already placed (rather low value actually) and you allow rotations in steps of 10 °, use N = 2 will multiply operations by 36 * 300! And for N = 3 by (36 * 300) ²!
- Keep original layer in output: If you do not trust the program (!), you can check this option, the original shapes will be kept as well as those placed, but placed in different layers. You will be able to check the work done.
- Select option for largest element placement: The first item can be placed where you want on the page. Usually at the top left, but the center also gives good results.
- Allow free rotation of paths, angle parameter not used : If this option is checked, the angle of rotation of each object will be chosen from 4 to make the edge coincide with one of the two of the vertex on which the object will be positioned. This option is economical in processing time (at most 4 tests) but can give less good results than the fixing of the angle of rotation. The results are worse when the segments are very short, if the input form is not a polygon for example.
- Try rotation by (0 no rotation allowed): This option is incompatible with the previous one, it is only valid if the previous option is NOT checked. In this case, the objects are positioned on discrete rotation steps. 0 means that rotations are prohibited, this is useful when the material is not homogeneous. For MDF, no restrictions, but for wood or even plywood, if you want to respect the direction of the wood, rotations are not advisable. Choose 180 ° in this case. Attention, low values ​​greatly increase the calculation time. With 10 °, there are 36 times more calculations than with 0 °! If the input shapes are rectangular, a value of 90 ° gives good results.
- Attach nested path to the bigger one :  In practice, we often deal with situations or related objects. For example, a plate with fixing holes. If this box is checked, the software checks if the path is included in another one, and if it is, it will not process it but link it to the larger path. Once it is placed, the same transformation (rotation / translation) will be applied to the "small" included object. Attention the software is not able to recover the space released in holes, you have to leave a little work anyway.
- Debug file generation: If this box is checked, a debug file (Debug_CutOptim.txt) is created in the inkscape extension directory. This can be used to understand what happened when it goes wrong.

### Comand line options

The program has the following options which are listed when typing CutOptim -h from the command line. For the explanations, please refer to the previous chapter, the command line options are the same as when used as an inkscape extension.

thierry@thierry-UX410UAR:~/Programmes/CutOptim$ CutOptim -h
 - example command line options
Usage:
  CutOptim [OPTION...] [optional args]

  -f, --file SVG Input File     File (default: TestPoly1.svg)
  -o, --output SVG Output File  Output file
      --positional arg          File to be processed
  -h, --help                    Print help
  -d, --distance 1.0            Min distance between paths to be cut
  -m, --max_length 1000.0       Max length of one segment, break than longer
  -l, --optimizing_level 1      Optimizing level, process list_size elements
                                together
      --debug_level 0           Level of debug info in specific debug file
      --debug_file              Generate debug info from inkscape (default:
                                true)
  -k, --original                Output Original layer
  -n, --nested                  Keep nested path together (default: true)
  -y, --layer_output 0          Output internal layers : 1 Input layer, 2
                                         Polygon, 4 Large polygon, 8 Hull layer, 16 Placed
                                Polygon layer, OR these values to output
                                multiple layers
  -a, --angle 90.0              Rotation step
  -r, --free_rot                allow free rotation (default: true)
  -p, --firstpos Position of largest object on the sheet
                                Position of largest object

## Next steps

The program could be improved, both in performance and optimization result. I will work loosely on both topics, but you want to contribute, you are welcome !

## References
As mentionned, this is a well known research topic, I have used to scientific papers when designing this software
1. Waste minimization in irregular stock cutting published in 2014 by Doraid Dalalah, Samir Khrais and Khaled Bataineh. I have used this paper as the main input of this work.
2. Jostle heuristics for the 2D-irregular shapes bin packing problems with free rotation published in 2018 by Ranga P. Abeysooriya, Julia A. Bennell and Antonio  Martinez-Sykora. 

## Libraries and other contributions

- For command line processing I have used cxxopts (https://github.com/jarro2783/cxxopts)
- For SVG processing I have used nanosvg (https://github.com/memononen/nanosvg)
- For convex hull computation I have used a result from tyhe project Nayuki (https://www.nayuki.io/page/convex-hull-algorithm)

Thanks a lot for sharing these contributions !
