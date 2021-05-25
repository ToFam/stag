#ifndef EDTYPES_H
#define EDTYPES_H

#include <memory.h>

enum GradientOperator {PREWITT_OPERATOR=101, SOBEL_OPERATOR=102, SCHARR_OPERATOR=103};

struct Pixel {int r, c;};

struct EdgeSegment {
  Pixel *pixels;       // Pointer to the pixels array
  int noPixels;        // # of pixels in the edge map
};

struct LineSegment {
  double a, b;          // y = a + bx (if invert = 0) || x = a + by (if invert = 1)
  int invert;

  double sx, sy;        // starting x & y coordinates
  double ex, ey;        // ending x & y coordinates

  int segmentNo;        // Edge segment that this line belongs to
  int firstPixelIndex;  // Index of the first pixel within the segment of pixels
  int len;              // No of pixels making up the line segment
};

///----------------------------------------------
/// Simple class to manipulate line segments
///
struct EDLines {
public:
  LineSegment *lines;
  int noLines;
  int capacity;

  // Thresholds used during line segment computation
  double LINE_ERROR;       // Mean square error during line fitting
  int MIN_LINE_LEN;        // Minimum length of the line segments

  // Temporary buffers used during line fitting
  double *x;
  double *y;

  // Used for timing :-)
  double edgeDetectionTime;
  double lineFitTime;
  double joinLineSegmentsTime;
  double lineValidationTime;
  double LUTComputationTime;      // Look Up Table (LUT) computation for line validation

public:
  /// Constructor:
  EDLines(int width, int height){
    int imageSize = width+height;
    capacity = imageSize*8;
    lines = new LineSegment[capacity];
    noLines = 0;

    LINE_ERROR = 1.0;   // in pixels
    MIN_LINE_LEN = 11;  // in pixels

    edgeDetectionTime = lineFitTime = joinLineSegmentsTime = 0;
    lineValidationTime = LUTComputationTime = 0;

    x = new double[imageSize*8];
    y = new double[imageSize*8];
  } //end-EDLines

  /// Destructor
  ~EDLines(){
    delete lines;
    delete x;
    delete y;
  } //end-EDLines

  /// clear
  void clear(){noLines = 0;}

  void expandCapacity(){
    capacity *= 2;
    LineSegment *newArr = new LineSegment[capacity];
    memcpy(newArr, lines, sizeof(LineSegment)*noLines);
    delete lines;
    lines = newArr;
  } //end-expandCapacity

  /// Append a new line to the end of the lines array
  void add(double a, double b, int invert, double sx, double sy, double ex, double ey, int segmentNo=0, int firstPixelIndex=0, int len=0){
    // Full array?
    if (noLines == capacity) expandCapacity();

    lines[noLines].a = a;
    lines[noLines].b = b;
    lines[noLines].invert = invert;
    lines[noLines].sx = sx;
    lines[noLines].sy = sy;
    lines[noLines].ex = ex;
    lines[noLines].ey = ey;

    lines[noLines].segmentNo = segmentNo;
    lines[noLines].firstPixelIndex = firstPixelIndex;
    lines[noLines].len = len;

    noLines++;
  } //end-add

  void add(LineSegment *ls){
    // Full array?
    if (noLines == capacity) expandCapacity();

    // Copy
    lines[noLines] = *ls;

    noLines++;
  } //end-add
};

struct EdgeMap {
public:
  int width, height;        // Width & height of the image
  unsigned char *edgeImg;   // BW edge map


  Pixel *pixels;            // Edge map in edge segment form
  EdgeSegment *segments;
  int noSegments;

public:
  // constructor
  EdgeMap(int w, int h){
    width = w;
    height = h;

    edgeImg = new unsigned char[width*height];

    pixels = new Pixel[width*height];
    segments = new EdgeSegment[width*height];
    noSegments = 0;
  } //end-EdgeMap

  // Destructor
  ~EdgeMap(){
    delete edgeImg;
    delete pixels;
    delete segments;
  } //end-~EdgeMap


  void ConvertEdgeSegments2EdgeImg(){
    memset(edgeImg, 0, width*height);

    for (int i=0; i<noSegments; i++){
      for (int j=0; j<segments[i].noPixels; j++){
        int r = segments[i].pixels[j].r;
        int c = segments[i].pixels[j].c;

        edgeImg[r*width+c] = 255;
      } //end-for
    } //end-for
  } //end-ConvertEdgeSegments2EdgeImg

  EdgeMap *clone(){
    EdgeMap *map2 = new EdgeMap(width, height);
    map2->noSegments = noSegments;

    Pixel *pix = map2->pixels;

    for (int i=0; i<noSegments; i++){
      map2->segments[i].noPixels = segments[i].noPixels;
      map2->segments[i].pixels = pix;
      pix += segments[i].noPixels;

      for (int j=0; j<segments[i].noPixels; j++){
        map2->segments[i].pixels[j] = segments[i].pixels[j];
      } //end-for
    } //end-for

    return map2;
  } //end-clone
};

///-----------------------------------------------------------------------------------
/// Uses the two end points (sx, sy)----(ex, ey) of the line segment
/// and computes the line that passes through these points (a, b, invert)
///
void UpdateLineParameters(LineSegment *ls);

#endif
