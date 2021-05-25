#include "EDLineDetect.h"

#include <algorithm>
#include <vector>

#include "utility.h"

#include <math.h>

/**

  Copied from the various parts of ED library
  to be removed once we can get edge segment numbers for the lines returned from opencv implementation


**/

#define PI 3.14159265358979323846

///----------------------------------------------
/// Fast arctan2 using a lookup table
///
#define MAX_LUT_SIZE 1024
static double LUT[MAX_LUT_SIZE+1];

double myAtan2(double yy, double xx){
  static bool tableInited = false;
  if (!tableInited){
    for (int i=0; i<=MAX_LUT_SIZE; i++){
      LUT[i] = atan((double)i/MAX_LUT_SIZE);
    } //end-for

    tableInited = true;
  } //end-if

  double y = fabs(yy);
  double x = fabs(xx);

#define EPSILON 0.0001
  if (x < EPSILON){
    if (y < EPSILON) return 0.0;
    else return PI/2;
  } // end-if

  bool invert = false;
  if (y > x){
    double t = x;
    x = y;
    y = t;
    invert = true;
  } //end-if

  double ratio = y/x;
  double angle = LUT[(int)(ratio*MAX_LUT_SIZE)];

  if (xx >= 0){
    if (yy >= 0){
      // I. quadrant
      if (invert) angle = PI/2 - angle;

    } else {
      // IV. quadrant
      if (invert == false) angle = PI - angle;
      else                 angle = PI/2 + angle;
    } //end-else

  } else {
    if (yy >= 0){
      /// II. quadrant
      if (invert == false) angle = PI - angle;
      else                 angle = PI/2 + angle;

    } else {
      /// III. quadrant
      if (invert) angle = PI/2 - angle;
    } //end-else
  } //end-else

  return angle;
} // end-myAtan2

///---------------------------------------------------------
/// Fast square root functions. Up to 6% error
///
float fastsqrt(float val)  {
        union
        {
                int tmp;
                float val;
        } u;
        u.val = val;
        u.tmp -= 1<<23; /* Remove last bit so 1.0 gives 1.0 */
        /* tmp is now an approximation to logbase2(val) */
        u.tmp >>= 1; /* divide by 2 */
        u.tmp += 1<<29; /* add 64 to exponent: (e+127)/2 =(e/2)+63, */
        /* that represents (e/2)-64 but want e/2 */
        return u.val;
} //end-fastsqrt


///------------------------------------------------------------------
/// Fast square root functions -- This gives a better approximation: 3.5% error
///
float fastsqrt2(float f){
  int *tmp = (int *)&f;
  (*tmp) = (1<<29) + ((*tmp) >> 1) - (1<<22) - 0x4C000;
//  (*tmp) = (1<<29) + ((*tmp) >> 1) - (1<<22);
  return f;
} //end-fastsqrt2



///------------------------------------------------------------------
/// Fast square root for a double variable
///
double fastsqrt (double y) {
  double x, z, tempf;
  unsigned long *tfptr = ((unsigned long *)&tempf) + 1;

    tempf = y;
    *tfptr = (0xbfcdd90a - *tfptr)>>1; /* estimate of 1/sqrt(y) */
    x =  tempf;
    z =  y*0.5;                        /* hoist out the �/2�    */
    x = (1.5*x) - (x*x)*(x*z);         /* iteration formula     */
    x = (1.5*x) - (x*x)*(x*z);
    x = (1.5*x) - (x*x)*(x*z);
    x = (1.5*x) - (x*x)*(x*z);
    x = (1.5*x) - (x*x)*(x*z);
    return x*y;
} //end-fastsqrt


/* Copyright (C) 1997 by Vesa Karvonen. All rights reserved.
    **
    ** Use freely as long as my copyright is retained.
    */


///-----------------------------------------------
/// Lookup table (LUT) for NFA computation
///
struct NFALUT {
public:
  int *LUT;
  int LUTSize;

  double prob;
  double logNT;

public:
  /// Constructor
  NFALUT(int size, double prob, double logNT);

  // Destructor
  ~NFALUT(){
    delete LUT;
  } //end-~NFALUT
};

///-------------------------------------------
/// Validation check width/without LUT
///
bool checkValidationByNFA(int n, int k, double prob, double logNT);
bool checkValidationByNFA(int n, int k, NFALUT *lut);

/// nfa function prototype
double nfa(int n, int k, double p, double logNT);

///----------------------------------------------
/// Look Up Table (LUT) for NFA Computation
///
NFALUT::NFALUT(int size, double _prob, double _logNT){
  LUTSize = size;
  LUT = new int[LUTSize];

  prob = _prob;
  logNT = _logNT;

  LUT[0] = 1;
  int j = 1;
  for (int i=1; i<LUTSize; i++){
    LUT[i] = LUTSize + 1;
    double ret = nfa(i, j, prob, logNT);
    if (ret < 0){
      while (j < i){
        j++;
        ret = nfa(i, j, prob, logNT);
        if (ret >= 0) break;
      } //end-while

      if (ret < 0) continue;
    } //end-if

    LUT[i] = j;
  } //end-for

#if 0
  fprintf(stderr, "================== ENTIRE TABLE ====================\n");
  for (int i=0; i<LUTSize; i++){
    fprintf(stderr, "n: %4d, k: %4d\n", i, LUT[i]);
  } //end-for
#endif
} //end-LUT

///-------------------------------------------
/// Validation check without the help of a LUT
///
bool checkValidationByNFA(int n, int k, double prob, double logNT){
  return nfa(n, k, prob, logNT) >= 0.0;
} //end-checkValidationByNFA

///-------------------------------------------
/// Validation check with the help of a LUT
///
bool checkValidationByNFA(int n, int k, NFALUT *lut){
  if (n >= lut->LUTSize) return nfa(n, k, lut->prob, lut->logNT) >= 0.0;
  else return k >= lut->LUT[n];
} // end-checkValidationByNFA

///------------------------------------------------------
/// The rest of this code used for NFA computation is directly taken from
/// the LSD code distribution. Hope they do not mind :-)
///
#define error(str) {printf(str); exit(1);}

/** ln(10) */
#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif /* !M_LN10 */

/** PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/** Label for pixels with undefined gradient. */
#define NOTDEF -1024.0

/** 3/2 pi */
#define M_3_2_PI 4.71238898038

/** 2 pi */
#define M_2__PI  6.28318530718

/** Label for pixels not used in yet. */
#define NOTUSED 0

/** Label for pixels already used in detection. */
#define USED    1

#define RELATIVE_ERROR_FACTOR 100.0

static int double_equal(double a, double b){
  double abs_diff,aa,bb,abs_max;

  /* trivial case */
  if( a == b ) return TRUE;

  abs_diff = fabs(a-b);
  aa = fabs(a);
  bb = fabs(b);
  abs_max = aa > bb ? aa : bb;

  /* DBL_MIN is the smallest normalized number, thus, the smallest
     number whose relative error is bounded by DBL_EPSILON. For
     smaller numbers, the same quantization steps as for DBL_MIN
     are used. Then, for smaller numbers, a meaningful "relative"
     error should be computed by dividing the difference by DBL_MIN. */
  if( abs_max < DBL_MIN ) abs_max = DBL_MIN;

  /* equal if relative error <= factor x eps */
  return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}

#define TABSIZE 100000

static double log_gamma_windschitl(double x)
{
  return 0.918938533204673 + (x-0.5)*log(x) - x
         + 0.5*x*log( x*sinh(1/x) + 1/(810.0*pow(x,6.0)) );
}

static double log_gamma_lanczos(double x)
{
  static double q[7] = { 75122.6331530, 80916.6278952, 36308.2951477,
                         8687.24529705, 1168.92649479, 83.8676043424,
                         2.50662827511 };
  double a = (x+0.5) * log(x+5.5) - (x+5.5);
  double b = 0.0;
  int n;

  for(n=0;n<7;n++)
    {
      a -= log( x + (double) n );
      b += q[n] * pow( x, (double) n );
    }
  return a + log(b);
}

#define log_gamma(x) ((x)>15.0?log_gamma_windschitl(x):log_gamma_lanczos(x))

double nfa(int n, int k, double p, double logNT){
  static double inv[TABSIZE];   /* table to keep computed inverse values */
  double tolerance = 0.1;       /* an error of 10% in the result is accepted */
  double log1term,term,bin_term,mult_term,bin_tail,err,p_term;
  int i;

  /* check parameters */
  if( n<0 || k<0 || k>n || p<=0.0 || p>=1.0 ) return -1.0;

  /* trivial cases */
  if( n==0 || k==0 ) return -logNT;
  if( n==k ) return -logNT - (double) n * log10(p);

  /* probability term */
  p_term = p / (1.0-p);

  /* compute the first term of the series */
  /*
     binomial_tail(n,k,p) = sum_{i=k}^n bincoef(n,i) * p^i * (1-p)^{n-i}
     where bincoef(n,i) are the binomial coefficients.
     But
       bincoef(n,k) = gamma(n+1) / ( gamma(k+1) * gamma(n-k+1) ).
     We use this to compute the first term. Actually the log of it.
   */
  log1term = log_gamma( (double) n + 1.0 ) - log_gamma( (double) k + 1.0 )
           - log_gamma( (double) (n-k) + 1.0 )
           + (double) k * log(p) + (double) (n-k) * log(1.0-p);
  term = exp(log1term);

  /* in some cases no more computations are needed */
  if (double_equal(term, 0.0)){              /* the first term is almost zero */
    if( (double) k > (double) n * p )     /* at begin or end of the tail?  */
      return -log1term / M_LN10 - logNT;  /* end: use just the first term  */
    else
      return -logNT;                      /* begin: the tail is roughly 1  */
  } //end-if

  /* compute more terms if needed */
  bin_tail = term;
  for (i=k+1; i<=n; i++){
    /*
       As
         term_i = bincoef(n,i) * p^i * (1-p)^(n-i)
       and
         bincoef(n,i)/bincoef(n,i-1) = n-1+1 / i,
       then,
         term_i / term_i-1 = (n-i+1)/i * p/(1-p)
       and
         term_i = term_i-1 * (n-i+1)/i * p/(1-p).
       1/i is stored in a table as they are computed,
       because divisions are expensive.
       p/(1-p) is computed only once and stored in 'p_term'.
     */
    bin_term = (double) (n-i+1) * ( i<TABSIZE ?
                 ( inv[i]!=0.0 ? inv[i] : ( inv[i] = 1.0 / (double) i ) ) :
                 1.0 / (double) i );

    mult_term = bin_term * p_term;
    term *= mult_term;
    bin_tail += term;

    if (bin_term<1.0){
      /* When bin_term<1 then mult_term_j<mult_term_i for j>i.
         Then, the error on the binomial tail when truncated at
         the i term can be bounded by a geometric series of form
         term_i * sum mult_term_i^j.                            */
      err = term * ( ( 1.0 - pow( mult_term, (double) (n-i+1) ) ) /
                     (1.0-mult_term) - 1.0 );

      /* One wants an error at most of tolerance*final_result, or:
         tolerance * abs(-log10(bin_tail)-logNT).
         Now, the error that can be accepted on bin_tail is
         given by tolerance*final_result divided by the derivative
         of -log10(x) when x=bin_tail. that is:
         tolerance * abs(-log10(bin_tail)-logNT) / (1/bin_tail)
         Finally, we truncate the tail if the error is less than:
         tolerance * abs(-log10(bin_tail)-logNT) * bin_tail        */
      if (err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail) break;
    } //end-if
  } //end-for

  return -log10(bin_tail) - logNT;
} // end-nfa


///------------------------------------------------------------------
/// Rounds a double number to its closest integer part.
/// E.g., 4.24-->4, 4.78-->5
///
int Round(double d){
  return (int)(d+0.5);
} //end-Round

///-----------------------------------------------------------------------------------------
/// Computes the minimum line length using the NFA formula given width & height values
///
int ComputeMinLineLength(int width, int height){
  // The reason we are dividing the theoretical minimum line length by 2 is because
  // we now test short line segments by a line support region rectangle having width=2.
  // This means that within a line support region rectangle for a line segment of length "l"
  // there are "2*l" many pixels. Thus, a line segment of length "l" has a chance of getting
  // validated by NFA.

  double logNT = 2.0*(log10((double)width) + log10((double)height));
  return Round((-logNT/log10(0.125))*0.5);
} //end-ComputeMinLineLength

//-----------------------------------------------------------------------------------
/// Fits a line of the form y=a+bx (invert == 0) OR x=a+by (invert == 1)
///
void LineFit(double *x, double *y, int count, double *a, double *b, double *e, int *invert) {
  if (count<2) return;

  double S=count, Sx=0.0, Sy=0.0, Sxx=0.0, Sxy=0.0;
  for (int i=0; i<count; i++) {
    Sx  += x[i];
    Sy  += y[i];
  } //end-for

  double mx = Sx/count;
  double my = Sy/count;

  double dx = 0.0;
  double dy = 0.0;
  for (int i=0; i < count; i++) {
    dx += (x[i] - mx)*(x[i] - mx);
    dy += (y[i] - my)*(y[i] - my);
  } //end-for

  if (dx < dy) {
    // Vertical line. Swap x & y, Sx & Sy
    *invert = 1;
    double *t = x;
    x = y;
    y = t;

    double d = Sx;
    Sx = Sy;
    Sy = d;

  } else {
    *invert = 0;
  } //end-else

  // Now compute Sxx & Sxy
  for (int i=0; i<count; i++) {
    Sxx += x[i] * x[i];
    Sxy += x[i] * y[i];
  } //end-for

  double D = S*Sxx - Sx*Sx;
  *a = (Sxx*Sy  - Sx*Sxy)/D;
  *b = (S  *Sxy - Sx* Sy)/D;

  if (*b==0.0) {
    // Vertical or horizontal line
    double error = 0.0;
    for (int i=0; i<count; i++) {
      error += fabs((*a) - y[i]);
    } //end-for
    *e = error/count;

  } else {
    double error = 0.0;
    for (int i=0; i<count; i++){
      // Let the line passing through (x[i], y[i]) that is perpendicular to a+bx be c+dx
      double d = -1.0/(*b);
      double c = y[i]-d*x[i];
      double x2 = ((*a)-c)/(d-(*b));
      double y2 = (*a)+(*b)*x2;

      double dist = (x[i]-x2)*(x[i]-x2) + (y[i]-y2)*(y[i]-y2);
      error += dist;
    } //end-for

    *e = sqrt(error/count);
  } //end-else
} // end LineFit

///---------------------------------------------------------------------------------
/// Given a point (x1, y1) and a line equation y=a+bx (invert=0) OR x=a+by (invert=1)
/// Computes the minimum distance of (x1, y1) to the line
///
double ComputeMinDistance(double x1, double y1, double a, double b, int invert){
  double x2, y2;

  if (invert == 0){
    if (b == 0){
      x2 = x1;
      y2 = a;

    } else {
      // Let the line passing through (x1, y1) that is perpendicular to a+bx be c+dx
      double d = -1.0/(b);
      double c = y1-d*x1;

      x2 = (a-c)/(d-b);
      y2 = a+b*x2;
    } //end-else

  } else {
    /// invert = 1
    if (b == 0){
      x2 = a;
      y2 = y1;

    } else {
      // Let the line passing through (x1, y1) that is perpendicular to a+by be c+dy
      double d = -1.0/(b);
      double c = x1-d*y1;

      y2 = (a-c)/(d-b);
      x2 = a+b*y2;
    } //end-else
  } //end-else

  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
} //end-ComputeMinDistance

///-----------------------------------------------------------------------------------
/// Fits a line of the form y=a+bx (invert == 0) OR x=a+by (invert == 1)
/// Assumes that the direction of the line is known by a previous computation
///
void LineFit(double *x, double *y, int count, double *a, double *b, int invert){
  if (count<2) return;

  double S=count, Sx=0.0, Sy=0.0, Sxx=0.0, Sxy=0.0;
  for (int i=0; i<count; i++) {
    Sx  += x[i];
    Sy  += y[i];
  } //end-for

  if (invert) {
    // Vertical line. Swap x & y, Sx & Sy
    double *t = x;
    x = y;
    y = t;

    double d = Sx;
    Sx = Sy;
    Sy = d;
  } //end-if

  // Now compute Sxx & Sxy
  for (int i=0; i<count; i++) {
    Sxx += x[i] * x[i];
    Sxy += x[i] * y[i];
  } //end-for

  double D = S*Sxx - Sx*Sx;
  *a = (Sxx*Sy  - Sx*Sxy)/D;
  *b = (S  *Sxy - Sx* Sy)/D;
} //end-LineFit

///---------------------------------------------------------------------------------
/// Given a point (x1, y1) and a line equation y=a+bx (invert=0) OR x=a+by (invert=1)
/// Computes the (x2, y2) on the line that is closest to (x1, y1)
///
void ComputeClosestPoint(double x1, double y1, double a, double b, int invert, double *xOut, double *yOut){
  double x2, y2;

  if (invert == 0){
    if (b == 0){
      x2 = x1;
      y2 = a;

    } else {
      // Let the line passing through (x1, y1) that is perpendicular to a+bx be c+dx
      double d = -1.0/(b);
      double c = y1-d*x1;

      x2 = (a-c)/(d-b);
      y2 = a+b*x2;
    } //end-else

  } else {
    /// invert = 1
    if (b == 0){
      x2 = a;
      y2 = y1;

    } else {
      // Let the line passing through (x1, y1) that is perpendicular to a+by be c+dy
      double d = -1.0/(b);
      double c = x1-d*y1;

      y2 = (a-c)/(d-b);
      x2 = a+b*y2;
    } //end-else
  } //end-else

  *xOut = x2;
  *yOut = y2;
} //end-ComputeClosestPoint


#define SS 0
#define SE 1
#define ES 2
#define EE 3
///-------------------------------------------------------------------------------
/// Computes the minimum distance between the end points of two lines
///
double ComputeMinDistanceBetweenTwoLines(LineSegment *ls1, LineSegment *ls2, int *pwhich){
  double dx = ls1->sx - ls2->sx;
  double dy = ls1->sy - ls2->sy;
  double d = sqrt(dx*dx + dy*dy);
  double min = d;
  int which = SS;

  dx = ls1->sx - ls2->ex;
  dy = ls1->sy - ls2->ey;
  d = sqrt(dx*dx + dy*dy);
  if (d < min){min = d; which = SE;}

  dx = ls1->ex - ls2->sx;
  dy = ls1->ey - ls2->sy;
  d = sqrt(dx*dx + dy*dy);
  if (d < min){min = d; which = ES;}

  dx = ls1->ex - ls2->ex;
  dy = ls1->ey - ls2->ey;
  d = sqrt(dx*dx + dy*dy);
  if (d < min){min = d; which = EE;}

  if (pwhich) *pwhich = which;
  return min;
} //end-ComputeMinDistanceBetweenTwoLines


///-----------------------------------------------------------------
/// Checks if the given line segments are collinear & joins them if they are
/// In case of a join, ls1 is updated. ls2 is NOT changed
/// Returns true if join is successful, false otherwise
///
bool TryToJoinTwoLineSegments(LineSegment *ls1, LineSegment *ls2, double MAX_DISTANCE_BETWEEN_TWO_LINES, double MAX_ERROR){
  // Must belong to the same segment
//  if (ls1->segmentNo != ls2->segmentNo) return false;

  int which;
  double dist = ComputeMinDistanceBetweenTwoLines(ls1, ls2, &which);
  if (dist > MAX_DISTANCE_BETWEEN_TWO_LINES) return false;

  // Compute line lengths. Use the longer one as the ground truth
  double dx = ls1->sx-ls1->ex;
  double dy = ls1->sy-ls1->ey;
  double prevLen = sqrt(dx*dx + dy*dy);

  dx = ls2->sx-ls2->ex;
  dy = ls2->sy-ls2->ey;
  double nextLen = sqrt(dx*dx + dy*dy);

  // Use the longer line as the ground truth
  LineSegment *shorter = ls1;
  LineSegment *longer = ls2;

  if (prevLen > nextLen){shorter = ls2; longer = ls1;}

#if 0
  // Use 5 points to check for collinearity
#define POINT_COUNT 5
  double decr = 1.0/(POINT_COUNT-1);
  double alpha = 1.0;
  dist = 0.0;

  while (alpha >= 0.0){
    double px = alpha*shorter->sx + (1.0-alpha)*shorter->ex;
    double py = alpha*shorter->sy + (1.0-alpha)*shorter->ey;

    dist += ComputeMinDistance(px, py, longer->a, longer->b, longer->invert);

    alpha -= decr;
  } //end-while

  dist /= POINT_COUNT;

#undef POINT_COUNT

#else
  // Just use 3 points to check for collinearity
  dist = ComputeMinDistance(shorter->sx, shorter->sy, longer->a, longer->b, longer->invert);
  dist += ComputeMinDistance((shorter->sx+shorter->ex)/2.0, (shorter->sy+shorter->ey)/2.0, longer->a, longer->b, longer->invert);
  dist += ComputeMinDistance(shorter->ex, shorter->ey, longer->a, longer->b, longer->invert);

  dist /= 3.0;
#endif

  if (dist > MAX_ERROR) return false;

#if 0
  // Update the end points of ls1
  if (which == 0) {       // SS
    ls1->sx = ls2->ex;
    ls1->sy = ls2->ey;

  } else if (which == 1){ // SE
    ls1->sx = ls2->sx;
    ls1->sy = ls2->sy;

  } else if (which == 2){ // ES
    ls1->ex = ls2->ex;
    ls1->ey = ls2->ey;

  } else {                // EE
    ls1->ex = ls2->sx;
    ls1->ey = ls2->sy;
  } //end-else

#else
  /// 4 cases: 1:(s1, s2), 2:(s1, e2), 3:(e1, s2), 4:(e1, e2)

  /// case 1: (s1, s2)
  dx = fabs(ls1->sx-ls2->sx);
  dy = fabs(ls1->sy-ls2->sy);
  double d = dx+dy;
  double max = d;
  which = 1;

  /// case 2: (s1, e2)
  dx = fabs(ls1->sx-ls2->ex);
  dy = fabs(ls1->sy-ls2->ey);
  d = dx+dy;
  if (d > max){
    max = d;
    which = 2;
  } //end-if

  /// case 3: (e1, s2)
  dx = fabs(ls1->ex-ls2->sx);
  dy = fabs(ls1->ey-ls2->sy);
  d = dx+dy;
  if (d > max){
    max = d;
    which = 3;
  } //end-if

  /// case 4: (e1, e2)
  dx = fabs(ls1->ex-ls2->ex);
  dy = fabs(ls1->ey-ls2->ey);
  d = dx+dy;
  if (d > max){
    max = d;
    which = 4;
  } //end-if

  if (which == 1){
    // (s1, s2)
    ls1->ex = ls2->sx;
    ls1->ey = ls2->sy;

  } else if (which == 2){
    // (s1, e2)
    ls1->ex = ls2->ex;
    ls1->ey = ls2->ey;

  } else if (which == 3){
    // (e1, s2)
    ls1->sx = ls2->sx;
    ls1->sy = ls2->sy;

  } else {
    // (e1, e2)
    ls1->sx = ls1->ex;
    ls1->sy = ls1->ey;

    ls1->ex = ls2->ex;
    ls1->ey = ls2->ey;
  } //end-else

#endif

  // Update the first line's parameters
  if (ls1->firstPixelIndex + ls1->len + 5 >= ls2->firstPixelIndex) ls1->len += ls2->len;
  else if (ls2->len > ls1->len){
    ls1->firstPixelIndex = ls2->firstPixelIndex;
    ls1->len = ls2->len;
  } //end-if

  UpdateLineParameters(ls1);

  return true;
} //end-TryToJoinTwoLineSegments






///-----------------------------------------------------------------
/// Given a full segment of pixels, splits the chain to lines
/// This code is used when we use the whole segment of pixels
///
void SplitSegment2Lines(double *x, double *y, int noPixels, int segmentNo, EDLines *lines){
  double LINE_ERROR = lines->LINE_ERROR;
  int MIN_LINE_LEN = lines->MIN_LINE_LEN;

  // First pixel of the line segment within the segment of points
  int firstPixelIndex = 0;

  while (noPixels >= MIN_LINE_LEN){
    // Start by fitting a line to MIN_LINE_LEN pixels
    bool valid = false;
    double lastA, lastB, error;
    int lastInvert;

    while (noPixels >= MIN_LINE_LEN){
      LineFit(x, y, MIN_LINE_LEN, &lastA, &lastB, &error, &lastInvert);
      if (error <= 0.5){valid = true; break;}

#if 1
      noPixels -= 1;   // Go slowly
      x+= 1; y+= 1;
      firstPixelIndex += 1;
#else
      noPixels -= 2;   // Go faster (for speed)
      x+= 2; y+= 2;
      firstPixelIndex += 2;
#endif
    } //end-while

    if (valid == false) return;

    // Now try to extend this line
    int index = MIN_LINE_LEN;
    int len = MIN_LINE_LEN;

    while (index < noPixels){
      int startIndex = index;
      int lastGoodIndex = index-1;
      int goodPixelCount = 0;
      int badPixelCount = 0;
      while (index < noPixels){

        double d = ComputeMinDistance(x[index], y[index], lastA, lastB, lastInvert);

        if (d <= LINE_ERROR){
          lastGoodIndex = index;
          goodPixelCount++;
          badPixelCount = 0;
        } else {
          badPixelCount++;
          if (badPixelCount >= 5) break;
        } //end-if

        //Burak - Fits a new line every 10 good pixel
          if (goodPixelCount % 10 == 0)
                LineFit(x, y, lastGoodIndex - startIndex + len + 1, &lastA, &lastB, lastInvert);
          //Burak - Fits a new line every 10 good pixel
        index++;
      } //end-while

      if (goodPixelCount >= 2){
        len += lastGoodIndex - startIndex + 1;
        LineFit(x, y, len, &lastA, &lastB, lastInvert);  // faster LineFit
        index = lastGoodIndex+1;
      } // end-if

      if (goodPixelCount < 2 || index >= noPixels){
        // End of a line segment. Compute the end points
        double sx, sy, ex, ey;

        int index = 0;
        while (ComputeMinDistance(x[index], y[index], lastA, lastB, lastInvert) > LINE_ERROR) index++;
        ComputeClosestPoint(x[index], y[index], lastA, lastB, lastInvert, &sx, &sy);
        int noSkippedPixels = index;

        index = lastGoodIndex;
        while (ComputeMinDistance(x[index], y[index], lastA, lastB, lastInvert) > LINE_ERROR) index--;
        ComputeClosestPoint(x[index], y[index], lastA, lastB, lastInvert, &ex, &ey);

        // Add the line segment to lines
        lines->add(lastA, lastB, lastInvert, sx, sy, ex, ey, segmentNo, firstPixelIndex+noSkippedPixels, index-noSkippedPixels+1);

        len = index+1;
        break;
      } //end-else
    } //end-while

    noPixels -= len;
    x += len;
    y += len;
    firstPixelIndex += len;
  } //end-while
} //end-SplitSegment2Lines

///------------------------------------------------------------------
/// Goes over the original line segments and combines collinear lines that belong to the same segment
///
void JoinCollinearLines(EDLines *lines, double MAX_DISTANCE_BETWEEN_TWO_LINES, double MAX_ERROR){
  int lastLineIndex = -1;   //Index of the last line in the joined lines
  int i = 0;
  while (i<lines->noLines){
    int segmentNo = lines->lines[i].segmentNo;

    lastLineIndex++;
    if (lastLineIndex != i) lines->lines[lastLineIndex] = lines->lines[i];  // copy the line
    int firstLineIndex = lastLineIndex;  // Index of the first line in this segment

    int count = 1;
    for (int j=i+1; j<lines->noLines; j++){
      if (lines->lines[j].segmentNo != segmentNo) break;

      // Try to combine this line with the previous line in this segment
      if (TryToJoinTwoLineSegments(&lines->lines[lastLineIndex], &lines->lines[j],
                                   MAX_DISTANCE_BETWEEN_TWO_LINES, MAX_ERROR) == false){
        lastLineIndex++;
        if (lastLineIndex != j) lines->lines[lastLineIndex] = lines->lines[j];  // copy the line
      } //end-if

      count++;
    } //end-for

    // Try to join the first & last line of this segment
    if (firstLineIndex != lastLineIndex){
      if (TryToJoinTwoLineSegments(&lines->lines[firstLineIndex], &lines->lines[lastLineIndex],
                                   MAX_DISTANCE_BETWEEN_TWO_LINES, MAX_ERROR)){
        lastLineIndex--;
      } //end-if
    } //end-if

    i += count;
  } //end-while

  lines->noLines = lastLineIndex+1;
} //end-JoinCollinearLines

//====================================================================================
/// The following rectangle enumeration code was taken from the LSD distribution
/// and adapted for our purposes. We hope that LSD guys are OK with this :-)
///
/// Enumerate the points within a rectangle of width 2 given its two end points
///
void EnumerateRectPoints(double sx, double sy, double ex, double ey,
                         int ptsx[], int ptsy[], int *pNoPoints){

  double vxTmp[4], vyTmp[4];
  double vx[4], vy[4];
  int n, offset;

  double x1 = sx;
  double y1 = sy;
  double x2 = ex;
  double y2 = ey;
  double width = 2;

  double dx = x2 - x1;
  double dy = y2 - y1;
  double vLen = sqrt(dx*dx + dy*dy);

  // make unit vector
  dx = dx/vLen;
  dy = dy/vLen;

  /* build list of rectangle corners ordered
     in a circular way around the rectangle */
  vxTmp[0] = x1 - dy * width / 2.0;
  vyTmp[0] = y1 + dx * width / 2.0;
  vxTmp[1] = x2 - dy * width / 2.0;
  vyTmp[1] = y2 + dx * width / 2.0;
  vxTmp[2] = x2 + dy * width / 2.0;
  vyTmp[2] = y2 - dx * width / 2.0;
  vxTmp[3] = x1 + dy * width / 2.0;
  vyTmp[3] = y1 - dx * width / 2.0;

  /* compute rotation of index of corners needed so that the first
     point has the smaller x.

     if one side is vertical, thus two corners have the same smaller x
     value, the one with the largest y value is selected as the first.
   */
  if      (x1 < x2 && y1 <= y2) offset = 0;
  else if (x1 >= x2 && y1 < y2) offset = 1;
  else if (x1 > x2 && y1 >= y2) offset = 2;
  else                          offset = 3;

  /* apply rotation of index. */
  for (n=0; n<4; n++){
    vx[n] = vxTmp[(offset+n)%4];
    vy[n] = vyTmp[(offset+n)%4];
  } //end-for

  /* Set a initial condition.

     The values are set to values that will cause 'ri_inc' (that will
     be called immediately) to initialize correctly the first 'column'
     and compute the limits 'ys' and 'ye'.

     'y' is set to the integer value of vy[0], the starting corner.

     'ys' and 'ye' are set to very small values, so 'ri_inc' will
     notice that it needs to start a new 'column'.

     The smaller integer coordinate inside of the rectangle is
     'ceil(vx[0])'. The current 'x' value is set to that value minus
     one, so 'ri_inc' (that will increase x by one) will advance to
     the first 'column'.
   */
  int x = (int) ceil(vx[0]) - 1;
  int y = (int) ceil(vy[0]);
  double ys = -DBL_MAX, ye = -DBL_MAX;

  int noPoints = 0;
  int maxNoOfPoints = (int)(fabs(sx-ex) + fabs(sy-ey))*4;
  while (noPoints < maxNoOfPoints){
//  while (1){
    /* if not at end of exploration,
       increase y value for next pixel in the 'column' */
    y++;

    /* if the end of the current 'column' is reached,
       and it is not the end of exploration,
       advance to the next 'column' */
    while (y > ye && x <= vx[2]){
      /* increase x, next 'column' */
      x++;

      /* if end of exploration, return */
      if (x > vx[2]) break;

      /* update lower y limit (start) for the new 'column'.

         We need to interpolate the y value that corresponds to the
         lower side of the rectangle. The first thing is to decide if
         the corresponding side is

           vx[0],vy[0] to vx[3],vy[3] or
           vx[3],vy[3] to vx[2],vy[2]

         Then, the side is interpolated for the x value of the
         'column'. But, if the side is vertical (as it could happen if
         the rectangle is vertical and we are dealing with the first
         or last 'columns') then we pick the lower value of the side
         by using 'inter_low'.
       */
      if ((double) x < vx[3]){
        /* interpolation */
        if (fabs(vx[0]-vx[3])<=0.01){
          if      (vy[0]<vy[3] ) ys = vy[0];
          else if (vy[0]>vy[3] ) ys = vy[3];
          else     ys = vy[0] + (x-vx[0]) * (vy[3]-vy[0]) / (vx[3]-vx[0]);
        } else
          ys = vy[0] + (x-vx[0]) * (vy[3]-vy[0]) / (vx[3]-vx[0]);

      } else {
        /* interpolation */
        if (fabs(vx[3]-vx[2])<=0.01){
          if      (vy[3]<vy[2] ) ys = vy[3];
          else if (vy[3]>vy[2] ) ys = vy[2];
          else     ys = vy[3] + (x-vx[3]) * (y2-vy[3]) / (vx[2]-vx[3]);
        } else
          ys = vy[3] + (x-vx[3]) * (vy[2]-vy[3]) / (vx[2]-vx[3]);
      } //end-else

      /* update upper y limit (end) for the new 'column'.

         We need to interpolate the y value that corresponds to the
         upper side of the rectangle. The first thing is to decide if
         the corresponding side is

           vx[0],vy[0] to vx[1],vy[1] or
           vx[1],vy[1] to vx[2],vy[2]

         Then, the side is interpolated for the x value of the
         'column'. But, if the side is vertical (as it could happen if
         the rectangle is vertical and we are dealing with the first
         or last 'columns') then we pick the lower value of the side
         by using 'inter_low'.
       */
      if ((double)x < vx[1] ){
        /* interpolation */
        if (fabs(vx[0]-vx[1])<=0.01){
          if      (vy[0]<vy[1] ) ye = vy[1];
          else if (vy[0]>vy[1] ) ye = vy[0];
          else     ye = vy[0] + (x-vx[0]) * (vy[1]-vy[0]) / (vx[1]-vx[0]);
        } else
          ye = vy[0] + (x-vx[0]) * (vy[1]-vy[0]) / (vx[1]-vx[0]);

      } else {
        /* interpolation */
        if (fabs(vx[1]-vx[2])<=0.01){
          if      (vy[1]<vy[2] ) ye = vy[2];
          else if (vy[1]>vy[2] ) ye = vy[1];
          else     ye = vy[1] + (x-vx[1]) * (vy[2]-vy[1]) / (vx[2]-vx[1]);
        } else
          ye = vy[1] + (x-vx[1]) * (vy[2]-vy[1]) / (vx[2]-vx[1]);
      } //end-else

      /* new y */
      y = (int)ceil(ys);
    } //end-while

    // Are we done?
    if (x > vx[2]) break;

    ptsx[noPoints] = x;
    ptsy[noPoints] = y;
    noPoints++;
  } //end-while

  *pNoPoints = noPoints;
} //end-EnumerateRectPoints

///----------------------------------------------------------------------
/// Using the rectangle code from LSD to validate a line segment
/// The idea is to form a rectangle of with 2 over the line segment,
/// and compute the angles of all pixels within the rectangle to validate
/// a line segment. Long line segments can easily be validated by just going
/// over a single chain of pixels, but short line segments must be validates
/// using a rectangle of width 2 because some short valid lines are rejected otherwise
/// The above pixel enumeration code is directly taken from LSD code.
/// Hope they are OK with it :-)
///
bool ValidateLineSegmentRect(unsigned char *srcImg, int width, int height, int *x, int *y, LineSegment *ls, NFALUT *LUT){

#define PI 3.14159265358979323846
#define PRECISON_ANGLE 22.5
  static double prec = (PRECISON_ANGLE/180)*PI;
  static double prob = 0.125;
#undef PRECISON_ANGLE

  // Compute Line's angle
  double lineAngle;

  if (ls->invert == 0){
    // y = a + bx
    lineAngle = atan(ls->b);

  } else {
    // x = a + by
    lineAngle = atan(1.0/ls->b);
  } //end-else

  if (lineAngle < 0) lineAngle += PI;

  int noPoints = 0;

  // Enumerate all pixels that fall within the bounding rectangle
  EnumerateRectPoints(ls->sx, ls->sy, ls->ex, ls->ey, x, y, &noPoints);

  int count = 0;
  int aligned = 0;

  for (int i=0; i<noPoints; i++){
    int r = y[i];
    int c = x[i];

    if (r<=0 || r>=height-1 || c<=0 || c>=width-1) continue;

    count++;

    // compute gx & gy using the simple [-1 -1 -1]
    //                                  [ 1  1  1]  filter in both directions
    // Faster method below
    // A B C
    // D x E
    // F G H
    // gx = (C-A) + (E-D) + (H-F)
    // gy = (F-A) + (G-B) + (H-C)
    //
    // To make this faster:
    // com1 = (H-A)
    // com2 = (C-F)
    // Then: gx = com1 + com2 + (E-D) = (H-A) + (C-F) + (E-D) = (C-A) + (E-D) + (H-F)
    //       gy = com2 - com1 + (G-B) = (H-A) - (C-F) + (G-B) = (F-A) + (G-B) + (H-C)
    //
    int com1 = srcImg[(r+1)*width+c+1] - srcImg[(r-1)*width+c-1];
    int com2 = srcImg[(r-1)*width+c+1] - srcImg[(r+1)*width+c-1];

    int gx = com1 + com2 + srcImg[r*width+c+1] - srcImg[r*width+c-1];
    int gy = com1 - com2 + srcImg[(r+1)*width+c] - srcImg[(r-1)*width+c];

    double pixelAngle = myAtan2((double)gx, (double)-gy);
    double diff = fabs(lineAngle - pixelAngle);

    if (diff <= prec || diff >= PI-prec) aligned++;
  } //end-for

  return checkValidationByNFA(count, aligned, LUT);
#undef PI
} //end-ValidateLineSegmentRect

///---------------------------------------------------------------------
/// Given a set of lines, validates them. This is the fast code
/// This code checks the directions of ONLY the pixels on the line
///
void ValidateLineSegments(EdgeMap *map, unsigned char *srcImg, EDLines *lines, EDLines *invalidLines){

  int width = map->width;
  int height = map->height;

  int *x = new int[(width+height)*4];
  int *y = new int[(width+height)*4];

#define PRECISON_ANGLE 22.5
  double prec = (PRECISON_ANGLE/180)*M_PI;
  double prob = 0.125;
#undef PRECISON_ANGLE

  double logNT = 2.0*(log10((double)width) + log10((double)height));

  /// Compute LUT for NFA computation
  ///
  NFALUT *LUT = new NFALUT((width+height)/8, prob, logNT);

  if (invalidLines) invalidLines->clear();

  int noValidLines = 0;
  for (int i=0; i<lines->noLines; i++){
    LineSegment *ls = &lines->lines[i];

    // Compute Line's angle
    double lineAngle;

    if (ls->invert == 0){
      // y = a + bx
      lineAngle = atan(ls->b);

    } else {
      // x = a + by
      lineAngle = atan(1.0/ls->b);
    } //end-else

    if (lineAngle < 0) lineAngle += M_PI;

    Pixel *pixels = &map->segments[ls->segmentNo].pixels[ls->firstPixelIndex];
    int noPixels = ls->len;

    bool valid = false;

    // Accept very long lines without testing. They are almost never invalidated.
    if (ls->len >= 80){
      valid = true;

    // Validate short line segments by a line support region rectangle having width=2
   } else if (ls->len <= 25){
      valid = ValidateLineSegmentRect(srcImg, width, height, x, y, ls, LUT);

    } else {
      // Longer line segments are first validated by a line support region rectangle having width=1 (for speed)
      // If the line segment is still invalid, then a line support region rectangle having width=2 is tried
      // If the line segment fails both tests, it is discarded
      int aligned = 0;
      int count = 0;
      for(int j=0; j<noPixels; j++){
        int r = pixels[j].r;
        int c = pixels[j].c;

        if (r<=0 || r>=height-1 || c<=0 || c>=width-1) continue;

        count++;

        // compute gx & gy using the simple [-1 -1 -1]
        //                                  [ 1  1  1]  filter in both directions
        // Faster method below
        // A B C
        // D x E
        // F G H
        // gx = (C-A) + (E-D) + (H-F)
        // gy = (F-A) + (G-B) + (H-C)
        //
        // To make this faster:
        // com1 = (H-A)
        // com2 = (C-F)
        // Then: gx = com1 + com2 + (E-D) = (H-A) + (C-F) + (E-D) = (C-A) + (E-D) + (H-F)
        //       gy = com2 - com1 + (G-B) = (H-A) - (C-F) + (G-B) = (F-A) + (G-B) + (H-C)
        //
        int com1 = srcImg[(r+1)*width+c+1] - srcImg[(r-1)*width+c-1];
        int com2 = srcImg[(r-1)*width+c+1] - srcImg[(r+1)*width+c-1];

        int gx = com1 + com2 + srcImg[r*width+c+1] - srcImg[r*width+c-1];
        int gy = com1 - com2 + srcImg[(r+1)*width+c] - srcImg[(r-1)*width+c];


        double pixelAngle = myAtan2((double)gx, (double)-gy);
        double diff = fabs(lineAngle - pixelAngle);

        if (diff <= prec || diff >= M_PI-prec) aligned++;
      } //end-for

      // Check validation by NFA computation (fast due to LUT)
      valid = checkValidationByNFA(count, aligned, LUT);
      if (valid == false) valid = ValidateLineSegmentRect(srcImg, width, height, x, y, ls, LUT);
    } //end-else

    if (valid){
      if (i != noValidLines) lines->lines[noValidLines] = lines->lines[i];
      noValidLines++;

    } else if (invalidLines){
      invalidLines->add(&lines->lines[i]);
    } //end-else
  } //end-for

  lines->noLines = noValidLines;

  delete LUT;

  delete x;
  delete y;

} // end-ValidateLineSegments


EDLines* detectLines(unsigned char *srcImg, EdgeMap* map, int width, int height)
{
    EDLines *lines = new EDLines(width, height);


    /*----------- FIT LINES ----------------*/

    // Now, go over the edge segments & fit lines
    lines->clear();

    lines->MIN_LINE_LEN = ComputeMinLineLength(width, height);

    // Too small?
  //  if (lines->MIN_LINE_LEN < 8) lines->MIN_LINE_LEN = 8;

    // Min line length 9 seems to give the closest results to LSD. Keep this?
    if (lines->MIN_LINE_LEN < 9) lines->MIN_LINE_LEN = 9;

    // A min line length of 10 gives less lines than LSD, but the lines are much cleaner
  //  if (lines->MIN_LINE_LEN < 10) lines->MIN_LINE_LEN = 10;

    double *x = lines->x;
    double *y = lines->y;

    // Use the whole segment
    for (int segmentNo=0; segmentNo<map->noSegments; segmentNo++){
      EdgeSegment *segment = &map->segments[segmentNo];

      for (int k=0; k<segment->noPixels; k++){
        x[k]=segment->pixels[k].c;
        y[k]=segment->pixels[k].r;
      } //end-for

      SplitSegment2Lines(x, y, segment->noPixels, segmentNo, lines);
    } //end-for


    EDLines *invalidLines = NULL;

    JoinCollinearLines(lines, 6.0, 1.50);
    ValidateLineSegments(map, srcImg, lines, invalidLines);

    return lines;
}
