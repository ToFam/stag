#include "EDTypes.h"

#include <math.h>


///-----------------------------------------------------------------------------------
/// Uses the two end points (sx, sy)----(ex, ey) of the line segment
/// and computes the line that passes through these points (a, b, invert)
///
void UpdateLineParameters(LineSegment *ls){
  double dx = ls->ex-ls->sx;
  double dy = ls->ey-ls->sy;

  if (fabs(dx) >= fabs(dy)){
    /// Line will be of the form y = a + bx
    ls->invert = 0;
    if (fabs(dy) < 1e-3){ls->b = 0; ls->a = (ls->sy+ls->ey)/2;}
    else {
      ls->b = dy/dx;
      ls->a = ls->sy-(ls->b)*ls->sx;
    } //end-else

  } else {
    /// Line will be of the form x = a + by
    ls->invert = 1;
    if (fabs(dx) < 1e-3){ls->b = 0; ls->a = (ls->sx+ls->ex)/2;}
    else {
      ls->b = dx/dy;
      ls->a = ls->sx-(ls->b)*ls->sy;
    } //end-else
  } //end-else
} //end-UpdateLineParameters

