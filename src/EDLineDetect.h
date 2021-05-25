#ifndef ED_LINE_DETECT_H
#define ED_LINE_DETECT_H

#include "EDTypes.h"


EDLines* detectLines(unsigned char *srcImg, EdgeMap* map, int width, int height);

#endif
