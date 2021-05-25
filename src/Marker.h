#ifndef MARKER_H
#define MARKER_H

#include "Quad.h"

#include "libstag_export.h"

class LIBSTAG_EXPORT Marker : public Quad
{
public:
	int id;
	cv::Mat C;

	Marker(const Quad &q, int inId);
	Marker(const Marker &m);
	void shiftCorners2(int shift);
};

#endif
