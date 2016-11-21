/* ========================================================================
    File: @(#)Gray2Vec_Grid.h
   ------------------------------------------------------------------------
    Main class for grayscale image vectorizer
    Copyright (C) 2016 Christoph Hormann <chris_hormann@gmx.de>
   ------------------------------------------------------------------------

    This file is part of gray2vec

    gray2vec is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gray2vec is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with gray2vec.  If not, see <http://www.gnu.org/licenses/>.

    Version history:

      0.1: initial public version, November 2016, 

   ========================================================================
 */

#ifndef _Gray2Vec_Grid_H
#define _Gray2Vec_Grid_H

#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#define cimg_display 0

#include "gdal_polygonize_mod.h"

#include "CImg.h"

using namespace cimg_library;

const static int x4[4] = { -1,0,1,0 };
const static int y4[4] = { 0,-1,0,1 };

class Gray2Vec_Grid
{
 public:
	Gray2Vec_Grid(const std::string file, const std::string file_c, const bool Complement, const bool Debug);
	~Gray2Vec_Grid() { };

	/// decides on the pixel classes based on the coverage fraction and the subpixel configuration
	void Analyze();
	/// adjust pixel classes based on neighbors
	void NeighborsAdjust();
	/// adjust pixel classes based on neighbors (part 2: switching corners to sides after tuning)
	void NeighborsAdjust2();
	/// resolve orientation conflicts between neighbor pixels
	void ResolveConflicts();
	/// initialize side fraction based on pixel coverage fractions
	void InitFractions();
	/// tune fractions after FractionsNeighborsAdj()
	void TuneFractions(const double max_error = -1.0);
	/// adjust fractions of matching neighbors
	void FractionsNeighborsAdj();
	/// vectorize the data and write polygons to an SQLite database
	bool Vectorize(const std::string file, const std::string layer, const bool Append);
	/// set x/y/z attributes to be written with the vector data
	void SetAttributes(const int Xc, const int Yc, const int Zc) { m_x = Xc; m_y = Yc; m_z = Zc; };

 protected:
	/// check if neighbourhood n covers direction d
	static bool check_cover(int n, int d);
	static void move_dir(int &px, int &py, const int dir);
	static int side1(const int dir);
	static int side2(const int dir);

	int share_sides(const int px1, const int py1, const int px2, const int py2);
	int sides_connected(const int px, const int py);
	int set_fraction(const int px, const int py);
	double pixel_error(const int px, const int py, const bool use_adjust);

	/// write out a polygon feature to the specified OGR layer
	bool EmitPolygonToLayer(OGRLayerH hOutLayer, RPolygon *poRPoly);
	/// vectorize the processed data using the helper grid img_h
	bool Polygonize(CImg<unsigned char> &img_h, OGRLayerH hOutLayer);

	double m_GeoTransform[6];
	OGRSpatialReferenceH m_SRS;

	bool m_debug;

	CImg<unsigned char> m_img;
	CImg<unsigned char> m_img_s;
	CImg<unsigned char> m_img_n;
	CImg<unsigned char> m_img_h;
	CImg<unsigned char> m_img_f1;
	CImg<unsigned char> m_img_f2;
	CImg<short> m_img_f3;

	int m_x;
	int m_y;
	int m_z;
};

#endif /* _Gray2Vec_Grid_H */
