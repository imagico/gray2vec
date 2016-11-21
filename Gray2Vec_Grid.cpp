/* ========================================================================
    File: @(#)Gray2Vec_Grid.cpp
   ------------------------------------------------------------------------
    Main class for grayscale image vectorizer
    Copyright (C) 2016 Christoph Hormann <chris_hormann@gmx.de>

    contains code

    Copyright (c) 2008, Frank Warmerdam
    Copyright (c) 2009-2011, Even Rouault <even dot rouault at mines-paris dot org>

    see specific comments in the code for more details.
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

// meaning of neighbourhood directions:
//  1 2 3
//  8 0 4
//  7 6 5
// 
// meaning of neighbourhood codes:
// 
// 1:
// 
//  -------
// |xxxx   |
// |xx     |
// |x      |
// |       |
//  -------
// 
// 2:
// 
//  -------
// |xxxxxxx|
// |xxxxxxx|
// |       |
// |       |
//  -------
// 
// 3:
// 
//  -------
// |   xxxx|
// |     xx|
// |      x|
// |       |
//  -------
// 
// .........
// 
// 6:
// 
//  -------
// |       |
// |       |
// |xxxxxxx|
// |xxxxxxx|
//  -------
// 
// .........
// 
// 11:
// 
//  -------
// |    xxx|
// |  xxxxx|
// |xxxxxxx|
// |xxxxxxx|
//  -------
//
// 
// 13:
// 
//  -------
// |xxx    |
// |xxxxx  |
// |xxxxxxx|
// |xxxxxxx|
//  -------
// 
// .........
//


#include <cmath>
#include <cstdlib>

#include "Gray2Vec_Grid.h"


bool Gray2Vec_Grid::check_cover(int n, int d)
{
	int n2;

	if (n == 0) return false;
	if (n == 255) return true;

	// check corners
	if ((n % 2) != 0)
	{
		// same corner and two neighbour sides
		if (d == n) return true;
		if (n > 10) if (n-10 == d) return true;
		n2 = n-1;
		if (n2 >= 10) n2 -= 10;
		if (n2 > 0)
			if (d == n2) return true;

		n2 = n+1;
		if (n2 > 10) n2 -= 10;
		if (n2 < 9)
			if (d == n2) return true;

		// additional corners for large coverage directions
		if (n > 10)
		{
			n2 = n-10-2;
			if (n2 < 0) n2 += 8;
			if (d == n2) return true;

			n2 = n-10+2;
			if (n2 > 8) n2 -= 8;
			if (d == n2) return true;
		}
		return false;
	}
	else // check sides
	{
		// same side and two neighbour corners
		if (d == n) return true;
		n2 = n-1;
		if (d == n2) return true;
		n2 = n+1;
		if (n2 > 8) n2 -= 8;
		if (d == n2) return true;
		return false;
	}
}

Gray2Vec_Grid::Gray2Vec_Grid(const std::string file, const std::string file_c, const bool Complement, const bool Debug)
	: m_debug(Debug)
{
	GDALAllRegister();
	OGRRegisterAll();

	m_x = -1;
	m_y = -1;
	m_z = -1;

	std::fprintf(stderr,"Loading image data...\n");

	GDALDataset  *poDataset;

	poDataset = static_cast<GDALDataset *>(GDALOpen( file.c_str(), GA_ReadOnly ));

	if (poDataset == NULL)
	{
		std::fprintf(stderr,"  opening file %s failed.\n\n", file.c_str());
		std::exit(1);
	}

	if (poDataset->GetGeoTransform( m_GeoTransform ) != CE_None)
	{
		std::fprintf(stderr,"  error reading coordinates from file %s.\n\n", file.c_str());
		std::exit(1);
	}

	m_SRS = OSRNewSpatialReference( poDataset->GetProjectionRef() );

	GDALRasterBand  *poBand;
	poBand = poDataset->GetRasterBand(1);

	int nXSize = poBand->GetXSize();
	int nYSize = poBand->GetYSize();

	m_img = CImg<unsigned char>(nXSize, nYSize, 1, 1);

	if (poBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize, m_img.data(), nXSize, nYSize, GDT_Byte, 0, 0) != CE_None)
	{
		std::fprintf(stderr,"  error reading image data from file %s.\n\n", file.c_str());
		std::exit(1);
	}

	std::fprintf(stderr,"input image: %s:\n", file.c_str());
	std::fprintf(stderr,"  %d x %d pixel\n", nXSize, nYSize);
	std::fprintf(stderr,"  (%d x %d pixel reduced)\n", m_img.width()/2, m_img.height()/2);
	
	if (m_debug)
	{
		std::fprintf(stderr,"coordinates:\n");
		float cx = m_GeoTransform[0] + m_GeoTransform[1] * 0 + m_GeoTransform[2] * 0;
		float cy = m_GeoTransform[3] + m_GeoTransform[4] * 0 + m_GeoTransform[5] * 0;
		std::fprintf(stderr," corner 1: %f/%f\n", cx, cy);
		cx = m_GeoTransform[0] + m_GeoTransform[1] * nXSize + m_GeoTransform[2] * 0;
		cy = m_GeoTransform[3] + m_GeoTransform[4] * nXSize + m_GeoTransform[5] * 0;
		std::fprintf(stderr," corner 2: %f/%f\n", cx, cy);
		cx = m_GeoTransform[0] + m_GeoTransform[1] * 0 + m_GeoTransform[2] * nYSize;
		cy = m_GeoTransform[3] + m_GeoTransform[4] * 0 + m_GeoTransform[5] * nYSize;
		std::fprintf(stderr," corner 3: %f/%f\n", cx, cy);
		cx = m_GeoTransform[0] + m_GeoTransform[1] * nXSize + m_GeoTransform[2] * nYSize;
		cy = m_GeoTransform[3] + m_GeoTransform[4] * nXSize + m_GeoTransform[5] * nYSize;
		std::fprintf(stderr," corner 4: %f/%f\n", cx, cy);
	}

	m_img_n = CImg<unsigned char>(m_img.width()/2, m_img.height()/2, 1, 1);
	m_img_s = CImg<unsigned char>(m_img.width()/2, m_img.height()/2, 1, 1);

	std::fprintf(stderr,"Averaging values...\n");

	cimg_forXY(m_img_s,px,py)
	{
		m_img_s(px,py) = 0.25*(m_img(px*2,py*2)+m_img(px*2+1,py*2)+m_img(px*2,py*2+1)+m_img(px*2+1,py*2+1));
	}

	if (!file_c.empty())
	{
		if (m_debug)  m_img_s.save("debug-so.tif");

		std::fprintf(stderr,"Loading combined image data...\n");

		GDALDataset  *poDataset2;

		poDataset2 = static_cast<GDALDataset *>(GDALOpen( file_c.c_str(), GA_ReadOnly ));

		if (poDataset2 == NULL)
		{
			std::fprintf(stderr,"  opening file %s failed.\n\n", file_c.c_str());
			std::exit(1);
		}

		GDALRasterBand  *poBand2;
		poBand2 = poDataset2->GetRasterBand(1);

		CImg<unsigned char> img_c = CImg<unsigned char>(nXSize, nYSize, 1, 1);

		if (poBand2->RasterIO(GF_Read, 0, 0, nXSize, nYSize, img_c.data(), nXSize, nYSize, GDT_Byte, 0, 0) != CE_None)
		{
			std::fprintf(stderr,"  error reading image data from file %s.\n\n", file_c.c_str());
			std::exit(1);
		}

		std::fprintf(stderr,"Processing partical pixels...\n");

		cimg_forXY(m_img_s,px,py)
		{
			if (Complement)
				m_img_s(px,py) = 0.25*(img_c(px*2,py*2)+img_c(px*2+1,py*2)+img_c(px*2,py*2+1)+img_c(px*2+1,py*2+1)) - m_img_s(px,py);

			// partial pixels are set to a value that - in combination with the rest
			// of the combined data approximate the target background value to avoid 
			// the background shining through with AGG type renderers
			if (m_img_s(px,py) != 0)
				if (m_img_s(px,py) != 255)
				{
					int fc = 0.25*(img_c(px*2,py*2)+img_c(px*2+1,py*2)+img_c(px*2,py*2+1)+img_c(px*2+1,py*2+1));
					if (fc > m_img_s(px,py))
						m_img_s(px,py) = 255*(1.0 - (1.0-fc/255.0)/(1.0-(fc-m_img_s(px,py))/255.0));
				}
		}
	}
}

void Gray2Vec_Grid::Analyze()
{
	std::fprintf(stderr,"Determining sides...\n");

	cimg_forXY(m_img_s,px,py)
	{
		m_img_n(px,py) = 0;

		if (m_img_s(px,py) != 0)
			if (m_img_s(px,py) == 255)
				m_img_n(px,py) = 255;
			else
			{
				// small fractions: use highest value corner
				if (m_img_s(px,py) < 255/3)
				{
					if (m_img(px*2,py*2) > m_img(px*2+1,py*2))
					{
						if (m_img(px*2,py*2) > m_img(px*2,py*2+1))
						{
							if (m_img(px*2,py*2) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 1;
							else
								m_img_n(px,py) = 5;
						}
						else
						{
							if (m_img(px*2,py*2+1) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 7;
							else
								m_img_n(px,py) = 5;
						}
					}
					else
					{
						if (m_img(px*2+1,py*2) > m_img(px*2,py*2+1))
						{
							if (m_img(px*2+1,py*2) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 3;
							else
								m_img_n(px,py) = 5;
						}
						else
						{
							if (m_img(px*2,py*2+1) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 7;
							else
								m_img_n(px,py) = 5;
						}
					}
				}
				// big fractions: use highest value corner
				else if (m_img_s(px,py) > 2*255/3)
				{
					if (m_img(px*2,py*2) > m_img(px*2+1,py*2))
					{
						if (m_img(px*2,py*2) > m_img(px*2,py*2+1))
						{
							if (m_img(px*2,py*2) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 11;
							else
								m_img_n(px,py) = 15;
						}
						else
						{
							if (m_img(px*2,py*2+1) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 17;
							else
								m_img_n(px,py) = 15;
						}
					}
					else
					{
						if (m_img(px*2+1,py*2) > m_img(px*2,py*2+1))
						{
							if (m_img(px*2+1,py*2) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 13;
							else
								m_img_n(px,py) = 15;
						}
						else
						{
							if (m_img(px*2,py*2+1) > m_img(px*2+1,py*2+1))
								m_img_n(px,py) = 17;
							else
								m_img_n(px,py) = 15;
						}
					}
				}
				// medium fractions: use highest value side
				else if (m_img_s(px,py) <= 2*255/3)
				{
					int s2 = m_img(px*2,py*2) + m_img(px*2+1,py*2);
					int s4 = m_img(px*2+1,py*2+1) + m_img(px*2+1,py*2);
					int s6 = m_img(px*2+1,py*2+1) + m_img(px*2,py*2+1);
					int s8 = m_img(px*2,py*2) + m_img(px*2,py*2+1);

					if (s2 > s4)
					{
						if (s2 > s6)
						{
							if (s2 > s8)
								m_img_n(px,py) = 2;
							else
								m_img_n(px,py) = 8;
						}
						else
						{
							if (s6 > s8)
								m_img_n(px,py) = 6;
							else
								m_img_n(px,py) = 8;
						}
					}
					else
					{
						if (s4 > s6)
						{
							if (s4 > s8)
								m_img_n(px,py) = 4;
							else
								m_img_n(px,py) = 8;
						}
						else
						{
							if (s6 > s8)
								m_img_n(px,py) = 6;
							else
								m_img_n(px,py) = 8;
						}
					}
				}
			}
	}
}

void Gray2Vec_Grid::NeighborsAdjust()
{
	std::fprintf(stderr,"Optimizing sides...\n");

	// change corners to sides depending on neighbors

	cimg_forXY(m_img_s,px,py)
	{
		if (px > 0)
			if (py > 0)
				if (px < m_img_s.width()-1)
					if (py < m_img_s.height()-1)
					{
						switch (m_img_n(px,py))
						{
							case 1:
								if (m_img_s(px,py) > 255/6)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 6) && Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 2;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 4) && Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 8;
								}
								break;
							case 3:
								if (m_img_s(px,py) > 255/6)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 6) && Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 2;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 8) && Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 4;
								}
								break;
							case 5:
								if (m_img_s(px,py) > 255/6)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 6) && Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 6;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 4) && Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 4;
								}
								break;
							case 7:
								if (m_img_s(px,py) > 255/6)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 4) && Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 8;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 2) && Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 6;
								}
								break;
						}
					}
	}

	std::fprintf(stderr,"Smoothing edges...\n");

	cimg_forXY(m_img_s,px,py)
	{
		if (px > 0)
			if (py > 0)
				if (px < m_img_s.width()-1)
					if (py < m_img_s.height()-1)
					{
						switch (m_img_n(px,py))
						{
							case 1:
								{
									if ((m_img_s(px-1,py) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 8;
									else if ((m_img_s(px,py-1) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 2;
								}
								break;
							case 3:
								{
									if ((m_img_s(px,py-1) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 2;
									else if ((m_img_s(px+1,py) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 4;
								}
								break;
							case 5:
								{
									if ((m_img_s(px+1,py) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 4;
									else if ((m_img_s(px,py+1) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 6;
								}
								break;
							case 7:
								{
									if ((m_img_s(px,py+1) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 6;
									else if ((m_img_s(px-1,py) == 255) && Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 8;
								}
								break;

							case 11:
								{
									if (!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 8;
									else if (!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 2;
								}
								break;
							case 13:
								{
									if (!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 2;
									else if (!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 4;
								}
								break;
							case 15:
								{
									if (!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 4;
									else if (!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 6;
								}
								break;
							case 17:
								{
									if (!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) &&
											!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 6;
									else if (!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) &&
													 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 8;
								}
								break;
						}
					}
	}
}


void Gray2Vec_Grid::ResolveConflicts()
{
	std::fprintf(stderr,"Resolving conflicts (1)...\n");

	size_t cnt1 = 0;

	cimg_forXY(m_img_s,px,py)
	{
		if (px > 0)
			if (py > 0)
				if (px < m_img_s.width()-1)
					if (py < m_img_s.height()-1)
					{
						int n = m_img_n(px,py);

						switch (n)
						{
							case 1:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
									m_img_n(px,py) = 7;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
									m_img_n(px,py) = 3;
								break;
							case 3:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
									m_img_n(px,py) = 1;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
									m_img_n(px,py) = 5;
								break;
							case 5:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
										Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
									m_img_n(px,py) = 3;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
									m_img_n(px,py) = 7;
								break;
							case 7:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
										Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
									m_img_n(px,py) = 5;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
									m_img_n(px,py) = 1;
								break;
								
							case 15:
								if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
									m_img_n(px,py) = 13;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
									m_img_n(px,py) = 17;
								break;
							case 17:
								if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
									m_img_n(px,py) = 15;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
									m_img_n(px,py) = 11;
								break;
							case 11:
								if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
									m_img_n(px,py) = 17;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
									m_img_n(px,py) = 13;
								break;
							case 13:
								if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
									m_img_n(px,py) = 11;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
									m_img_n(px,py) = 15;
								break;
							case 2:
								if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7)) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1)))
									m_img_n(px,py) = 4;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5)) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 3)))
									m_img_n(px,py) = 8;

								break;
							case 6:
								if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7)) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1)))
									m_img_n(px,py) = 4;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5)) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 3)))
									m_img_n(px,py) = 8;
								break;
							case 4:
								if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3)) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1)))
									m_img_n(px,py) = 6;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5)) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7)))
									m_img_n(px,py) = 2;
								break;
							case 8:
								if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3)) && 
										(Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1)))
									m_img_n(px,py) = 6;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5)) && 
												 (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7)))
									m_img_n(px,py) = 2;
								break;
						}
						if (n != m_img_n(px,py))
						{
							if (m_img_f1.width() > 0)
								set_fraction(px,py);
							cnt1++;
						}
					}
	}

	std::fprintf(stderr,"Resolving conflicts (2)...\n");

	size_t cnt2 = 0;

	cimg_forXY(m_img_s,px,py)
	{
		if (px > 0)
			if (py > 0)
				if (px < m_img_s.width()-1)
					if (py < m_img_s.height()-1)
					{
						int n = m_img_n(px,py);
						switch (n)
						{
							case 1:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
									m_img_n(px,py) = 7;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
									m_img_n(px,py) = 3;
								break;
							case 3:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
									m_img_n(px,py) = 1;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
									m_img_n(px,py) = 5;
								break;
							case 5:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
										Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
									m_img_n(px,py) = 3;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
									m_img_n(px,py) = 7;
								break;
							case 7:
								if (!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
									m_img_n(px,py) = 5;
								else if (!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
									m_img_n(px,py) = 1;
								break;
								
							case 15:
								if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
									m_img_n(px,py) = 13;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
									m_img_n(px,py) = 17;
								break;
							case 17:
								if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
									m_img_n(px,py) = 15;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
									m_img_n(px,py) = 11;
								break;
							case 11:
								if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
										Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
									m_img_n(px,py) = 17;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
									m_img_n(px,py) = 13;
								break;
							case 13:
								if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
										!Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
										Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
									m_img_n(px,py) = 11;
								else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
												 Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
									m_img_n(px,py) = 15;
								break;

							case 2:
								if ((Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3)) ||
										(Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1)))
									if (!((Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5)) ||
												(Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))))
										m_img_n(px,py) = 6;
								break;
							case 6:
								if ((Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5)) ||
										(Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7)))
									if (!((Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3)) ||
												(Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))))
										m_img_n(px,py) = 2;
								break;
							case 4:
								if ((Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5)) ||
										(Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3)))
									if (!((Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7)) ||
												(Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))))
										m_img_n(px,py) = 8;
								break;
							case 8:
								if ((Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7)) ||
										(Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) && 
										 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))) 
									if (!((Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5)) ||
												(Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) && 
												 !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))))
										m_img_n(px,py) = 4;
								break;

						}
						if (n != m_img_n(px,py))
						{
							if (m_img_f1.width() > 0)
								set_fraction(px,py);
							cnt2++;
						}
					}
	}

	std::fprintf(stderr,"  changed %ld + %ld pixels\n", cnt1, cnt2);
}


void Gray2Vec_Grid::NeighborsAdjust2()
{
	std::fprintf(stderr,"Adjusting pixel types...\n");

	// change corners to sides depending on neighbors

	size_t cnt = 0;

	cimg_forXY(m_img_s,px,py)
	{
		if (px > 0)
			if (py > 0)
				if (px < m_img_s.width()-1)
					if (py < m_img_s.height()-1)
					{
						int n = m_img_n(px,py);
						switch (n)
						{
							case 1:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 2;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 8;									
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 8;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 2;
								}
								break;
							case 3:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 4;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 2;
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 2;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 4;
								}
								break;
							case 5:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 6;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 4;
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 4;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 6;
								}
								break;
							case 7:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 8;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 6;
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 6;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 8;
								}
								break;

							case 2:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 13;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 11;									
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 11;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 13;
								}
								else if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 1;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 3;									
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 3;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 1;
								}
								break;
							case 4:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 15;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 13;									
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 13;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 15;
								}
								else if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 3;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 5;									
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 5;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 3;
								}
								break;
							case 6:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 17;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 15;									
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 15;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 17;
								}
								else if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 5;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 7;									
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 7;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 5;
								}
								break;
							case 8:
								if (m_img_f1(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 11;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 17;									
								}
								else if (m_img_f2(px,py) == 255)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 17;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 11;
								}
								else if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 7;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 1;									
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 1;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 7;
								}
								break;

							case 11:
								if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 8;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 2;
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3))
										m_img_n(px,py) = 2;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7))
										m_img_n(px,py) = 8;
								}
								break;
							case 13:
								if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 2;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 4;
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5))
										m_img_n(px,py) = 4;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 3) || !Gray2Vec_Grid::check_cover(m_img_n(px,py+1), 1))
										m_img_n(px,py) = 2;
								}
								break;
							case 15:
								if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 4;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 6;
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7))
										m_img_n(px,py) = 6;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 5) || !Gray2Vec_Grid::check_cover(m_img_n(px-1,py), 3))
										m_img_n(px,py) = 4;
								}
								break;
							case 17:
								if (m_img_f1(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 6;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 8;
								}
								else if (m_img_f2(px,py) == 0)
								{
									if (Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px+1,py), 1))
										m_img_n(px,py) = 8;
									else if (Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 7) || !Gray2Vec_Grid::check_cover(m_img_n(px,py-1), 5))
										m_img_n(px,py) = 6;
								}
								break;
						}
						if (n != m_img_n(px,py))
						{
							set_fraction(px,py);
							cnt++;
						}
					}
	}

	std::fprintf(stderr,"  changed %ld pixels\n", cnt);
}

int Gray2Vec_Grid::set_fraction(const int px, const int py)
{
	m_img_f1(px,py) = 0;
	m_img_f2(px,py) = 0;
	m_img_f3(px,py) = -1;

	switch (m_img_n(px,py))
	{
		case 1:
		case 3:
		case 5:
		case 7:
			m_img_f1(px,py) = std::sqrt(m_img_s(px,py)*255.0*2);
			m_img_f2(px,py) = m_img_f1(px,py);
			return 1;
			break;
		case 2:
		case 4:
		case 6:
		case 8:
			m_img_f1(px,py) = m_img_s(px,py);
			m_img_f2(px,py) = m_img_s(px,py);
			return 2;
			break;
		case 11:
		case 13:
		case 15:
		case 17:
			m_img_f1(px,py) = 255-std::sqrt(((255-m_img_s(px,py))*255.0)*2);
			m_img_f2(px,py) = m_img_f1(px,py);
			return 3;
			break;
	}

	return 0;
}

double Gray2Vec_Grid::pixel_error(const int px, const int py, const bool use_adjust)
{
	switch (m_img_n(px,py))
	{
		case 1:
		case 3:
		case 5:
		case 7:
			if (use_adjust)
				return 0.5*m_img_f1(px,py)*m_img_f2(px,py)*m_img_f3(px,py)/(255*255) - m_img_s(px,py);
			else
				return 0.5*m_img_f1(px,py)*m_img_f2(px,py)/255 - m_img_s(px,py);
			break;
		case 2:
		case 4:
		case 6:
		case 8:
			if (use_adjust)
				return (2.0*m_img_f3(px,py)+m_img_f1(px,py)+m_img_f2(px,py))*0.25 - m_img_s(px,py);
			else
				return (m_img_f1(px,py)+m_img_f2(px,py))*0.5 - m_img_s(px,py);
			break;
		case 11:
		case 13:
		case 15:
		case 17:
			if (use_adjust)
				return 0.5*(255-m_img_f1(px,py))*(255-m_img_f2(px,py))*m_img_f3(px,py)/(255*255) - (255-m_img_s(px,py));
			else
				return 0.5*(255-m_img_f1(px,py))*(255-m_img_f2(px,py))/255 - (255-m_img_s(px,py));
			break;
	}
}


void Gray2Vec_Grid::InitFractions()
{
	std::fprintf(stderr,"Determining initial fractions...\n");

	// fractions apply clockwise to the pixel sides
	m_img_f1 = CImg<unsigned char>(m_img.width()/2, m_img.height()/2, 1, 1);
	m_img_f2 = CImg<unsigned char>(m_img.width()/2, m_img.height()/2, 1, 1);
	m_img_f3 = CImg<short>(m_img.width()/2, m_img.height()/2, 1, 1);

	size_t cnt_all = 0;
	size_t cnt_corner = 0;
	size_t cnt_side = 0;
	size_t cnt_corner2 = 0;

	cimg_forXY(m_img_s,px,py)
	{
		cnt_all++;

		switch (set_fraction(px,py))
		{
			case 1:
				cnt_corner++;
				break;
			case 2:
				cnt_side++;
				break;
			case 3:
				cnt_corner2++;
				break;
		}
	}

	std::fprintf(stderr,"  %ld side, %ld + %ld corner of %ld pixels\n", cnt_side, cnt_corner, cnt_corner2, cnt_all);
}

void Gray2Vec_Grid::TuneFractions(const double max_error)
{
	std::fprintf(stderr,"Tuning fractions...\n");

	double f;
	double df;

	double df_max = 0.0;
	double df_sum = 0.0;
	size_t df_cnt = 0;

	size_t cnt_all = 0;
	size_t cnt_corner = 0;
	size_t cnt_corner_n = 0;
	size_t cnt_corner_s = 0;
	size_t cnt_side = 0;
	size_t cnt_side_n = 0;
	size_t cnt_side_s = 0;
	size_t cnt_corner2 = 0;
	size_t cnt_corner2_n = 0;
	size_t cnt_corner2_s = 0;
	size_t cnt_changed = 0;
	size_t cnt_changed2 = 0;
	size_t cnt_changed_type = 0;

	CImg<unsigned char> img_e = CImg<unsigned char>(m_img.width()/2, m_img.height()/2, 1, 1);

	cimg_forXY(m_img_s,px,py)
	{
		cnt_all++;

		df = pixel_error(px, py, false);

		switch (m_img_n(px,py))
		{
			case 1:
			case 3:
			case 5:
			case 7:
				if ((max_error > 0) && (std::abs(df) > max_error*255))
				{
					if (m_img_f1(px,py)*m_img_f2(px,py) > 0)
						f = double(m_img_s(px,py)*255)/(0.5*m_img_f1(px,py)*m_img_f2(px,py));
					else
						f = 1.0;

					double m = 255.0/std::max(m_img_f1(px,py), m_img_f2(px,py));
					m_img_f3(px,py) = std::min(std::max(f, 0.0), m)*255;

					df = pixel_error(px, py, true);

					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					img_e(px,py) = std::abs(df);

					cnt_changed++;
				}
				else
				{
					int s = sides_connected(px,py);

					if (s == 0)
					{
						m_img_f1(px,py) = std::sqrt(m_img_s(px,py)*255.0*2);
						m_img_f2(px,py) = m_img_f1(px,py);
						cnt_corner_n++;
					}
					else if (s == 1)
					{
						if (m_img_f1(px,py) > 0)
							f = m_img_s(px,py)*255*2/m_img_f1(px,py);
						else
							f = m_img_s(px,py)*255*2;
						if (f > 255.0)
						{
							m_img_f2(px,py) = 255;
						}
						else
						{
							m_img_f2(px,py) = f;
						}
						cnt_corner_s++;
					}
					else if (s == 2)
					{
						if (m_img_f2(px,py) > 0)
							f = m_img_s(px,py)*255*2/m_img_f2(px,py);
						else
							f = m_img_s(px,py)*255*2;
						if (f > 255.0)
						{
							m_img_f1(px,py) = 255;
						}
						else
						{
							m_img_f1(px,py) = f;
						}
						cnt_corner_s++;
					}
					else
					{
						m_img_f1(px,py) = 0.75*m_img_f1(px,py) + 0.25*std::sqrt(m_img_s(px,py)*255*2);
						m_img_f2(px,py) = 0.75*m_img_f2(px,py) + 0.25*std::sqrt(m_img_s(px,py)*255*2);
						cnt_corner++;
					}

					df = pixel_error(px, py, false);

					m_img_f3(px,py) = -1;

					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					img_e(px,py) = std::abs(df);

				}
					
				break;
			case 2:
			case 4:
			case 6:
			case 8:
				if ((max_error > 0) && (std::abs(df) > max_error*255))
				{
					f = 2.0*(m_img_s(px,py) - 0.25*m_img_f1(px,py)-0.25*m_img_f2(px,py));
					m_img_f3(px,py) = std::min(std::max(f, 0.0), 255.0);

					df = pixel_error(px, py, true);

					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					img_e(px,py) = std::abs(df);

					cnt_changed++;
				}
				else
				{
					int s = sides_connected(px,py);

					if (s == 0)
					{
						m_img_f1(px,py) = m_img_s(px,py);
						m_img_f2(px,py) = m_img_s(px,py);
						cnt_side_n++;
					}
					else if (s == 1)
					{
						f = 2.0*(m_img_s(px,py)-0.5*m_img_f1(px,py));
						if (f < 0.0)
						{
							m_img_f2(px,py) = 0;
						}
						else if (f > 255.0)
						{
							m_img_f2(px,py) = 255;
						}
						else
						{
							m_img_f2(px,py) = f;
						}
						cnt_side_s++;
					}
					else if (s == 2)
					{
						f = 2.0*(m_img_s(px,py)-0.5*m_img_f2(px,py));
						if (f < 0.0)
						{
							m_img_f1(px,py) = 0.0;
						}
						else if (f > 255.0)
						{
							m_img_f1(px,py) = 255;
						}
						else
						{
							m_img_f1(px,py) = f;
						}
						cnt_side_s++;
					}
					else
					{
						m_img_f1(px,py) = 0.75*m_img_f1(px,py) + 0.25*m_img_s(px,py);
						m_img_f2(px,py) = 0.75*m_img_f2(px,py) + 0.25*m_img_s(px,py);
						cnt_side++;
					}

					df = pixel_error(px, py, false);

					m_img_f3(px,py) = -1;
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					img_e(px,py) = std::abs(df);
				}

				break;
			case 11:
			case 13:
			case 15:
			case 17:
				if ((max_error > 0) && (std::abs(df) > max_error*255))
				{
					if ((255-m_img_f1(px,py))*(255-m_img_f2(px,py)) > 0)
						f = double((255-m_img_s(px,py))*255)/(0.5*(255-m_img_f1(px,py))*(255-m_img_f2(px,py)));
					else
						f = 1.0;

					double m = 255.0/std::max((255-m_img_f1(px,py)), (255-m_img_f2(px,py)));
					m_img_f3(px,py) = std::min(std::max(f, 0.0), m)*255;

					df = pixel_error(px, py, true);

					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					img_e(px,py) = std::abs(df);

					cnt_changed2++;
				}
				else
				{
					int s = sides_connected(px,py);

					if (s == 0)
					{
						m_img_f1(px,py) = 255-std::sqrt(((255-m_img_s(px,py))*255.0)*2);
						m_img_f2(px,py) = m_img_f1(px,py);
						cnt_corner2_n++;
					}
					else if (s == 1)
					{
						if (m_img_f1(px,py) < 255)
							f = 255.0 - (255-m_img_s(px,py))*255*2/(255-m_img_f1(px,py));
						else
							f = 255.0 - (255-m_img_s(px,py))*255*2;
						if (f < 0.0)
						{
							m_img_f2(px,py) = 0.0;
						}
						else
						{
							m_img_f2(px,py) = f;
						}
						cnt_corner2_s++;
					}
					else if (s == 2)
					{
						if (m_img_f2(px,py) < 255)
							f = 255.0 - (255-m_img_s(px,py))*255*2/(255-m_img_f2(px,py));
						else
							f = 255.0 - (255-m_img_s(px,py))*255*2;
						if (f < 0.0)
						{
							m_img_f1(px,py) = 0.0;
						}
						else
						{
							m_img_f1(px,py) = f;
						}
						cnt_corner2_s++;
					}
					else
					{
						m_img_f1(px,py) = 0.75*m_img_f1(px,py) + 0.25*(255-std::sqrt(((255-m_img_s(px,py))*255.0*2)));
						m_img_f2(px,py) = 0.75*m_img_f2(px,py) + 0.25*(255-std::sqrt(((255-m_img_s(px,py))*255.0*2)));
						cnt_corner2++;
					}

					df = pixel_error(px, py, false);

					m_img_f3(px,py) = -1;
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					img_e(px,py) = std::abs(df);
				}
				break;
			default:
				img_e(px,py) = 0;
				break;
		}
	}

	if (df_cnt > 0)
	{
		std::fprintf(stderr,"  %ld + %ld + %ld side, %ld + %ld + %ld small corner , %ld + %ld + %ld large corner\n", cnt_side, cnt_side_n, cnt_side_s, cnt_corner, cnt_corner_n, cnt_corner_s, cnt_corner2, cnt_corner2_n, cnt_corner2_s);

		std::fprintf(stderr,"  maximum error: %.1f, average: %.1f (%.3f), %ld + %ld + %ld changed\n", df_max, df_sum/df_cnt, (df_sum/df_cnt)/255, cnt_changed_type, cnt_changed, cnt_changed2);
	}

	if (m_debug)  img_e.save("debug-e.tif");
}

void Gray2Vec_Grid::move_dir(int &px, int &py, const int dir)
{
	switch (dir)
	{
		case 1:
			px -= 1;
			py -= 1;
			break;
		case 2:
			py -= 1;
			break;
		case 3:
			px += 1;
			py -= 1;
			break;
		case 4:
			px += 1;
			break;
		case 5:
			px += 1;
			py += 1;
			break;
		case 6:
			py += 1;
			break;
		case 7:
			px -= 1;
			py += 1;
			break;
		case 8:
			px -= 1;
			break;
	}
}


int Gray2Vec_Grid::side1(const int dir)
{
	switch (dir)
	{
		case 1:
			return 2;
			break;
		case 2:
			return 4;
			break;
		case 3:
			return 4;
			break;
		case 4:
			return 6;
			break;
		case 5:
			return 6;
			break;
		case 6:
			return 8;
			break;
		case 7:
			return 8;
			break;
		case 8:
			return 2;
			break;
		case 11:
			return 4;
			break;
		case 13:
			return 6;
			break;
		case 15:
			return 8;
			break;
		case 17:
			return 2;
			break;
	}
	return 0;
}

int Gray2Vec_Grid::side2(const int dir)
{
	switch (dir)
	{
		case 1:
			return 8;
			break;
		case 2:
			return 8;
			break;
		case 3:
			return 2;
			break;
		case 4:
			return 2;
			break;
		case 5:
			return 4;
			break;
		case 6:
			return 4;
			break;
		case 7:
			return 6;
			break;
		case 8:
			return 6;
			break;
		case 11:
			return 6;
			break;
		case 13:
			return 8;
			break;
		case 15:
			return 2;
			break;
		case 17:
			return 4;
			break;
	}
	return 0;
}


int Gray2Vec_Grid::share_sides(const int px1, const int py1, const int px2, const int py2)
{
	if (px1 < 0) return 0;
	if (py1 < 0) return 0;
	if (px2 < 0) return 0;
	if (py2 < 0) return 0;

	if (px1 >= m_img_n.width()) return 0;
	if (py1 >= m_img_n.height()) return 0;
	if (px2 >= m_img_n.width()) return 0;
	if (py2 >= m_img_n.height()) return 0;

	int n1 = m_img_n(px1, py1);
	int n2 = m_img_n(px2, py2);

	int s11 = Gray2Vec_Grid::side1(n1);
	int s12 = Gray2Vec_Grid::side2(n1);
	int s21 = Gray2Vec_Grid::side1(n2);
	int s22 = Gray2Vec_Grid::side2(n2);

	int pxn11 = px1;
	int pyn11 = py1;

	int pxn21 = px2;
	int pyn21 = py2;

	int pxn12 = px1;
	int pyn12 = py1;

	int pxn22 = px2;
	int pyn22 = py2;

	Gray2Vec_Grid::move_dir(pxn11,pyn11,s11);
	Gray2Vec_Grid::move_dir(pxn21,pyn21,s21);
	Gray2Vec_Grid::move_dir(pxn12,pyn12,s12);
	Gray2Vec_Grid::move_dir(pxn22,pyn22,s22);

	if ((pxn11 == px2) && (pyn11 == py2) && (pxn21 == px1) && (pyn21 == py1)) return 11;
	if ((pxn11 == px2) && (pyn11 == py2) && (pxn22 == px1) && (pyn22 == py1)) return 12;

	if ((pxn12 == px2) && (pyn12 == py2) && (pxn21 == px1) && (pyn21 == py1)) return 21;
	if ((pxn12 == px2) && (pyn12 == py2) && (pxn22 == px1) && (pyn22 == py1)) return 22;

	return 0;
}

int Gray2Vec_Grid::sides_connected(const int px, const int py)
{
	if (px < 0) return 0;
	if (py < 0) return 0;

	if (px >= m_img_n.width()) return 0;
	if (py >= m_img_n.height()) return 0;

	int n = m_img_n(px, py);

	int s1 = Gray2Vec_Grid::side1(n);
	int s2 = Gray2Vec_Grid::side2(n);

	int pxn1 = px;
	int pyn1 = py;

	int pxn2 = px;
	int pyn2 = py;

	Gray2Vec_Grid::move_dir(pxn1,pyn1,s1);
	Gray2Vec_Grid::move_dir(pxn2,pyn2,s2);

	int s11 = 0;
	int s21 = 0;
	int s12 = 0;
	int s22 = 0;

	if (pxn1 < 0) s11 = -100;
	if (pyn1 < 0) s11 = -100;

	if (pxn1 >= m_img_n.width()) s11 = -100;
	if (pyn1 >= m_img_n.height()) s11 = -100;

	if (pxn2 < 0) s21 = -100;
	if (pyn2 < 0) s21 = -100;

	if (pxn2 >= m_img_n.width()) s21 = -100;
	if (pyn2 >= m_img_n.height()) s21 = -100;

	int res = 0;

	if (s11 == 0)
	{
		s11 = Gray2Vec_Grid::side1(m_img_n(pxn1, pyn1));
		s12 = Gray2Vec_Grid::side2(m_img_n(pxn1, pyn1));
		if (std::abs(s1-s11) == 4) res += 1;
		else if (std::abs(s1-s12) == 4) res += 1;
	}

	if (s21 == 0)
	{
		s21 = Gray2Vec_Grid::side1(m_img_n(pxn2, pyn2));
		s22 = Gray2Vec_Grid::side2(m_img_n(pxn2, pyn2));
		if (std::abs(s2-s21) == 4) res += 2;
		else if (std::abs(s2-s22) == 4) res += 2;
	}

	return res;
}


void Gray2Vec_Grid::FractionsNeighborsAdj()
{
	std::fprintf(stderr,"Adjusting neighbor fractions to match...\n");

	// this adjust fractions of neighboring pixels to match
	// this resolves cases where neighboring pixels have
	// exactly opposite fractions at their respective sides

	double df;
	double df_max = 0.0;
	double df_sum = 0.0;
	size_t df_cnt = 0;

	cimg_forXY(m_img_s,px,py)
	{
		int nx, ny, nn, avg;

		for (int i = 0; i < 4; i++)
		{
			nx = px + x4[i];
			ny = py + y4[i];

			nn = Gray2Vec_Grid::share_sides(px, py, nx, ny);

			switch (nn)
			{
				case 11:
					avg = (m_img_f1(px,py) + m_img_f1(nx,ny))/2;

					df = avg - m_img_f1(px,py);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f1(px,py) = (m_img_f1(px,py)+avg)/2;

					df = avg - m_img_f1(nx,ny);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f1(nx,ny) = (m_img_f1(nx,ny)+avg)/2;
					break;

				case 12:
					avg = (m_img_f1(px,py) + m_img_f2(nx,ny))/2;

					df = avg - m_img_f1(px,py);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f1(px,py) = (m_img_f1(px,py)+avg)/2;

					df = avg - m_img_f2(nx,ny);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f2(nx,ny) = (m_img_f2(nx,ny)+avg)/2;
					break;
				case 21:
					avg = (m_img_f2(px,py) + m_img_f1(nx,ny))/2;

					df = avg - m_img_f2(px,py);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f2(px,py) = (m_img_f2(px,py)+avg)/2;

					df = avg - m_img_f1(nx,ny);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;
				
					m_img_f1(nx,ny) = (m_img_f1(nx,ny)+avg)/2;
					break;
				case 22:
					avg = (m_img_f2(px,py) + m_img_f2(nx,ny))/2;

					df = avg - m_img_f2(px,py);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f2(px,py) = (m_img_f2(px,py)+avg)/2;

					df = avg - m_img_f2(nx,ny);
					df_max = std::max(std::abs(df), df_max);
					df_sum += std::abs(df);
					df_cnt++;

					m_img_f2(nx,ny) = (m_img_f2(nx,ny)+avg)/2;
					break;
			}
		}
	}

	if (df_cnt > 0)
		std::fprintf(stderr,"  %ld pairs, maximum error: %.3f, average: %.3f\n", df_cnt/2, df_max, df_sum/df_cnt);
}

bool Gray2Vec_Grid::Vectorize(const std::string file, const std::string layer, const bool Append)
{
	std::fprintf(stderr,"Generating subgrid...\n");

	CImg<unsigned char> img_h = CImg<unsigned char>(m_img.width(), m_img.height(), 1, 1);

	if (m_debug)
	{
		m_img_s.save("debug-s.tif");
		m_img_n.save("debug-n.tif");
		m_img_f1.save("debug-f1.tif");
		m_img_f2.save("debug-f2.tif");
		m_img_f3.save("debug-f3.tif");
	}

	cimg_forXY(m_img_s,px,py)
	{
		if (m_img_s(px,py) == 0)
		{
			img_h(px*2,py*2) = 0;
			img_h(px*2+1,py*2) = 0;
			img_h(px*2,py*2+1) = 0;
			img_h(px*2+1,py*2+1) = 0;
		}
		else if (m_img_s(px,py) == 255)
		{
			img_h(px*2,py*2) = 255;
			img_h(px*2+1,py*2) = 255;
			img_h(px*2,py*2+1) = 255;
			img_h(px*2+1,py*2+1) = 255;
		}
		else
		{
			switch (m_img_n(px,py))
			{
				case 1:
					img_h(px*2,py*2) = 255;
					img_h(px*2+1,py*2) = 0;
					img_h(px*2,py*2+1) = 0;
					img_h(px*2+1,py*2+1) = 0;
					break;
				case 2:
					img_h(px*2,py*2) = 255;
					img_h(px*2+1,py*2) = 255;
					img_h(px*2,py*2+1) = 0;
					img_h(px*2+1,py*2+1) = 0;
					break;
				case 3:
					img_h(px*2,py*2) = 0;
					img_h(px*2+1,py*2) = 255;
					img_h(px*2,py*2+1) = 0;
					img_h(px*2+1,py*2+1) = 0;
					break;
				case 4:
					img_h(px*2,py*2) = 0;
					img_h(px*2+1,py*2) = 255;
					img_h(px*2,py*2+1) = 0;
					img_h(px*2+1,py*2+1) = 255;
					break;
				case 5:
					img_h(px*2,py*2) = 0;
					img_h(px*2+1,py*2) = 0;
					img_h(px*2,py*2+1) = 0;
					img_h(px*2+1,py*2+1) = 255;
					break;
				case 6:
					img_h(px*2,py*2) = 0;
					img_h(px*2+1,py*2) = 0;
					img_h(px*2,py*2+1) = 255;
					img_h(px*2+1,py*2+1) = 255;
					break;
				case 7:
					img_h(px*2,py*2) = 0;
					img_h(px*2+1,py*2) = 0;
					img_h(px*2,py*2+1) = 255;
					img_h(px*2+1,py*2+1) = 0;
					break;
				case 8:
					img_h(px*2,py*2) = 255;
					img_h(px*2+1,py*2) = 0;
					img_h(px*2,py*2+1) = 255;
					img_h(px*2+1,py*2+1) = 0;
					break;
				case 11:
					img_h(px*2,py*2) = 255;
					img_h(px*2+1,py*2) = 255;
					img_h(px*2,py*2+1) = 255;
					img_h(px*2+1,py*2+1) = 0;
					break;
				case 13:
					img_h(px*2,py*2) = 255;
					img_h(px*2+1,py*2) = 255;
					img_h(px*2,py*2+1) = 0;
					img_h(px*2+1,py*2+1) = 255;
					break;
				case 15:
					img_h(px*2,py*2) = 0;
					img_h(px*2+1,py*2) = 255;
					img_h(px*2,py*2+1) = 255;
					img_h(px*2+1,py*2+1) = 255;
					break;
				case 17:
					img_h(px*2,py*2) = 255;
					img_h(px*2+1,py*2) = 0;
					img_h(px*2,py*2+1) = 255;
					img_h(px*2+1,py*2+1) = 255;
					break;
			}
		}

	}

	std::fprintf(stderr,"Preparing vector file...\n");

	const char *pszDriverName = "SQLite";
	const char *Options[] = { "SPATIALITE=TRUE", "INIT_WITH_EPSG=no", NULL };
	GDALDriverH hDriver;
	GDALDatasetH hDS;
	OGRLayerH hLayer;
	OGRFieldDefnH hFieldDefn;
	double x, y;
	bool CreateLayer = false;

	GDALAllRegister();

	CPLSetConfigOption("OGR_SQLITE_SYNCHRONOUS", "OFF");

	if (Append)
	{
		std::fprintf(stderr,"  Opening %s to append...\n", file.c_str());

#if GDAL_VERSION_MAJOR >= 2
		hDS = GDALOpenEx(file.c_str(), GDAL_OF_VECTOR || GDAL_OF_UPDATE, NULL, NULL, NULL );
#else
		hDS = OGROpen(file.c_str(), TRUE, NULL);
#endif

		if (hDS == NULL)
		{
			fprintf(stderr, "Opening output file failed.\n");
			return false;
		}

#if GDAL_VERSION_MAJOR >= 2
		hLayer = GDALDatasetGetLayerByName(hDS, layer.c_str() );
#else
		hLayer = OGR_DS_GetLayerByName(hDS, layer.c_str() );
#endif

		if( hLayer == NULL ) CreateLayer = true;
	}
	else
	{

#if GDAL_VERSION_MAJOR >= 2
		hDriver = GDALGetDriverByName( pszDriverName );
#else
		hDriver = OGRGetDriverByName( pszDriverName );
#endif

		if (hDriver == NULL)
		{
			fprintf(stderr, "%s driver not available.\n", pszDriverName);
			return false;
		}

#if GDAL_VERSION_MAJOR >= 2
		hDS = GDALCreate( hDriver, file.c_str(), 0, 0, 0, GDT_Unknown, const_cast<char**>(Options));
#else
		hDS = OGR_Dr_CreateDataSource( hDriver, file.c_str(), const_cast<char**>(Options));
#endif

		if (hDS == NULL)
		{
			fprintf(stderr, "Creation of output file failed.\n");
			return false;
		}

		CreateLayer = true;
	}

	if (CreateLayer)
	{
#if GDAL_VERSION_MAJOR >= 2
		hLayer = GDALDatasetCreateLayer(hDS, layer.c_str(), m_SRS, wkbPolygon, NULL);
#else
		hLayer = OGR_DS_CreateLayer(hDS, layer.c_str(), m_SRS, wkbPolygon, NULL);
#endif

		if (hLayer == NULL)
		{
			fprintf(stderr, "Layer creation failed.\n");
			return false;
		}

		if (m_x >= 0)
		{
			hFieldDefn = OGR_Fld_Create( "x", OFTInteger );
			OGR_Fld_SetWidth( hFieldDefn, 12);
			if( OGR_L_CreateField( hLayer, hFieldDefn, TRUE ) != OGRERR_NONE )
			{
				fprintf(stderr, "Creating attribute field x failed.\n" );
				return false;
			}

			OGR_Fld_Destroy(hFieldDefn);
		}

		if (m_y >= 0)
		{
			hFieldDefn = OGR_Fld_Create( "y", OFTInteger );
			OGR_Fld_SetWidth( hFieldDefn, 12);
			if( OGR_L_CreateField( hLayer, hFieldDefn, TRUE ) != OGRERR_NONE )
			{
				fprintf(stderr, "Creating attribute field y failed.\n" );
				return false;
			}

			OGR_Fld_Destroy(hFieldDefn);
		}

		if (m_z >= 0)
		{
			hFieldDefn = OGR_Fld_Create( "z", OFTInteger );
			OGR_Fld_SetWidth( hFieldDefn, 12);
			if( OGR_L_CreateField( hLayer, hFieldDefn, TRUE ) != OGRERR_NONE )
			{
				fprintf(stderr, "Creating attribute field z failed.\n" );
				return false;
			}

			OGR_Fld_Destroy(hFieldDefn);
		}
	}

	std::fprintf(stderr,"Vectorizing grid...\n");

	Polygonize(img_h, hLayer);

#if GDAL_VERSION_MAJOR >= 2
	GDALClose(hDS);
#else
	OGR_DS_Destroy(hDS);
#endif

}

/*
 * This method is derived from polygonize.cpp from the gdal source package
 * which comes with the following copyright notice:
 *
 ******************************************************************************
 * Copyright (c) 2008, Frank Warmerdam
 * Copyright (c) 2009-2011, Even Rouault <even dot rouault at mines-paris dot org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/
bool Gray2Vec_Grid::EmitPolygonToLayer(OGRLayerH hOutLayer, RPolygon *poRPoly)
{
	OGRFeatureH hFeat;
	OGRGeometryH hPolygon;

	/* -------------------------------------------------------------------- */
	/*      Turn bits of lines into coherent rings.                         */
	/* -------------------------------------------------------------------- */
	poRPoly->Coalesce();

	/* -------------------------------------------------------------------- */
	/*      Create the polygon geometry.                                    */
	/* -------------------------------------------------------------------- */
	size_t iString;

	hPolygon = OGR_G_CreateGeometry( wkbPolygon );

	for( iString = 0; iString < poRPoly->aanXY.size(); iString++ )
	{
		std::vector<int> &anString = poRPoly->aanXY[iString];
		OGRGeometryH hRing = OGR_G_CreateGeometry( wkbLinearRing );

		int iVert;

		// we go last to first to ensure the linestring is allocated to
		// the proper size on the first try.
		for (iVert = 0; iVert < anString.size()/2; iVert++ )
		{
			double dfX, dfY;
			int    nPixelX, nPixelY;

			nPixelX = anString[iVert*2];
			nPixelY = anString[iVert*2+1];

			double fx = nPixelX;
			double fy = nPixelY;

			int px = nPixelX/2;
			int py = nPixelY/2;

			double fA = -1.0;
			double fB = -1.0;
			double f = 0.5;

			if ((nPixelX % 2) == 0)
			{
				if ((nPixelY % 2) != 0)
				{
					// vertical middle

					// right side
					switch (m_img_n(px,py))
					{
						case 1:
						case 2:
						case 13:
							fA = m_img_f2(px,py);
							break;
						case 15:
						case 6:
						case 7:
							fA = 255-m_img_f1(px,py);
							break;
					}

					// left side
					if (px > 0)
					switch (m_img_n(px-1,py))
					{
						case 11:
						case 2:
						case 3:
							fB = m_img_f1(px-1,py);
							break;
						case 5:
						case 6:
						case 17:
							fB = 255-m_img_f2(px-1,py);
							break;
					}

					// average both sides
					if (fA >= 0.0)
					{
						if (fB >= 0.0)
							f = 0.5*(fA+fB)/255;
						else
							f = fA/255;
					}
					else
						if (fB >= 0.0)
							f = fB/255;

					// apply offset to fractional coordinate
					fy += 2.0*(f-0.5);
				}
			}
			else
			{
				if ((nPixelY % 2) == 0)
				{
					// horizontal middle

					// bottom side
					switch (m_img_n(px,py))
					{
						case 17:
						case 8:
						case 1:
							fA = m_img_f1(px,py);
							break;
						case 3:
						case 4:
						case 15:
							fA = 255-m_img_f2(px,py);
							break;
					}

					// top side
					if (py > 0)
					switch (m_img_n(px,py-1))
					{
						case 7:
						case 8:
						case 11:
							fB = m_img_f2(px,py-1);
							break;
						case 13:
						case 4:
						case 5:
							fB = 255-m_img_f1(px,py-1);
							break;
					}

					// average both sides
					if (fA >= 0.0)
					{
						if (fB >= 0.0)
							f = 0.5*(fA+fB)/255;
						else
							f = fA/255;
					}
					else
						if (fB >= 0.0)
							f = fB/255;

					// apply offset to fractional coordinate
					fx += 2.0*(f-0.5);
				}
				else
				{
					double df;

					// middle in both directions: unnecessary point - unless necessary to limit error
					switch (m_img_n(px,py))
					{
						case 1:
						case 3:
						case 5:
						case 7:
							if (m_img_f3(px,py) >= 0)
							{
								if (m_img_n(px,py) == 1)
								{
									fx += 2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)-1;
									fy += 2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)-1;
								}
								else if (m_img_n(px,py) == 3)
								{
									fx += -2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)+1;
									fy += 2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)-1;
								}
								else if (m_img_n(px,py) == 5)
								{
									fx += -2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)+1;
									fy += -2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)+1;
								}
								else if (m_img_n(px,py) == 7)
								{
									fx += 2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)-1;
									fy += -2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)+1;
								}
								//std::fprintf(stderr," Precalculated error compenation corner at (%d/%d): %d, %d\n", px, py, m_img_s(px,py), m_img_f3(px,py));
							}
							else
								continue;

							break;
						case 2:
						case 4:
						case 6:
						case 8:
							if (m_img_f3(px,py) >= 0)
							{
								if (m_img_n(px,py) == 2)
								{
									fy += 2.0*double(m_img_f3(px,py))/255-1;
								}
								else if (m_img_n(px,py) == 4)
								{
									fx += -2.0*double(m_img_f3(px,py))/255+1;
								}
								else if (m_img_n(px,py) == 6)
								{
									fy += -2.0*double(m_img_f3(px,py))/255+1;
								}
								else if (m_img_n(px,py) == 8)
								{
									fx += 2.0*double(m_img_f3(px,py))/255-1;
								}
							}
							else
								continue;

							break;
						case 11:
						case 13:
						case 15:
						case 17:
							if (m_img_f3(px,py) >= 0)
							{
								if (m_img_n(px,py) == 1)
								{
									fx += 2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)-1;
									fy += 2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)-1;
								}
								else if (m_img_n(px,py) == 3)
								{
									fx += -2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)+1;
									fy += 2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)-1;
								}
								else if (m_img_n(px,py) == 5)
								{
									fx += -2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)+1;
									fy += -2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)+1;
								}
								else if (m_img_n(px,py) == 7)
								{
									fx += 2.0*double(m_img_f1(px,py)*m_img_f3(px,py))/(255*255)-1;
									fy += -2.0*double(m_img_f2(px,py)*m_img_f3(px,py))/(255*255)+1;
								}
								//std::fprintf(stderr," Precalculated error compenation corner at (%d/%d): %d, %d\n", px, py, m_img_s(px,py), m_img_f3(px,py));
							}
							else
								continue;
							break;
					}
					//continue;
				}
			}

			dfX = m_GeoTransform[0]
				+ fx * m_GeoTransform[1]
				+ fy * m_GeoTransform[2];
			dfY = m_GeoTransform[3]
				+ fx * m_GeoTransform[4]
				+ fy * m_GeoTransform[5];

			OGR_G_AddPoint_2D(hRing, dfX, dfY);
		}

		OGR_G_AddGeometryDirectly( hPolygon, hRing );
	}

	/* -------------------------------------------------------------------- */
	/*      Create the feature object.                                      */
	/* -------------------------------------------------------------------- */
	hFeat = OGR_F_Create( OGR_L_GetLayerDefn( hOutLayer ) );

	if (m_x >= 0)
		OGR_F_SetFieldInteger( hFeat, OGR_F_GetFieldIndex(hFeat, "x"), m_x );
	if (m_y >= 0)
		OGR_F_SetFieldInteger( hFeat, OGR_F_GetFieldIndex(hFeat, "y"), m_y );
	if (m_z >= 0)
		OGR_F_SetFieldInteger( hFeat, OGR_F_GetFieldIndex(hFeat, "z"), m_z );

	OGR_G_CloseRings( hPolygon );

	OGR_F_SetGeometryDirectly( hFeat, hPolygon );

	/* -------------------------------------------------------------------- */
	/*      Write the to the layer.                                         */
	/* -------------------------------------------------------------------- */
	bool Res = true;

	if( OGR_L_CreateFeature( hOutLayer, hFeat ) != OGRERR_NONE ) Res = false;

	OGR_F_Destroy( hFeat );

	return Res;
}

/*
 * This method is derived from polygonize.cpp from the gdal source package
 * which comes with the following copyright notice:
 *
 ******************************************************************************
 * Copyright (c) 2008, Frank Warmerdam
 * Copyright (c) 2009-2011, Even Rouault <even dot rouault at mines-paris dot org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/
bool Gray2Vec_Grid::Polygonize(CImg<unsigned char> &img_h, OGRLayerH hOutLayer)
{
	int nConnectedness = 4;

	if( !OGR_L_TestCapability( hOutLayer, OLCSequentialWrite ) )
	{
		fprintf(stderr, "Output feature layer does not appear to support creation\nof features in GDALPolygonize().\n");
		return false;
	}

	/* -------------------------------------------------------------------- */
	/*      Allocate working buffers.                                       */
	/* -------------------------------------------------------------------- */
	bool Res = true;
	int nXSize = img_h.width();
	int nYSize = img_h.height();

#if GDAL_VERSION_MAJOR >= 2
	int *panLastLineVal = (int *) VSI_MALLOC2_VERBOSE(sizeof(int),nXSize + 2);
	int *panThisLineVal = (int *) VSI_MALLOC2_VERBOSE(sizeof(int),nXSize + 2);
	GInt32 *panLastLineId =  (GInt32 *) VSI_MALLOC2_VERBOSE(sizeof(GInt32),nXSize + 2);
	GInt32 *panThisLineId =  (GInt32 *) VSI_MALLOC2_VERBOSE(sizeof(GInt32),nXSize + 2);
#else
	GInt32 *panLastLineVal = (GInt32 *) VSIMalloc2(sizeof(GInt32),nXSize + 2);
	int *panThisLineVal = (int *) VSIMalloc2(sizeof(int),nXSize + 2);
	GInt32 *panLastLineId =  (GInt32 *) VSIMalloc2(sizeof(GInt32),nXSize + 2);
	GInt32 *panThisLineId =  (GInt32 *) VSIMalloc2(sizeof(GInt32),nXSize + 2);
#endif
	GByte *pabyMaskLine = NULL;
	if (panLastLineVal == NULL || panThisLineVal == NULL ||
			panLastLineId == NULL || panThisLineId == NULL)
	{
		CPLFree( panThisLineId );
		CPLFree( panLastLineId );
		CPLFree( panThisLineVal );
		CPLFree( panLastLineVal );
		CPLFree( pabyMaskLine );
		return false;
	}

	/* -------------------------------------------------------------------- */
	/*      The first pass over the raster is only used to build up the     */
	/*      polygon id map so we will know in advance what polygons are     */
	/*      what on the second pass.                                        */
	/* -------------------------------------------------------------------- */
	int iY;

	GDALRasterPolygonEnumeratorT<int, IntEqualityTest> oFirstEnum(nConnectedness);

	for( iY = 0; Res && iY < nYSize; iY++ )
	{

		for (int iX = 0; iX < nXSize; iX++)
		{
			if (img_h(iX,iY) > 0)
				panThisLineVal[iX] = img_h(iX,iY);
			else
				panThisLineVal[iX] = GP_NODATA_MARKER;
		}

		if( iY == 0 )
			oFirstEnum.ProcessLine(
														 NULL, panThisLineVal, NULL, panThisLineId, nXSize );
		else
			oFirstEnum.ProcessLine(
														 panLastLineVal, panThisLineVal,
														 panLastLineId,  panThisLineId,
														 nXSize );

		// swap lines
		int *panTmpVal = panLastLineVal;
		panLastLineVal = panThisLineVal;
		panThisLineVal = panTmpVal;

		GInt32* panTmp = panThisLineId;
		panThisLineId = panLastLineId;
		panLastLineId = panTmp;
	}

	/* -------------------------------------------------------------------- */
	/*      Make a pass through the maps, ensuring every polygon id         */
	/*      points to the final id it should use, not an intermediate       */
	/*      value.                                                          */
	/* -------------------------------------------------------------------- */
	oFirstEnum.CompleteMerges();

	/* -------------------------------------------------------------------- */
	/*      Initialize ids to -1 to serve as a nodata value for the         */
	/*      previous line, and past the beginning and end of the            */
	/*      scanlines.                                                      */
	/* -------------------------------------------------------------------- */
	int iX;

	panThisLineId[0] = -1;
	panThisLineId[nXSize+1] = -1;

	for( iX = 0; iX < nXSize+2; iX++ )
		panLastLineId[iX] = -1;

	/* -------------------------------------------------------------------- */
	/*      We will use a new enumerator for the second pass primarily      */
	/*      so we can preserve the first pass map.                          */
	/* -------------------------------------------------------------------- */
	GDALRasterPolygonEnumeratorT<int, IntEqualityTest> oSecondEnum(nConnectedness);

	RPolygon **papoPoly = (RPolygon **)
		CPLCalloc(sizeof(RPolygon*),oFirstEnum.nNextPolygonId);

	/* ==================================================================== */
	/*      Second pass during which we will actually collect polygon       */
	/*      edges as geometries.                                            */
	/* ==================================================================== */
	for( iY = 0; Res && iY < nYSize+1; iY++ )
	{
		/* -------------------------------------------------------------------- */
		/*      Read the image data.                                            */
		/* -------------------------------------------------------------------- */
		if( iY < nYSize )
		{
			for (int iX = 0; iX < nXSize; iX++)
				if (img_h(iX,iY) > 0)
					panThisLineVal[iX] = img_h(iX,iY);
				else
					panThisLineVal[iX] = GP_NODATA_MARKER;
		}

		/* -------------------------------------------------------------------- */
		/*      Determine what polygon the various pixels belong to (redoing    */
		/*      the same thing done in the first pass above).                   */
		/* -------------------------------------------------------------------- */
		if( iY == nYSize )
		{
			for( iX = 0; iX < nXSize+2; iX++ )
				panThisLineId[iX] = -1;
		}
		else if( iY == 0 )
			oSecondEnum.ProcessLine(
															NULL, panThisLineVal, NULL, panThisLineId+1, nXSize );
		else
			oSecondEnum.ProcessLine(
															panLastLineVal, panThisLineVal,
															panLastLineId+1,  panThisLineId+1,
															nXSize );

		/* -------------------------------------------------------------------- */
		/*      Add polygon edges to our polygon list for the pixel             */
		/*      boundaries within and above this line.                          */
		/* -------------------------------------------------------------------- */
		for( iX = 0; iX < nXSize+1; iX++ )
		{
			AddEdges( panThisLineId, panLastLineId,
								oFirstEnum.panPolyIdMap, oFirstEnum.panPolyValue,
								papoPoly, iX, iY );
		}

		/* -------------------------------------------------------------------- */
		/*      Periodically we scan out polygons and write out those that      */
		/*      haven't been added to on the last line as we can be sure        */
		/*      they are complete.                                              */
		/* -------------------------------------------------------------------- */
		if( iY % 8 == 7 )
		{
			for( iX = 0; Res && iX < oSecondEnum.nNextPolygonId; iX++ )
			{
				if( papoPoly[iX] && papoPoly[iX]->nLastLineUpdated < iY-1 )
				{
					Res = EmitPolygonToLayer(hOutLayer, papoPoly[iX]);
					delete papoPoly[iX];
					papoPoly[iX] = NULL;
				}
			}
		}

		/* -------------------------------------------------------------------- */
		/*      Swap pixel value, and polygon id lines to be ready for the      */
		/*      next line.                                                      */
		/* -------------------------------------------------------------------- */
		int *panTmpVal = panLastLineVal;
		panLastLineVal = panThisLineVal;
		panThisLineVal = panTmpVal;

		GInt32* panTmp = panThisLineId;
		panThisLineId = panLastLineId;
		panLastLineId = panTmp;
	}

	/* -------------------------------------------------------------------- */
	/*      Make a cleanup pass for all unflushed polygons.                 */
	/* -------------------------------------------------------------------- */
	for( iX = 0; Res && iX < oSecondEnum.nNextPolygonId; iX++ )
  {
		if( papoPoly[iX] )
		{
			Res = EmitPolygonToLayer(hOutLayer, papoPoly[iX]);
			delete papoPoly[iX];
			papoPoly[iX] = NULL;
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Cleanup                                                         */
	/* -------------------------------------------------------------------- */
	CPLFree( panThisLineId );
	CPLFree( panLastLineId );
	CPLFree( panThisLineVal );
	CPLFree( panLastLineVal );
	CPLFree( pabyMaskLine );
	CPLFree( papoPoly );

	return Res;
}

