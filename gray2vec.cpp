/* ========================================================================
    File: @(#)gray2vec.cpp
   ------------------------------------------------------------------------
    vectorizes grayscale image into polygons
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

const char PROGRAM_TITLE[] = "gray2vec version 0.1";

#include <cstdlib>
#include <algorithm>
#include <stack>
#include <cmath>

#include "Gray2Vec_Grid.h"

int main(int argc,char **argv)
{
	std::fprintf(stderr,"%s\n", PROGRAM_TITLE);
	std::fprintf(stderr,"-------------------------------------------------------\n");
	std::fprintf(stderr,"Copyright (C) 2016 Christoph Hormann\n");
	std::fprintf(stderr,"This program comes with ABSOLUTELY NO WARRANTY;\n");
	std::fprintf(stderr,"This is free software, and you are welcome to redistribute\n");
	std::fprintf(stderr,"it under certain conditions; see COPYING for details.\n");

	cimg_usage("Usage: gray2vec [options]");

	// --- Read command line parameters ---

	// Files
	const std::string file_o = cimg_option("-o","","output vector file");
	const std::string file_i = cimg_option("-i","","input image");
	const std::string file_c = cimg_option("-c","","combined input image");

	const std::string Layer = cimg_option("-l","polygons","output layer");

	const int Xc = cimg_option("-x",-1,"x attribute to apply to generated polygons");
	const int Yc = cimg_option("-y",-1,"y attribute to apply to generated polygons");
	const int Zc = cimg_option("-z",-1,"z attribute to apply to generated polygons");

	const bool Complement = cimg_option("-complement",false,"process complement of input");

	const bool Append = cimg_option("-append",false,"append to existing vector file");

	const double MaxError = cimg_option("-me",0.05,"maximum error to accept for pixel coverage fraction");

	const bool Debug = cimg_option("-debug",false,"generate debug output");

	if (file_i.empty() || file_o.empty())
	{
		std::fprintf(stderr,"You must specify input and output files (try '%s -h').\n\n",argv[0]);
		std::exit(1);
	}

	Gray2Vec_Grid g2v(file_i, file_c, Complement, Debug);

	if ((Xc >= 0) || (Yc >= 0) || (Zc >= 0))
		g2v.SetAttributes(Xc, Yc, Zc);

	g2v.Analyze();
	g2v.NeighborsAdjust();
	g2v.ResolveConflicts();
	g2v.InitFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.NeighborsAdjust2();
	g2v.ResolveConflicts();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.NeighborsAdjust2();
	g2v.ResolveConflicts();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions();
	g2v.FractionsNeighborsAdj();
	g2v.TuneFractions(MaxError);

	g2v.Vectorize(file_o, Layer, Append);
}
