/* @(#)gdal_polygonize_mod.h

  This code is derived from polygonize.cpp from the gdal source package

 * $Id: polygonize.cpp 33757 2016-03-20 20:22:33Z goatbar $
 * Project:  GDAL
 * Purpose:  Raster to Polygon Converter
 * Author:   Frank Warmerdam, warmerdam@pobox.com
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

#ifndef _gdal_polygonize_mod_H
#define _gdal_polygonize_mod_H

#include <gdal_alg_priv.h>
#include <cpl_conv.h>
#include <cpl_string.h>
#include <vector>


#ifndef GP_NODATA_MARKER
    #define GP_NODATA_MARKER -51502112
#endif

// we use the newer version of this from GDAL 2.x in case of GDAL 1.x
#if GDAL_VERSION_MAJOR < 2

struct IntEqualityTest
{
    bool operator()(GInt32 a, GInt32 b) { return a == b; }
};

template<class DataType, class EqualityTest> class GDALRasterPolygonEnumeratorT
{
private:
    void     MergePolygon( int nSrcId, int nDstId );
    int      NewPolygon( DataType nValue );

public:  // these are intended to be readonly.

    GInt32   *panPolyIdMap;
    DataType   *panPolyValue;

    int      nNextPolygonId;
    int      nPolyAlloc;

    int      nConnectedness;

public:
             GDALRasterPolygonEnumeratorT( int nConnectedness=4 );
            ~GDALRasterPolygonEnumeratorT();

    void     ProcessLine( DataType *panLastLineVal, DataType *panThisLineVal,
                          GInt32 *panLastLineId,  GInt32 *panThisLineId,
                          int nXSize );

    void     CompleteMerges();

    void     Clear();
};
#endif

/************************************************************************/
/* ==================================================================== */
/*                               RPolygon                               */
/*									*/
/*      This is a helper class to hold polygons while they are being    */
/*      formed in memory, and to provide services to coalesce a much    */
/*      of edge sections into complete rings.                           */
/* ==================================================================== */
/************************************************************************/

class RPolygon {
public:
    RPolygon(  double dfValue ) { dfPolyValue = dfValue; nLastLineUpdated = -1; }

    double              dfPolyValue;
    int              nLastLineUpdated;

    std::vector< std::vector<int> > aanXY;

    void             AddSegment( int x1, int y1, int x2, int y2 );
    void             Dump();
    void             Coalesce();
    void             Merge( int iBaseString, int iSrcString, int iDirection );
};


/************************************************************************/
/*                              AddEdges()                              */
/*                                                                      */
/*      Examine one pixel and compare to its neighbour above            */
/*      (previous) and right.  If they are different polygon ids        */
/*      then add the pixel edge to this polygon and the one on the      */
/*      other side of the edge.                                         */
/************************************************************************/

template<class DataType>
static void AddEdges( GInt32 *panThisLineId, GInt32 *panLastLineId,
                      GInt32 *panPolyIdMap, DataType *panPolyValue,
                      RPolygon **papoPoly, int iX, int iY )

{
    int nThisId = panThisLineId[iX];
    int nRightId = panThisLineId[iX+1];
    int nPreviousId = panLastLineId[iX];
    int iXReal = iX - 1;

    if( nThisId != -1 )
        nThisId = panPolyIdMap[nThisId];
    if( nRightId != -1 )
        nRightId = panPolyIdMap[nRightId];
    if( nPreviousId != -1 )
        nPreviousId = panPolyIdMap[nPreviousId];

    if( nThisId != nPreviousId )
    {
        if( nThisId != -1 )
        {
            if( papoPoly[nThisId] == NULL )
                papoPoly[nThisId] = new RPolygon( panPolyValue[nThisId] );

            papoPoly[nThisId]->AddSegment( iXReal, iY, iXReal+1, iY );
        }
        if( nPreviousId != -1 )
        {
            if( papoPoly[nPreviousId] == NULL )
                papoPoly[nPreviousId] = new RPolygon(panPolyValue[nPreviousId]);

            papoPoly[nPreviousId]->AddSegment( iXReal, iY, iXReal+1, iY );
        }
    }

    if( nThisId != nRightId )
    {
        if( nThisId != -1 )
        {
            if( papoPoly[nThisId] == NULL )
                papoPoly[nThisId] = new RPolygon(panPolyValue[nThisId]);

            papoPoly[nThisId]->AddSegment( iXReal+1, iY, iXReal+1, iY+1 );
        }

        if( nRightId != -1 )
        {
            if( papoPoly[nRightId] == NULL )
                papoPoly[nRightId] = new RPolygon(panPolyValue[nRightId]);

            papoPoly[nRightId]->AddSegment( iXReal+1, iY, iXReal+1, iY+1 );
        }
    }
}

/************************************************************************/
/*                          GPMaskImageData()                           */
/*                                                                      */
/*      Mask out image pixels to a special nodata value if the mask     */
/*      band is zero.                                                   */
/************************************************************************/

template<class DataType>
static CPLErr
GPMaskImageData( GDALRasterBandH hMaskBand, GByte* pabyMaskLine, int iY, int nXSize,
                 DataType *panImageLine )

{
    CPLErr eErr;

    eErr = GDALRasterIO( hMaskBand, GF_Read, 0, iY, nXSize, 1,
                         pabyMaskLine, nXSize, 1, GDT_Byte, 0, 0 );
    if( eErr == CE_None )
    {
        int i;
        for( i = 0; i < nXSize; i++ )
        {
            if( pabyMaskLine[i] == 0 )
                panImageLine[i] = GP_NODATA_MARKER;
        }
    }

    return eErr;
}

#endif /* _gdal_polygonize_mod_H */
