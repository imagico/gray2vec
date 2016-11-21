/* @(#)gdal_polygonize_mod.cpp

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

#include "gdal_polygonize_mod.h"

/************************************************************************/
/*                                Dump()                                */
/************************************************************************/
void RPolygon::Dump()
{
    size_t iString;

    printf( "RPolygon: Value=%g, LastLineUpdated=%d\n",
            dfPolyValue, nLastLineUpdated );

    for( iString = 0; iString < aanXY.size(); iString++ )
    {
        std::vector<int> &anString = aanXY[iString];
        size_t iVert;

        printf( "  String %d:\n", (int) iString );
        for( iVert = 0; iVert < anString.size(); iVert += 2 )
        {
            printf( "    (%d,%d)\n", anString[iVert], anString[iVert+1] );
        }
    }
}

/************************************************************************/
/*                              Coalesce()                              */
/************************************************************************/

void RPolygon::Coalesce()

{
    size_t iBaseString;

/* -------------------------------------------------------------------- */
/*      Iterate over loops starting from the first, trying to merge     */
/*      other segments into them.                                       */
/* -------------------------------------------------------------------- */
    for( iBaseString = 0; iBaseString < aanXY.size(); iBaseString++ )
    {
        std::vector<int> &anBase = aanXY[iBaseString];
        int bMergeHappened = TRUE;

/* -------------------------------------------------------------------- */
/*      Keep trying to merge the following strings into our target      */
/*      "base" string till we have tried them all once without any      */
/*      mergers.                                                        */
/* -------------------------------------------------------------------- */
        while( bMergeHappened )
        {
            size_t iString;

            bMergeHappened = FALSE;

/* -------------------------------------------------------------------- */
/*      Loop over the following strings, trying to find one we can      */
/*      merge onto the end of our base string.                          */
/* -------------------------------------------------------------------- */
            for( iString = iBaseString+1;
                 iString < aanXY.size();
                 iString++ )
            {
                std::vector<int> &anString = aanXY[iString];

                if( anBase[anBase.size()-2] == anString[0]
                    && anBase[anBase.size()-1] == anString[1] )
                {
                    Merge( static_cast<int>(iBaseString), static_cast<int>(iString), 1 );
                    bMergeHappened = TRUE;
                }
                else if( anBase[anBase.size()-2] == anString[anString.size()-2]
                         && anBase[anBase.size()-1] == anString[anString.size()-1] )
                {
                    Merge( static_cast<int>(iBaseString), static_cast<int>(iString), -1 );
                    bMergeHappened = TRUE;
                }
            }
        }

        /* At this point our loop *should* be closed! */

        CPLAssert( anBase[0] == anBase[anBase.size()-2]
                   && anBase[1] == anBase[anBase.size()-1] );
    }

}

/************************************************************************/
/*                               Merge()                                */
/************************************************************************/

void RPolygon::Merge( int iBaseString, int iSrcString, int iDirection )

{
    std::vector<int> &anBase = aanXY[iBaseString];
    std::vector<int> &anString = aanXY[iSrcString];
    int iStart, iEnd, i;

    if( iDirection == 1 )
    {
        iStart = 1;
        iEnd = static_cast<int>(anString.size()) / 2;
    }
    else
    {
        iStart = static_cast<int>(anString.size()) / 2 - 2;
        iEnd = -1;
    }

    for( i = iStart; i != iEnd; i += iDirection )
    {
        anBase.push_back( anString[i*2+0] );
        anBase.push_back( anString[i*2+1] );
    }

    if( iSrcString < ((int) aanXY.size())-1 )
        aanXY[iSrcString] = aanXY[aanXY.size()-1];

    size_t nSize = aanXY.size();
    aanXY.resize(nSize-1);
}

/************************************************************************/
/*                             AddSegment()                             */
/************************************************************************/

void RPolygon::AddSegment( int x1, int y1, int x2, int y2 )

{
    nLastLineUpdated = MAX(y1, y2);

/* -------------------------------------------------------------------- */
/*      Is there an existing string ending with this?                   */
/* -------------------------------------------------------------------- */
    size_t iString;

    for( iString = 0; iString < aanXY.size(); iString++ )
    {
        std::vector<int> &anString = aanXY[iString];
        size_t nSSize = anString.size();

        if( anString[nSSize-2] == x1
            && anString[nSSize-1] == y1 )
        {
            int nTemp;

            nTemp = x2;
            x2 = x1;
            x1 = nTemp;

            nTemp = y2;
            y2 = y1;
            y1 = nTemp;
        }

        if( anString[nSSize-2] == x2
            && anString[nSSize-1] == y2 )
        {
            // We are going to add a segment, but should we just extend
            // an existing segment already going in the right direction?

            int nLastLen = MAX(ABS(anString[nSSize-4]-anString[nSSize-2]),
                               ABS(anString[nSSize-3]-anString[nSSize-1]));
						/*
            if( nSSize >= 4
                && (anString[nSSize-4] - anString[nSSize-2]
                    == (anString[nSSize-2] - x1)*nLastLen)
                && (anString[nSSize-3] - anString[nSSize-1]
                    == (anString[nSSize-1] - y1)*nLastLen) )
            {
                anString.pop_back();
                anString.pop_back();
            }
						*/
            anString.push_back( x1 );
            anString.push_back( y1 );
            return;
        }
    }

/* -------------------------------------------------------------------- */
/*      Create a new string.                                            */
/* -------------------------------------------------------------------- */
    size_t nSize = aanXY.size();
    aanXY.resize(nSize + 1);
    std::vector<int> &anString = aanXY[nSize];

    anString.push_back( x1 );
    anString.push_back( y1 );
    anString.push_back( x2 );
    anString.push_back( y2 );

    return;
}

/************************************************************************/
/* ==================================================================== */
/*     End of RPolygon                                                  */
/* ==================================================================== */
/************************************************************************/

