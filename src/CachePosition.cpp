#include "CachePosition.h"
#include <cstring>

//  numVertex : number of vertex of this Polygon
//  num_rot : possible number of rotations
//  numPLacedVertex : number of vertex of all polygons before this one

CachePosition::CachePosition(int numVertex, int num_rot, int numPlacedVertex)
{
    sizeCacheEntry1 = numVertex;                           //  Used later to compute offsets
    sizeCacheEntry2 = numVertex*numPlacedVertex;           //  Used later to compute offsets
    currentIdxCache = -1;                                  //  Init state, nothing in cache
    memCache = numPlacedVertex * num_rot * numVertex;       //  Memory claim
    CachePlaced = new Polygon *[memCache];
    memset(CachePlaced, 0, memCache*sizeof(Polygon *));     //  Init all entries, such as impossible !
    FixedPerRot = new int[num_rot];
    memset(FixedPerRot, -1, num_rot*sizeof(int));             //  Init all entries, such as not already computed !
}

CachePosition::~CachePosition()
{
    for ( int i = 0; i < currentIdxCache; i++ )
    {
        if ( CachePlaced[i] != NULL) delete CachePlaced[i];
    }
    delete [] CachePlaced;
    delete [] FixedPerRot;
}

//  Return status of position polygon such as Polygon is Rotated iRot and vertex #iVertex is placed on Placed Polygon vertex iPlacedVertex
//  Return -1 if impossible, > 0 OK and 0 don't know
//  If OK, update PlacedPoly parameter with the vertexes placed at the right position
//  If OK return the current value of MAxFixedPoly +1

int CachePosition::isOKPlaced(int iVertex, int iRot, int iPlacedVertex, Polygon **PlacedPoly)
{
int idx = iRot * sizeCacheEntry2 + iPlacedVertex * sizeCacheEntry1 + iVertex;

    if ( idx > currentIdxCache )
        return 0;       //  Not yet computed
    if ( iPlacedVertex > FixedPerRot[iRot])
        return 0;       //  Not yet computed
    if ( CachePlaced[idx] != NULL )
    {
        *PlacedPoly = CachePlaced[idx];
        return FixedPerRot[iRot] + 1;
    }
    return -1;
}

//  Add the polygon corresponding to the position where Polygon Poly is rotated with an index iRot and placed on vertex iPlacedVertex of the fixed polygons list
//  If this entry is invalid, Poly should be NULL

void CachePosition::addOKPlaced(int iVertex, int iRot, int iPlacedVertex, Polygon *Poly)
{
    if ( iPlacedVertex > FixedPerRot[iRot] )
        FixedPerRot[iRot] = iPlacedVertex;
int idx = iRot * sizeCacheEntry2 + iPlacedVertex * sizeCacheEntry1 + iVertex;
    if ( currentIdxCache < idx)
        currentIdxCache = idx;
    CachePlaced[idx] = Poly;
    if ( Poly != NULL )
    {
        Poly->TypePoly = CachedPoly;
        nbCreatedCachedPoly++;
        nbCachedPoly++;
        nbTranslatedPoly--;
    }
}

//  Return status of position polygon such as Polygon is Rotated iRot and vertex #iVertex is placed on Placed Polygon vertex iPlacedVertex
//  Return -1 if impossible, > 0 OK and 0 don't know
//  If OK, update PlacedPoly parameter with the vertexes placed at the right position
//  If OK return the current value of MAxFixedPoly +1

int CachePosition::isOKPlacedFreeRot(int iVertex, int iRot, int iPlacedVertex, Polygon **PlacedPoly)
{
int idx = (iPlacedVertex * sizeCacheEntry1 + iVertex)*4 + iRot;

    if ( idx > currentIdxCache )
        return 0;       //  Not yet computed
    if ( CachePlaced[idx] != NULL )
    {
        *PlacedPoly = CachePlaced[idx];
        return 1;
    }
    return -1;
}

//  Add the polygon corresponding to the position where Polygon Poly is rotated with an index iRot and placed on vertex iPlacedVertex of the fixed polygons list
//  If this entry is invalid, Poly should be NULL

void CachePosition::addOKPlacedFreeRot(int iVertex, int iRot, int iPlacedVertex, Polygon *Poly)
{
int idx = ( iPlacedVertex * sizeCacheEntry1 + iVertex) * 4 + iRot;
    if ( currentIdxCache < idx)
        currentIdxCache = idx;
    CachePlaced[idx] = Poly;
    if ( Poly != NULL )
    {
        Poly->TypePoly = CachedPoly;
        nbCreatedCachedPoly++;
        nbCachedPoly++;
        nbTranslatedPoly--;
    }
}
