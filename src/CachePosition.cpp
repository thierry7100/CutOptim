#include "CachePosition.h"
#include <cstring>

/*
*   This class implements the caching used to speed of computations.
*   When the level of optimization is greater than 2 (and only if greater than 2), we will eventually recompute if a current floating polygon
*   intersects the fixed ones.
*   So we will store the result to skip the next computations.
*   Beware, only computations against the fixed polygons are meaningful
*   This will happen when LevelOptimize is at least 2, because when LevelOptimize is 1, the floating polygon will move or rotate at each iteration
*   But with level 2, when there are for example 5 fixed polygons. We will try all positions and rotations for polygon 6 and then all positions
*   and rotations for polygon 7. But we will test polygon 7 positioned on vertices belonging to fixed polygons for each position/rotation of
*   polygon 6. And if there is an impossibility with polygon 7 and vertex X of fixed polygon N, this will always be the case.
*/

//  Constructor.
//  The object is built once (for each floating polygon) and remains active for all iterations for the placement of the biggest remaining polygon
//  It stays alive until this polygon is placed and become fixed.
//  This save time, as when an entry is impossible with X fixed polygons it is also impossible with X+1 fixed polygons
//  numVertex : number of vertex of this Polygon
//  num_rot : possible number of rotations
//  numPLacedVertex : number of vertex of all fixed polygons before this one
//  Allocate a (large) array of pointers to Polygon objects
//  The size of the array is numPlacedVertex * num_rot * numVertex, one entry per vertex of current polygon with every rotation and every vertices already placed.
//  This array will be accessed as memCache[idxRot][iPlacedVertex][numVertex]

CachePosition::CachePosition(int numVertex, int num_rot, int numPlacedVertex)
{
    sizeCacheEntry1 = numVertex;                           //  Used later to compute offsets for iPlacedVertex so multiply by numVertex
    sizeCacheEntry2 = numVertex*numPlacedVertex;           //  Used later to compute offsets for idxRot so multiply by numVertex*numPlacedVertex
    currentIdxCache = -1;                                  //  Init state, nothing in cache
    memCache = numPlacedVertex * num_rot * numVertex;      //  Memory claim
    CachePlaced = new Polygon *[memCache];
    //  The fake value  0xCACACACACACACA is impossible as a pointer
    //  It will be used to detect entries which are not yet processed.
    memset(CachePlaced, 0xCA, memCache*sizeof(Polygon *));    //  Init all entries, such as impossible !
    FixedPerRot = new int[num_rot];
    memset(FixedPerRot, -1, num_rot*sizeof(int));          //  Init all entries, such as not already computed !
}

//  Just release the memory

CachePosition::~CachePosition()
{
    for ( int i = 0; i < currentIdxCache; i++ )
    {
        uint64_t val = (uint64_t )CachePlaced[i];       //  default value, skip, DO NOT try to free this fake memory pointer
        if ( val == 0xCACACACACACACACA )
            continue;
        if ( CachePlaced[i] != NULL) delete CachePlaced[i];
    }
    delete [] CachePlaced;
    delete [] FixedPerRot;
}

//  Return status of position polygon such as Polygon is Rotated iRot and vertex #iVertex is placed on Placed Polygon vertex iPlacedVertex
//  Return -1 if impossible, > 0 OK and 0 don't know
//  If OK, update PlacedPoly parameter with the vertexes placed at the right position
//  If OK return the current value of MaxFixedPoly +1


int CachePosition::isOKPlaced(int iVertex, int iRot, int iPlacedVertex, Polygon **PlacedPoly)
{
int idx = iRot * sizeCacheEntry2 + iPlacedVertex * sizeCacheEntry1 + iVertex;

    if ( idx > currentIdxCache )
        return 0;       //  Not yet computed, return don't know
    if ( iPlacedVertex > FixedPerRot[iRot])
        return 0;       //  Not yet computed
    uint64_t val = (uint64_t )CachePlaced[idx];
    if ( val == 0xCACACACACACACACA )            //  Check if default value (not yet processed) ?
    {
        return 0;           //  If so, return don't know
    }
    if ( CachePlaced[idx] != NULL )             //  Value not null, cache is OK
    {
        *PlacedPoly = CachePlaced[idx];         //  Update pointer to polygon to avoid recompute
        return FixedPerRot[iRot] + 1;           //  Return positive value
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
