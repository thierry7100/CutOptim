#ifndef CACHEPOSITION_H
#define CACHEPOSITION_H

#include "Geometry.h"

class CachePosition
{
    public:
        CachePosition(int numVertex, int num_rot, int numPlacedVertex);
        virtual ~CachePosition();
        int isOKPlaced(int iVertex, int iRot, int iPlacedVertex, Polygon **PlacedPoly);
        void addOKPlaced(int iVertex, int iRot, int iPlacedVertex, Polygon *Poly);
        int isOKPlacedFreeRot(int iVertex, int iRot, int iPlacedVertex, Polygon **PlacedPoly);
        void addOKPlacedFreeRot(int iVertex, int iRot, int iPlacedVertex, Polygon *Poly);
    protected:
        size_t sizeCacheEntry1;
        size_t sizeCacheEntry2;
        size_t memCache;
        int currentIdxCache;
        Polygon **CachePlaced;          //  Polygon cached, if impossible set to 0
        int *FixedPerRot;               //  Index of cached value per rotation. Used when vertices are added (new polygons)

    private:
};

#endif // CACHEPOSITION_H
