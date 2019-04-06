#ifndef SAVEPLACEDPIECE_H
#define SAVEPLACEDPIECE_H

#include "Geometry.h"

struct PiecePosition {
    unsigned short idx_Rot1;
    unsigned short idx_Vertex1;
    unsigned short idxPoly2;
    unsigned short idx_Rot2;
    unsigned short idx_Vertex2;
};


class SavePlacedPiece
{
    public:
        SavePlacedPiece(Polygon *Poly, int idxPoly, int numVertex);
        virtual ~SavePlacedPiece();
        int isOK(int idx_rot1, int idx_vertex1, int idxPoly2, int idx_rot2, int idx_vertex2);

    protected:
        unsigned char CachedOK

    private:
};

#endif // SAVEPLACEDPIECE_H
