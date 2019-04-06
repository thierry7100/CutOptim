#ifndef PLACEDPIECE_H
#define PLACEDPIECE_H


class PlacedPiece
{
    public:
        PlacedPiece(int nRot, int nEdges);
        virtual ~PlacedPiece();

        int isPlaceOK(iRot, iEdge, iPoly, jRot, jEdge, jPoly);

    protected:
        Polygon **PlacedPoly;           //  Array used to store the placed polygon with rotation and translation.
        int numRot;
        int numEdges;
        size_t size_PlacedPoly;
    private:
};

#endif // PLACEDPIECE_H
