#ifndef SYMMETRY_SQUARE_H
#define SYMMETRY_SQUARE_H

//---------------------------- Symmetry Square ---------------------------
class SymmetrySquare : public Move<SquareLattice> {
public:
    SymmetrySquare(int seed) : Move<SquareLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.dE = 0.0;
        moveArgs.symmetry_axis = SymmetryAxis::None;
    }

    bool apply(Configuration<SquareLattice>& config) override;

    void generateMoveArgs(Configuration<SquareLattice>& config) override;

    void checkMove(Configuration<SquareLattice>& config) override;

    double getdE() { return moveArgs.dE; }

    MoveStatus getStatus() { return moveArgs.status; }

    void setBackboneDist(const int b_min, const int b_max) {
        dist_backbone = std::uniform_int_distribution<int>{b_min, b_max};
    }

    void setAxisRange(const int min_axis, const int max_axis) {
        dist_axis = std::uniform_int_distribution<int>{min_axis, max_axis};
    }

private:
    SymMoveArgs<SquareLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_axis{1, 5};   // symmetry axis selection
};

#include "symmetry_square.tpp"

#endif // SYMMETRY_SQUARE_H
