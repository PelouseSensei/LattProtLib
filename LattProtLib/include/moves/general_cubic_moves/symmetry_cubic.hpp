#ifndef SYMMETRY_CUBIC_H
#define SYMMETRY_CUBIC_H

//---------------------------- Symmetry Cubic ---------------------------
class SymmetryCubic : public Move<CubicLattice> {
public:
    SymmetryCubic(int seed) : Move<CubicLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.dE = 0.0;
        moveArgs.symmetry_axis = SymmetryAxis3D::None;
    }

    bool apply(Configuration<CubicLattice>& config) override;

    void generateMoveArgs(Configuration<CubicLattice>& config) override;

    void checkMove(Configuration<CubicLattice>& config) override;

    double getdE() { return moveArgs.dE; }

    MoveStatus getStatus() { return moveArgs.status; }

    void setBackboneDist(const int b_min, const int b_max) {
        dist_backbone = std::uniform_int_distribution<int>{b_min, b_max};
    }

    void setAxisRange(const int min_axis, const int max_axis) {
        dist_axis = std::uniform_int_distribution<int>{min_axis, max_axis};
    }

private:
    SymMoveArgs3D<CubicLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_axis{1, 5};   // symmetry axis selection
};

#include "symmetry_cubic.tpp"

#endif // SYMMETRY_CUBIC_H
