#ifndef ROTATION_CUBIC_H
#define ROTATION_CUBIC_H

//---------------------------- Rotation Cubic ---------------------------
class RotationCubic : public Move<CubicLattice> {
public:
    RotationCubic(int seed) : Move<CubicLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.dE = 0.0;
        moveArgs.rotation_angle = 0;
        moveArgs.rotation_axis = RotationAxis::X;
    }

    bool apply(Configuration<CubicLattice>& config) override;

    void generateMoveArgs(Configuration<CubicLattice>& config) override;

    void checkMove(Configuration<CubicLattice>& config) override;

    double getdE() { return moveArgs.dE; }

    MoveStatus getStatus() { return moveArgs.status; }

    void setBackboneDist(const int b_min, const int b_max) {
        dist_backbone = std::uniform_int_distribution<int>{b_min, b_max};
    }

    void setAngleRange(const int min_angle, const int max_angle) {
        dist_angle = std::uniform_int_distribution<int>{min_angle, max_angle};
    }

    void setDistMono(const int min, const int max) {
        dist_mono = std::uniform_int_distribution<int>{min, max};
    }

private:
    RotMoveArgs3D<CubicLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_angle{0, 3};     // rotation angle selection
    std::uniform_int_distribution<int> dist_axis{0, 2};     // rotation angle selection
    std::uniform_int_distribution<int> dist_mono{0, 0};   // pivot point offset
};

#include "rotation_cubic.tpp"

#endif // ROTATION_CUBIC_H
