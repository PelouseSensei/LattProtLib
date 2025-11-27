#ifndef ROTATION_SQUARE_H
#define ROTATION_SQUARE_H

//---------------------------- Rotation Square ---------------------------
class RotationSquare : public Move<SquareLattice> {
public:
    RotationSquare(int seed) : Move<SquareLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.dE = 0.0;
        moveArgs.rotation_angle = 0;
    }

    bool apply(Configuration<SquareLattice>& config) override;

    void generateMoveArgs(Configuration<SquareLattice>& config) override;

    void checkMove(Configuration<SquareLattice>& config) override;

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
    RotMoveArgs<SquareLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_angle{0, 3};     // rotation angle selection
    std::uniform_int_distribution<int> dist_mono{0, 0};   // pivot point offset
};

#include "rotation_square.tpp"

#endif // ROTATION_SQUARE_H
