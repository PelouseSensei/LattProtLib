#ifndef PULL_END_SQUARE_H
#define PULL_END_SQUARE_H

#include "../base_move.hpp"
#include "../custom_structs.hpp"

//-------------------------- Pull End move -----------------------------
class PullEndSquare : public Move<SquareLattice> {
public:
    PullEndSquare(int seed) : Move<SquareLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.num_monomer_moved = 0;
        moveArgs.dE = 0.0;
        forward = true;
        b_length = 0;
    }

    bool apply(Configuration<SquareLattice>& config) override;
    void generateMoveArgs(Configuration<SquareLattice>& config) override;
    void checkMove(Configuration<SquareLattice>& config) override;
    
    double getdE() { return moveArgs.dE; }
    MoveStatus getStatus() { return moveArgs.status; }
    
    void setBackboneDist(const int b_min, const int b_max) {
        dist_backbone = std::uniform_int_distribution<int>{b_min, b_max};
    }

private:
    PullMoveArgs<SquareLattice> moveArgs;
    bool forward;
    int b_length;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_end{0, 1}; // choose 0 or N-1
};

#include "pull_end_square.tpp"

#endif // PULL_END_SQUARE_H
