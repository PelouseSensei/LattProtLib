#ifndef PULL_SQUARE_H
#define PULL_SQUARE_H

#include "../base_move.hpp"
#include "../custom_structs.hpp"

//-------------------------- Pull move -----------------------------
template<bool Forward>
class PullSquare : public Move<SquareLattice> {
public:
    PullSquare(int seed) : Move<SquareLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.num_monomer_moved = 0;
        moveArgs.dE = 0.0;

        dist_backbone = std::uniform_int_distribution<int>{0, 0};
        dist_start = std::uniform_int_distribution<int>{0, 0};
        dist_start_backward = std::uniform_int_distribution<int>{0, 0};
    }

    // Core move functions
    void checkMove(Configuration<SquareLattice>& config);
    bool apply(Configuration<SquareLattice>& config);
    void generateMoveArgs(Configuration<SquareLattice>& config);

    // Setter functions for distributions
    void setBackboneDist(int min, int max) {
        dist_backbone = std::uniform_int_distribution<int>(min, max);
    }

    void setMonomerDist(int min, int max) {
        dist_start = std::uniform_int_distribution<int>(min, max);
        dist_start_backward = std::uniform_int_distribution<int>(min, max);
    }

    double getdE() { return moveArgs.dE; }
    MoveStatus getStatus() { return moveArgs.status; }

private:
    PullMoveArgs<SquareLattice> moveArgs;

    std::mt19937 rng;
    std::uniform_int_distribution<int> dist_backbone;
    std::uniform_int_distribution<int> dist_start;
    std::uniform_int_distribution<int> dist_start_backward;
};

// Aliases for convenience
using PullForwardSquare = PullSquare<true>;
using PullBackwardSquare = PullSquare<false>;

#include "pull_square.tpp"

#endif // PULL_SQUARE_H
