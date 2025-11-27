#ifndef PULL_CUBIC_H
#define PULL_CUBIC_H

#include "../base_move.hpp"
#include "../custom_structs.hpp"

//-------------------------- Pull move -----------------------------



template<bool Forward>
class PullCubic : public Move<CubicLattice> {
public:
    PullCubic(int seed) : Move<CubicLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.num_monomer_moved = 0;
        moveArgs.dE = 0.0;

        dist_backbone = std::uniform_int_distribution<int>{0, 0};
        dist_start = std::uniform_int_distribution<int>{0, 0};
        dist_start_backward = std::uniform_int_distribution<int>{0, 0};
    }
    

    // Core move functions
    void checkMove(Configuration<CubicLattice>& config);
    bool apply(Configuration<CubicLattice>& config);
    void generateMoveArgs(Configuration<CubicLattice>& config);

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
    PullMoveArgs<CubicLattice> moveArgs;

    std::uniform_int_distribution<int> dist_backbone;
    std::uniform_int_distribution<int> dist_start;
    std::uniform_int_distribution<int> dist_perp{0,3};
    std::uniform_int_distribution<int> dist_start_backward;
    
    // For each axis-aligned unit vector, store its 4 perpendicular moves
	const Position perp_table[6][4] = {
		// dx=+1
		{ {0,  1, 0}, {0, -1, 0}, {0, 0,  1}, {0, 0, -1} },
		// dx=-1
		{ {0,  1, 0}, {0, -1, 0}, {0, 0,  1}, {0, 0, -1} },
		// dy=+1
		{ { 1, 0, 0}, {-1, 0, 0}, {0, 0,  1}, {0, 0, -1} },
		// dy=-1
		{ { 1, 0, 0}, {-1, 0, 0}, {0, 0,  1}, {0, 0, -1} },
		// dz=+1
		{ { 1, 0, 0}, {-1, 0, 0}, {0,  1, 0}, {0, -1, 0} },
		// dz=-1
		{ { 1, 0, 0}, {-1, 0, 0}, {0,  1, 0}, {0, -1, 0} }
	};
};

// Aliases for convenience
using PullForwardCubic = PullCubic<true>;
using PullBackwardCubic = PullCubic<false>;

#include "pull_cubic.tpp"

#endif // PULL_CUBIC_H
