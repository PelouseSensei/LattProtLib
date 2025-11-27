#ifndef PULL_END_CUBIC_H
#define PULL_END_CUBIC_H

#include "../base_move.hpp"
#include "../custom_structs.hpp"

//-------------------------- Pull End move -----------------------------



class PullEndCubic : public Move<CubicLattice> {
public:
    PullEndCubic(int seed) : Move<CubicLattice>(seed) {
        moveArgs.backbone_id = 0;
        moveArgs.status = MoveStatus::NotPossible;
        moveArgs.num_monomer_moved = 0;
        moveArgs.dE = 0.0;
        forward = true;
        b_length = 0;
    }

    bool apply(Configuration<CubicLattice>& config) override;
    void generateMoveArgs(Configuration<CubicLattice>& config) override;
    void checkMove(Configuration<CubicLattice>& config) override;
    
    double getdE() { return moveArgs.dE; }
    MoveStatus getStatus() { return moveArgs.status; }
    
    void setBackboneDist(const int b_min, const int b_max) {
        dist_backbone = std::uniform_int_distribution<int>{b_min, b_max};
    }
    

private:
    struct EndMove3D {
		Position A_offset;
		Position B_offset;
	};

	struct EndMovesForDir {
		std::array<EndMove3D, 21> moves;
		size_t count = 21;
	};
    
    std::array<Position,6> dirs {{
		Position{1,0,0}, Position{-1,0,0},
		Position{0,1,0}, Position{0,-1,0},
		Position{0,0,1}, Position{0,0,-1}
	}};


	std::pair<Position, Position> makePerps(const Position& d) {
		if (d[0] != 0) return { Position{0,1,0}, Position{0,0,1} };
		if (d[1] != 0) return { Position{1,0,0}, Position{0,0,1} };
		return { Position{1,0,0}, Position{0,1,0} };
	}


	EndMovesForDir makeEndMovesForDir(const Position& d) {
		auto [perp1, perp2] = makePerps(d);
		Position perp3 = {-perp1[0], -perp1[1], -perp1[2]};
		Position perp4 = {-perp2[0], -perp2[1], -perp2[2]};

		return EndMovesForDir{{
			EndMove3D{d, d},
			{perp1, d}, {perp2, d}, {perp3, d}, {perp4, d},
			{perp1, perp1}, {perp1, perp2}, {perp1, perp4},
			{perp2, perp1}, {perp2, perp2}, {perp2, perp3},
			{perp3, perp2}, {perp3, perp3}, {perp3, perp4},
			{perp4, perp1}, {perp4, perp3}, {perp4, perp4},
			{d, perp1}, {d, perp2}, {d, perp3}, {d, perp4}
		}};
	}

	std::array<EndMovesForDir,6> EndMoves3D {{
		makeEndMovesForDir(dirs[0]),
		makeEndMovesForDir(dirs[1]),
		makeEndMovesForDir(dirs[2]),
		makeEndMovesForDir(dirs[3]),
		makeEndMovesForDir(dirs[4]),
		makeEndMovesForDir(dirs[5])
	}};
	
	PullMoveArgs<CubicLattice> moveArgs;
    bool forward;
    int b_length;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_end{0, 1}; // choose 0 or N-1
    std::uniform_int_distribution<int> dist_move{0, 20};
};

#include "pull_end_cubic.tpp"

#endif // PULL_END_CUBIC_H
