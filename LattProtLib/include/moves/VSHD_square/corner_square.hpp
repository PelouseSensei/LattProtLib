#ifndef CORNER_SQUARE_H
#define CORNER_SQUARE_H

//--------------------------- Corner Square ----------------------------
class CornerSquare : public Move<SquareLattice> {
public :
	
	CornerSquare(int seed) : Move<SquareLattice>(seed){
		moveArgs.backbone_id = 0;
		moveArgs.monomer_id = 0;
		moveArgs.status= MoveStatus::NotPossible;
		moveArgs.type = 0;
		moveArgs.dE = 0.0;
	}
	
	bool apply(Configuration<SquareLattice>& config) override;
	
	void generateMoveArgs(Configuration<SquareLattice>& config) override;
	
	void checkMove(Configuration<SquareLattice>& config) override;
	
	double getdE() { return moveArgs.dE; }
	
	MoveStatus getStatus() { return moveArgs.status; }
	
	void setBackboneDist(const int b_min, const int b_max) { dist_backbone = std::uniform_int_distribution<int>{b_min, b_max}; }
	
	void setMonomerDist(const int m_min, const int m_max) { dist_monomer = std::uniform_int_distribution<int>{m_min, m_max}; }
		
private:
    OneMonoMoveArgs<SquareLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{1, 1};
    std::uniform_int_distribution<int> dist_monomer{1, 3};
};

#include "corner_square.tpp"

#endif // CORNER_SQUARE_H
