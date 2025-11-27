#ifndef TAIL_SQUARE_H
#define TAIL_SQUARE_H

//---------------------------- Tail Square -----------------------------
class TailSquare : public Move<SquareLattice> {
public :
	
	TailSquare(int seed) : Move<SquareLattice>(seed){
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
		
private:
    OneMonoMoveArgs<SquareLattice> moveArgs;
    std::uniform_int_distribution<int> dist_dim{0, this->Dim-1};
    std::uniform_int_distribution<int> dist_backbone{1, 1};
    std::uniform_int_distribution<int> dist_monomer{0, 1};
};

#include "tail_square.tpp"

#endif // TAIL_SQUARE_H
