#ifndef CRANKSHAFT_SQUARE_H
#define CRANKSHAFT_SQUARE_H

//---------------------------- Crankshaft Square -----------------------------
class CrankshaftSquare : public Move<SquareLattice> {
public :
	
	CrankshaftSquare(const int seed) : Move<SquareLattice>(seed) {
		moveArgs.backbone_id1 = 0;
		moveArgs.backbone_id2 = 0;
		moveArgs.monomer_id1 = 0;
		moveArgs.monomer_id2 = 0;
		moveArgs.status= MoveStatus::NotPossible;
		moveArgs.type1 = 0;
		moveArgs.type2 = 0;
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
    TwoMonoMoveArgs<SquareLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{1, 1};
    std::uniform_int_distribution<int> dist_monomer{0,0};
};

#include "crankshaft_square.tpp"

#endif // CRANKSHAFT_SQUARE_H
