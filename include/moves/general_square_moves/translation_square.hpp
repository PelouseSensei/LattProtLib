#ifndef TRANSLATION_SQUARE_H
#define TRANSLATION_SQUARE_H

//---------------------------- Translation Square ---------------------------
class TranslationSquare : public Move<SquareLattice> {
public :
	
	TranslationSquare(int seed) : Move<SquareLattice>(seed){
		moveArgs.backbone_id = 0;
		moveArgs.status= MoveStatus::NotPossible;
		moveArgs.dE = 0.0;
	}
	
	bool apply(Configuration<SquareLattice>& config) override;
	
	void generateMoveArgs(Configuration<SquareLattice>& config) override;
	
	void checkMove(Configuration<SquareLattice>& config) override;
	
	double getdE() { return moveArgs.dE; }
	
	MoveStatus getStatus() { return moveArgs.status; }
	
	void setBackboneDist(const int b_min, const int b_max) { dist_backbone = std::uniform_int_distribution<int>{b_min, b_max}; }
	
	void setTransRange(const int min, const int max) { dist_position = std::uniform_int_distribution<int>{min, max}; }
		
private:
    TranslationMoveArgs<SquareLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_position{-5, 5};
};

#include "translation_square.tpp"

#endif // TRANSLATION_SQUARE_H
