#ifndef TRANSLATION_CUBIC_H
#define TRANSLATION_CUBIC_H

//---------------------------- Translation Cubic ---------------------------
class TranslationCubic : public Move<CubicLattice> {
public :
	
	TranslationCubic(int seed) : Move<CubicLattice>(seed){
		moveArgs.backbone_id = 0;
		moveArgs.status= MoveStatus::NotPossible;
		moveArgs.dE = 0.0;
	}
	
	bool apply(Configuration<CubicLattice>& config) override;
	
	void generateMoveArgs(Configuration<CubicLattice>& config) override;
	
	void checkMove(Configuration<CubicLattice>& config) override;
	
	double getdE() { return moveArgs.dE; }
	
	MoveStatus getStatus() { return moveArgs.status; }
	
	void setBackboneDist(const int b_min, const int b_max) { dist_backbone = std::uniform_int_distribution<int>{b_min, b_max}; }
	
	void setTransRange(const int min, const int max) { dist_position = std::uniform_int_distribution<int>{min, max}; }
		
private:
    TranslationMoveArgs<CubicLattice> moveArgs;
    std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_position{-5, 5};
};

#include "translation_cubic.tpp"

#endif // TRANSLATION_CUBIC_H
