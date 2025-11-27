#ifndef TWO_POINT_SWITCH_H
#define TWO_POINT_SWITCH_H

//---------------------------- Two Point Type Switch Move ---------------------------
template<typename LatticeType>
class TwoPointSwitch : public Move<LatticeType> {
public:
	
	TwoPointSwitch(int seed) : Move<LatticeType>(seed){
		moveArgs.backbone_id1 = 0;
		moveArgs.backbone_id2 = 0;
		moveArgs.dE = 0;
		moveArgs.old_type1 = 0;
		moveArgs.old_type2 = 0;
	}
	
	bool apply(Configuration<LatticeType>& config) override;
	
	void generateMoveArgs(Configuration<LatticeType>& config) override;
	
	void checkMove(Configuration<LatticeType>& config) override;
	
	double getdE() { return moveArgs.dE; }
	
	MoveStatus getStatus() { return moveArgs.status; }
	
	void setBackboneDist(const int b_min, const int b_max) { dist_backbone = std::uniform_int_distribution<int>{b_min, b_max}; }
	
	void setMonomerDist (const int min, const int max) { dist_monomer = std::uniform_int_distribution<int>{min, max}; }
	

private:
	TwoTypeMoveArgs<LatticeType> moveArgs;
	std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_monomer{0, 0};
    

};

#include "two_point_switch.tpp"

#endif // TWO_POINT_SWITCH_H
