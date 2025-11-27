#ifndef ONE_POINT_MUTATION_H
#define ONE_POINT_MUTATION_H

//---------------------------- One point Mutation Move ---------------------------
template<typename LatticeType>
class OnePointMutation : public Move<LatticeType> {
public:
	
	OnePointMutation(int seed) : Move<LatticeType>(seed){
		moveArgs.backbone_id = 0;
		moveArgs.dE = 0;
		moveArgs.new_type = 0;
		moveArgs.old_type = 0;
	}
	
	bool apply(Configuration<LatticeType>& config) override;
	
	void generateMoveArgs(Configuration<LatticeType>& config) override;
	
	void checkMove(Configuration<LatticeType>& config) override;
	
	void revert(Configuration<LatticeType>& config);
	
	double getdE() { return moveArgs.dE; }
	
	MoveStatus getStatus() { return moveArgs.status; }
	
	void setBackboneDist(const int b_min, const int b_max) { dist_backbone = std::uniform_int_distribution<int>{b_min, b_max}; }
	
	void setMonomerDist (const int min, const int max) { dist_monomer = std::uniform_int_distribution<int>{min, max}; }
	
	void setTypeDist (const int min, const int max) { dist_type = std::uniform_int_distribution<int>{min, max}; }
	
	int getMoveOldType() { return moveArgs.old_type; }
	
	int getMoveNewType() { return moveArgs.new_type; }
	
	int getMoveMono() { return moveArgs.monomer_id; }
	
	int getMoveBackbone() { return moveArgs.backbone_id; }


private:
	OneTypeMoveArgs<LatticeType> moveArgs;
	std::uniform_int_distribution<int> dist_backbone{0, 0};
    std::uniform_int_distribution<int> dist_monomer{0, 0};
    std::uniform_int_distribution<int> dist_type{0, 0};
    

};

#include "one_point_mutation.tpp"

#endif // ONE_POINT_MUTATION_H
