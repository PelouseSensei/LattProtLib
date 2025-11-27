#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------ One Point Mutation -------------------------**/
template<typename LatticeType>
void OnePointMutation<LatticeType>::checkMove(Configuration<LatticeType>& config) {
	
    auto global_id = config.computeGlobalID(moveArgs.backbone_id, moveArgs.monomer_id);
	
	double E_old = config.computeLocalEnergy(moveArgs.position, global_id, moveArgs.old_type);
	double E_new = config.computeLocalEnergy(moveArgs.position, global_id, moveArgs.new_type);

    moveArgs.dE = E_new - E_old;
    moveArgs.status = MoveStatus::Possible;
}

template<typename LatticeType>
bool OnePointMutation<LatticeType>::apply(Configuration<LatticeType>& config) {
	config.mutateMonomer(moveArgs.backbone_id, moveArgs.monomer_id, moveArgs.new_type, config.getEnergy() + moveArgs.dE);
	return true;
}

template<typename LatticeType>
void OnePointMutation<LatticeType>::generateMoveArgs(Configuration<LatticeType>& config) {
	
	const size_t backbone_id = dist_backbone(this->rng);
	moveArgs.backbone_id = backbone_id;
	
	const int monomer_id = dist_monomer(this->rng);
	moveArgs.monomer_id = monomer_id;
	
	auto backbone = config.getBackbone_ptr(backbone_id);
	moveArgs.old_type = backbone->getType(monomer_id);
	moveArgs.position = backbone->getPosition(monomer_id);
	
	moveArgs.new_type = dist_type(this->rng);
	
}

template<typename LatticeType>
void OnePointMutation<LatticeType>::revert(Configuration<LatticeType>& config) {
	
	auto global_id = config.computeGlobalID(moveArgs.backbone_id, moveArgs.monomer_id);
	auto backbone = config.getBackbone_ptr(moveArgs.backbone_id);
	auto position = backbone->getPosition(moveArgs.monomer_id);
	
	double E_old = config.computeLocalEnergy(position, global_id, moveArgs.new_type);
	double E_new = config.computeLocalEnergy(position, global_id, moveArgs.old_type);

    double dE = E_new - E_old;
    //~ moveArgs.status = MoveStatus::Possible;
	
	config.mutateMonomer(moveArgs.backbone_id, moveArgs.monomer_id, moveArgs.old_type, config.getEnergy() + dE);
}
