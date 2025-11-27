#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**-------------------------- Type switch ---------------------------**/
template<typename LatticeType>
void TwoPointSwitch<LatticeType>::checkMove(Configuration<LatticeType>& config) {
	
    auto global_id1 = config.computeGlobalID(moveArgs.backbone_id1, moveArgs.monomer_id1);
    auto global_id2 = config.computeGlobalID(moveArgs.backbone_id2, moveArgs.monomer_id2);
	
	double E_old = config.computeLocalEnergy(moveArgs.position1, global_id1, moveArgs.old_type1)
				 + config.computeLocalEnergy(moveArgs.position2, global_id2, moveArgs.old_type2);
	
	config.getBackbone_ptr(moveArgs.backbone_id1)->setType(moveArgs.monomer_id1, moveArgs.old_type2);
	config.getBackbone_ptr(moveArgs.backbone_id2)->setType(moveArgs.monomer_id2, moveArgs.old_type1);
	
	double E_new = config.computeLocalEnergy(moveArgs.position1, global_id1, moveArgs.old_type2)
				 + config.computeLocalEnergy(moveArgs.position2, global_id2, moveArgs.old_type1);
	
	config.getBackbone_ptr(moveArgs.backbone_id1)->setType(moveArgs.monomer_id1, moveArgs.old_type1);
	config.getBackbone_ptr(moveArgs.backbone_id2)->setType(moveArgs.monomer_id2, moveArgs.old_type2);
	
    moveArgs.dE = E_new - E_old;
    //~ std::cout << "dE: " << moveArgs.dE << std::endl;
    moveArgs.status = MoveStatus::Possible;
}

template<typename LatticeType>
bool TwoPointSwitch<LatticeType>::apply(Configuration<LatticeType>& config) {
	//~ std::cout << "Energy before: " << config.getEnergy() << std::endl;
	config.mutateMonomer(moveArgs.backbone_id1, moveArgs.monomer_id1, moveArgs.old_type2, config.getEnergy() + moveArgs.dE);
	//~ std::cout << "Energy in between: " << config.getEnergy() << std::endl;
	config.mutateMonomer(moveArgs.backbone_id2, moveArgs.monomer_id2, moveArgs.old_type1, config.getEnergy());
	//~ std::cout << "Energy after: " << config.getEnergy() << std::endl;
	return true;
}

template<typename LatticeType>
void TwoPointSwitch<LatticeType>::generateMoveArgs(Configuration<LatticeType>& config) {
	
	const size_t backbone_id1 = dist_backbone(this->rng);
	const size_t backbone_id2 = dist_backbone(this->rng);
	moveArgs.backbone_id1 = backbone_id1;
	moveArgs.backbone_id2 = backbone_id2;
	
	const int monomer_id1 = dist_monomer(this->rng);
	const int monomer_id2 = dist_monomer(this->rng);
	moveArgs.monomer_id1 = monomer_id1;
	moveArgs.monomer_id2 = monomer_id2;
	
	auto backbone1 = config.getBackbone_ptr(backbone_id1);
	auto backbone2 = config.getBackbone_ptr(backbone_id2);
	moveArgs.old_type1 = backbone1->getType(monomer_id1);
	moveArgs.old_type2 = backbone2->getType(monomer_id2);
	moveArgs.position1 = backbone1->getPosition(monomer_id1);
	moveArgs.position2 = backbone2->getPosition(monomer_id2);
	
	//~ std::cout << "move : \n";
	//~ std::cout << "Monomers: " << moveArgs.monomer_id1 << " " << moveArgs.monomer_id2 << std::endl;
	//~ std::cout << "Types: " << moveArgs.old_type1 << " " << moveArgs.old_type2 << std::endl;
	
}
