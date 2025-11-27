#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------- Crankshaft Square ---------------------------**/

void CrankshaftSquare::checkMove(Configuration<SquareLattice>& config) {//, const std::any moveArgs) {

	if (!config.areAdjacent(moveArgs.pos_mono_0, moveArgs.pos_mono_3) || config.isSiteOccupied(config.applyBoundary(moveArgs.new_position1)) || config.isSiteOccupied(config.applyBoundary(moveArgs.new_position2))) {
		moveArgs.status = MoveStatus::NotPossible;
		return;
	}

    auto global_id1 = config.computeGlobalID(moveArgs.backbone_id1, moveArgs.monomer_id1);
	auto global_id2 = config.computeGlobalID(moveArgs.backbone_id2, moveArgs.monomer_id2);

    // old energy and old contacts
    moveArgs.old_count1 = 0;
    moveArgs.old_count2 = 0;
    double E_old = config.computeLocalEnergyAndContacts(moveArgs.old_position1, global_id1, moveArgs.type1, moveArgs.old_contacts1, moveArgs.old_count1)
				 + config.computeLocalEnergyAndContacts(moveArgs.old_position2, global_id2, moveArgs.type2, moveArgs.old_contacts2, moveArgs.old_count2);
	
    // new energy and new contacts
    moveArgs.new_count1 = 0;
    moveArgs.new_count2 = 0;
    double E_new = config.computeLocalEnergyAndContacts(moveArgs.new_position1, global_id1, moveArgs.type1, moveArgs.new_contacts1, moveArgs.new_count1)
				 + config.computeLocalEnergyAndContacts(moveArgs.new_position2, global_id2, moveArgs.type2, moveArgs.new_contacts2, moveArgs.new_count2);

    moveArgs.dE = E_new - E_old;
    
    moveArgs.status = MoveStatus::Possible;
}

bool CrankshaftSquare::apply(Configuration<SquareLattice>& config) {//, const std::any moveArgs) {
	config.moveMonomerMinimal(moveArgs.backbone_id1, moveArgs.monomer_id1, moveArgs.old_position1, moveArgs.new_position1, config.getEnergy() + moveArgs.dE);
	config.moveMonomerOnly(moveArgs.backbone_id2, moveArgs.monomer_id2, moveArgs.old_position2, moveArgs.new_position2);
	
	// Update contact matrix using packed upper triangle
    const int self_bb = moveArgs.backbone_id1;

    // decrement old contacts
    for (int i = 0; i < moveArgs.old_count1; ++i) {
        int nb = moveArgs.old_contacts1[i];
        config.contact_matrix.dec(self_bb, nb);
    }
    for (int i = 0; i < moveArgs.old_count2; ++i) {
        int nb = moveArgs.old_contacts2[i];
        config.contact_matrix.dec(self_bb, nb);
    }

    // increment new contacts
    for (int i = 0; i < moveArgs.new_count1; ++i) {
        int nb = moveArgs.new_contacts1[i];
        config.contact_matrix.inc(self_bb, nb);
    }
    for (int i = 0; i < moveArgs.new_count2; ++i) {
        int nb = moveArgs.new_contacts2[i];
        config.contact_matrix.inc(self_bb, nb);
    }
	
	return true;
}


void CrankshaftSquare::generateMoveArgs(Configuration<SquareLattice>& config) {
	
	const size_t backbone_id = dist_backbone(this->rng);
	const size_t monomer_id = dist_monomer(this->rng);
	
	auto b_ptr = config.getBackbone_ptr(backbone_id);
	
	const Position& curr = b_ptr->getPosition(monomer_id);
	const Position& next = b_ptr->getPosition(monomer_id+1);
	const Position& nnext = b_ptr->getPosition(monomer_id+2);
	const Position& nnnext = b_ptr->getPosition(monomer_id+3);
	
	moveArgs.old_position1 = next;
	moveArgs.old_position2 = nnext;
	
	moveArgs.new_position1[0] = 2 * curr[0] - next[0];
	moveArgs.new_position1[1] = 2 * curr[1] - next[1];
	
	moveArgs.new_position2[0] = 2 * nnnext[0] - nnext[0];
	moveArgs.new_position2[1] = 2 * nnnext[1] - nnext[1];
	
	
	moveArgs.backbone_id1 = backbone_id;
	moveArgs.backbone_id2 = backbone_id;
	moveArgs.monomer_id1 = monomer_id+1;
	moveArgs.monomer_id2 = monomer_id+2;
	moveArgs.type1 = b_ptr->getType(monomer_id+1);
	moveArgs.type2 = b_ptr->getType(monomer_id+2);
	moveArgs.pos_mono_0 = curr;
	moveArgs.pos_mono_3 = nnnext;
}
