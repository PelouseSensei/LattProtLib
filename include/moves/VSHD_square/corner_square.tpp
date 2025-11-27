#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------- Corner Square --------------------------**/

void CornerSquare::checkMove(Configuration<SquareLattice>& config) {
    if (config.isSiteOccupied(config.applyBoundary(moveArgs.new_position))) {
        moveArgs.status = MoveStatus::NotPossible;
        return;
    }

    auto global_id = config.computeGlobalID(moveArgs.backbone_id, moveArgs.monomer_id);

    // old energy and old contacts
    moveArgs.old_count = 0;
    double E_old = config.computeLocalEnergyAndContacts(
        moveArgs.old_position,
        global_id,
        moveArgs.type,
        moveArgs.old_contacts,
        moveArgs.old_count
    );

    // new energy and new contacts
    moveArgs.new_count = 0;
    double E_new = config.computeLocalEnergyAndContacts(
        moveArgs.new_position,
        global_id,
        moveArgs.type,
        moveArgs.new_contacts,
        moveArgs.new_count
    );

    moveArgs.dE = E_new - E_old;
    moveArgs.status = MoveStatus::Possible;
}

bool CornerSquare::apply(Configuration<SquareLattice>& config) {
    // Move monomer in lattice
    config.moveMonomerMinimal(moveArgs.backbone_id,
                              moveArgs.monomer_id,
                              moveArgs.old_position,
                              moveArgs.new_position,
                              config.getEnergy() + moveArgs.dE);

    // Update contact matrix using packed upper triangle
    const int self_bb = moveArgs.backbone_id;

    // decrement old contacts
    for (int i = 0; i < moveArgs.old_count; ++i) {
        int nb = moveArgs.old_contacts[i];
        config.contact_matrix.dec(self_bb, nb);  // automatically handles i>j
    }

    // increment new contacts
    for (int i = 0; i < moveArgs.new_count; ++i) {
        int nb = moveArgs.new_contacts[i];
        config.contact_matrix.inc(self_bb, nb);
    }

    return true;
}


void CornerSquare::generateMoveArgs(Configuration<SquareLattice>& config) {
	
	const size_t backbone_id = dist_backbone(this->rng);
	const size_t monomer_id = dist_monomer(this->rng);
	
	auto b_ptr = config.getBackbone_ptr(backbone_id);
	
	const Position& prev = b_ptr->getPosition(monomer_id - 1);
	const Position& next = b_ptr->getPosition(monomer_id + 1);
	const Position& curr = b_ptr->getPosition(monomer_id);

	
	moveArgs.old_position = curr;
	
	moveArgs.new_position[0] = prev[0] + next[0] - curr[0];
	moveArgs.new_position[1] = prev[1] + next[1] - curr[1];

	
	moveArgs.backbone_id = backbone_id;
	moveArgs.monomer_id = monomer_id;
	moveArgs.type = b_ptr->getType(monomer_id);
}
