#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------- Tail Square ---------------------------**/

void TailSquare::checkMove(Configuration<SquareLattice>& config) {
    
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
        //~ moveArgs.backbone_id,
        moveArgs.old_contacts,
        moveArgs.old_count
    );

    // new energy and new contacts
    moveArgs.new_count = 0;
    double E_new = config.computeLocalEnergyAndContacts(
        moveArgs.new_position,
        global_id,
        moveArgs.type,
        //~ moveArgs.backbone_id,
        moveArgs.new_contacts,
        moveArgs.new_count
    );

    moveArgs.dE = E_new - E_old;
    moveArgs.status = MoveStatus::Possible;   
}


bool TailSquare::apply(Configuration<SquareLattice>& config) {
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


void TailSquare::generateMoveArgs(Configuration<SquareLattice>& config) {
	
	const size_t backbone_id = dist_backbone(this->rng);
	const size_t monomer_id = dist_monomer(this->rng);
	
	auto b_ptr = config.getBackbone_ptr(backbone_id);
	
	if ( monomer_id == 0) {
		moveArgs.old_position = b_ptr->getPosition(monomer_id);
		moveArgs.new_position = b_ptr->getPosition(monomer_id + 1);
		
		int dim = dist_dim(this->rng);
		
		int dir = dist_monomer(this->rng);
		
		moveArgs.new_position[dim] = moveArgs.new_position[dim] + 2 * ( dir - 0.5) ;
		
		moveArgs.backbone_id = backbone_id;
		moveArgs.monomer_id = 0;
		moveArgs.type = b_ptr->getType(monomer_id);
		
	} else {
		auto b_length = b_ptr->getLength();
		moveArgs.old_position = b_ptr->getPosition(b_length-1);
		moveArgs.new_position = b_ptr->getPosition(b_length -2);
		
		int dim = dist_dim(this->rng);
		
		int dir = dist_monomer(this->rng);
		
		moveArgs.new_position[dim] = moveArgs.new_position[dim] + 2 * ( dir - 0.5) ; // Maybe having a random dist in {0, 1} is faster ?
		
		moveArgs.backbone_id = backbone_id;
		moveArgs.monomer_id = b_length-1;
		moveArgs.type = b_ptr->getType(b_length-1);
	}
	
	
}
