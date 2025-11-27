#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------- Translation Square --------------------------**/

void TranslationSquare::checkMove(Configuration<SquareLattice>& config) {
    auto backbone = config.getBackbone_ptr(moveArgs.backbone_id);
    //~ auto global_id = config.computeGlobalID(moveArgs.backbone_id, 0);

    double dE = 0.0;

    // --- Loop over all monomers in the backbone ---
    for (int i = 0; i < moveArgs.b_length; i++) {
        const auto& old_pos = moveArgs.old_positions[i];
        const auto& new_pos = moveArgs.new_positions[i];

        // Check for occupation at the new site (ignoring self)
        if (config.isSiteOccupiedIgnoringBackbonei(config.applyBoundary(new_pos),
                                                   moveArgs.backbone_id)) {
            moveArgs.status = MoveStatus::NotPossible;
            return;
        }

        //~ int g_id = global_id + i;
        int type = backbone->getType(i);

        // Reset per-monomer contact counts
        moveArgs.old_counts[i] = 0;
        moveArgs.new_counts[i] = 0;

        // Compute local energy + contacts for old position
        dE -= config.computeLocalEnergyAndContactsIgnoringBackbonei(old_pos, type, moveArgs.old_contacts[i], moveArgs.old_counts[i], moveArgs.backbone_id);

        // Compute local energy + contacts for new position
        dE += config.computeLocalEnergyAndContactsIgnoringBackbonei(new_pos, type, moveArgs.new_contacts[i], moveArgs.new_counts[i], moveArgs.backbone_id);
    }

    moveArgs.dE = dE;
    moveArgs.status = MoveStatus::Possible;
}


bool TranslationSquare::apply(Configuration<SquareLattice>& config) {
    int bb_id = moveArgs.backbone_id;

    // --- 1. Remove old contacts ---
    for (int i = 0; i < moveArgs.b_length; i++) {
        for (int c = 0; c < moveArgs.old_counts[i]; c++) {
            int neighbor_bb = moveArgs.old_contacts[i][c];
            config.contact_matrix.dec(bb_id, neighbor_bb);
        }
    }

    // --- 2. Add new contacts ---
    for (int i = 0; i < moveArgs.b_length; i++) {
        for (int c = 0; c < moveArgs.new_counts[i]; c++) {
            int neighbor_bb = moveArgs.new_contacts[i][c];
            config.contact_matrix.inc(bb_id, neighbor_bb);
        }
    }
    config.moveBackboneMinimal(bb_id, moveArgs.old_positions, moveArgs.new_positions, moveArgs.b_length, config.getEnergy() + moveArgs.dE);

    return true;
}




void TranslationSquare::generateMoveArgs(Configuration<SquareLattice>& config) {
    const size_t bb = dist_backbone(this->rng);
    moveArgs.backbone_id = bb;

    auto backbone = config.getBackbone_ptr(bb);
    moveArgs.b_length = backbone->getLength();

	moveArgs.translation_vec[0] = dist_position(this->rng);
	moveArgs.translation_vec[1] = dist_position(this->rng);

    for (int i = 0; i < moveArgs.b_length; ++i) {
        const auto& curr = backbone->getPosition(i);
        moveArgs.old_positions[i] = curr;
        
        Position t = curr;
        t[0] += moveArgs.translation_vec[0];
        t[1] += moveArgs.translation_vec[1];
        moveArgs.new_positions[i] = t;
    }
}
