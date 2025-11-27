#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------- Rotation Square --------------------------**/

void RotationSquare::checkMove(Configuration<SquareLattice>& config) {
    auto backbone = config.getBackbone_ptr(moveArgs.backbone_id);

    double dE = 0.0;

    // --- Loop over all monomers in the backbone ---
    for (int i = 0; i < moveArgs.b_length; i++) {
        const auto& old_pos = moveArgs.old_positions[i];
        const auto& new_pos = moveArgs.new_positions[i];

        // Check for occupation at the new site (ignoring self)
        if (config.isSiteOccupiedIgnoringBackbonei(config.applyBoundary(new_pos),
                                                   moveArgs.backbone_id)) {
            moveArgs.status = MoveStatus::NotPossible;
            //~ std::cout << "test" << std::endl;
            return;
        }

        int type = backbone->getType(i);

        // Reset per-monomer contact counts
        moveArgs.old_counts[i] = 0;
        moveArgs.new_counts[i] = 0;

        // Compute local energy + contacts for old position
        dE -= config.computeLocalEnergyAndContactsIgnoringBackbonei(
            old_pos, type, moveArgs.old_contacts[i], moveArgs.old_counts[i], moveArgs.backbone_id);

        // Compute local energy + contacts for new position
        dE += config.computeLocalEnergyAndContactsIgnoringBackbonei(
            new_pos, type, moveArgs.new_contacts[i], moveArgs.new_counts[i], moveArgs.backbone_id);
    }

    moveArgs.dE = dE;
    moveArgs.status = MoveStatus::Possible;
}


bool RotationSquare::apply(Configuration<SquareLattice>& config) {
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

    config.moveBackboneMinimal(bb_id,
                               moveArgs.old_positions,
                               moveArgs.new_positions,
                               moveArgs.b_length,
                               config.getEnergy() + moveArgs.dE);

    return true;
}


void RotationSquare::generateMoveArgs(Configuration<SquareLattice>& config) {
    const size_t bb = dist_backbone(this->rng);
    moveArgs.backbone_id = bb;

    auto backbone = config.getBackbone_ptr(bb);
    moveArgs.b_length = backbone->getLength();

    // Pick rotation angle (multiples of 90 degrees for square lattice)
    moveArgs.rotation_angle = 90 * dist_angle(this->rng);

    // Pick a pivot monomer index inside the backbone
    moveArgs.pivot_index = dist_mono(this->rng);
    //~ std::cout << "pivot index : " << moveArgs.pivot_index << std::endl;

    // The pivot *position* is just that monomer’s position
    moveArgs.rotation_center = backbone->getPosition(moveArgs.pivot_index);
    //~ std::cout << "Rotation center : " << moveArgs.rotation_center[0] << " " << moveArgs.rotation_center[1] << std::endl;

    // Build old/new positions
    for (int i = 0; i < moveArgs.b_length; ++i) {
        const auto& curr = backbone->getPosition(i);
        moveArgs.old_positions[i] = curr;

        // translate relative to pivot
        int x = curr[0] - moveArgs.rotation_center[0];
        int y = curr[1] - moveArgs.rotation_center[1];

        // rotate around pivot
        Position r;
        switch (moveArgs.rotation_angle % 360) {
            case 90:  r = { -y,  x }; break;
            case 180: r = { -x, -y }; break;
            case 270: r = {  y, -x }; break;
            default:  r = {  x,  y }; break; // 0° fallback
        }

        // translate back to lattice coordinates
        r[0] += moveArgs.rotation_center[0];
        r[1] += moveArgs.rotation_center[1];

        moveArgs.new_positions[i] = r;
    }
    
    //~ for (int i = 0; i<moveArgs.b_length; i++) {
		//~ std::cout << "Old position : " << moveArgs.old_positions[i][0] << " " << moveArgs.old_positions[i][1] << std::endl;
		//~ std::cout << "New position : " << moveArgs.new_positions[i][0] << " " << moveArgs.new_positions[i][1] << std::endl;
	//~ }
}
