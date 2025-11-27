#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**-------------------------- Pull End move ---------------------------**/
void PullEndCubic::checkMove(Configuration<CubicLattice>& config) {
    if (moveArgs.status == MoveStatus::NotPossible) return;

    const int bb = moveArgs.backbone_id;
    auto backbone = config.getBackbone_ptr(bb);
    double dE = 0.0;

    // Precompute types
    for (int k = 0; k < moveArgs.num_monomer_moved; ++k) {
        int idx = forward ? moveArgs.start_index + k : moveArgs.start_index - k;
        moveArgs.types[k] = backbone->getType(idx);
    }

    // --- Remove old energy ---
    dE -= config.computeBlockEnergyAndContacts(
        bb,
        moveArgs.start_index,
        moveArgs.num_monomer_moved,
        moveArgs.old_positions,
        moveArgs.types,
        moveArgs.old_contacts,
        moveArgs.old_counts,
        forward
    );

    // --- Temporarily move ---
    if (forward) {
        config.moveBackboneMinimalPullForward(bb, moveArgs.old_positions,
            moveArgs.new_positions, moveArgs.start_index, moveArgs.num_monomer_moved, config.getEnergy());
    } else {
        config.moveBackboneMinimalPullBackward(bb, moveArgs.old_positions,
            moveArgs.new_positions, moveArgs.start_index, moveArgs.num_monomer_moved, config.getEnergy());
    }

    // --- Add new energy ---
    dE += config.computeBlockEnergyAndContacts(
        bb,
        moveArgs.start_index,
        moveArgs.num_monomer_moved,
        moveArgs.new_positions,
        moveArgs.types,
        moveArgs.new_contacts,
        moveArgs.new_counts,
        forward
    );

    // --- Restore lattice ---
    if (forward) {
        config.moveBackboneMinimalPullForward(bb, moveArgs.new_positions,
            moveArgs.old_positions, moveArgs.start_index, moveArgs.num_monomer_moved, config.getEnergy());
    } else {
        config.moveBackboneMinimalPullBackward(bb, moveArgs.new_positions,
            moveArgs.old_positions, moveArgs.start_index, moveArgs.num_monomer_moved, config.getEnergy());
    }

    moveArgs.dE = dE;
    moveArgs.status = MoveStatus::Possible;
}


bool PullEndCubic::apply(Configuration<CubicLattice>& config) {
    if (moveArgs.status != MoveStatus::Possible) return false;

    const int bb = moveArgs.backbone_id;

    // Remove old contacts
    for (int k = 0; k < moveArgs.num_monomer_moved; ++k)
        for (int c = 0; c < moveArgs.old_counts[k]; ++c)
            config.contact_matrix.dec(bb, moveArgs.old_contacts[k][c]);

    // Add new contacts
    for (int k = 0; k < moveArgs.num_monomer_moved; ++k)
        for (int c = 0; c < moveArgs.new_counts[k]; ++c)
            config.contact_matrix.inc(bb, moveArgs.new_contacts[k][c]);

    // Apply move
    if (forward) {
        config.moveBackboneMinimalPullForward(bb, moveArgs.old_positions,
            moveArgs.new_positions, moveArgs.start_index, moveArgs.num_monomer_moved, config.getEnergy() + moveArgs.dE);
    } else {
        config.moveBackboneMinimalPullBackward(bb, moveArgs.old_positions,
            moveArgs.new_positions, moveArgs.start_index, moveArgs.num_monomer_moved, config.getEnergy() + moveArgs.dE);
    }

    return true;
}



void PullEndCubic::generateMoveArgs(Configuration<CubicLattice>& config) {
    const size_t bb = dist_backbone(this->rng);
    moveArgs.backbone_id = bb;
    auto backbone = config.getBackbone_ptr(bb);
    const int b_length = backbone->getLength();

    // need at least two beads
    if (b_length < 2) {
        moveArgs.status = MoveStatus::NotPossible;
        return;
    }

    // choose end (forward/backward)
    forward = (dist_end(this->rng) == 0);
    moveArgs.start_index = forward ? 0 : b_length - 1;
    const int neighbor_idx = forward ? 1 : b_length - 2;

    const Position end_pos = backbone->getPosition(moveArgs.start_index);
    const Position neighbor_pos = backbone->getPosition(neighbor_idx);

    // bond vector d = end - neighbor (unit vector)
    Position d = { end_pos[0]-neighbor_pos[0], end_pos[1]-neighbor_pos[1], end_pos[2]-neighbor_pos[2] };

	// find index 0..5
	int d_idx = 0;
	for (int i=0; i<6; ++i)
		if (dirs[i][0]==d[0] && dirs[i][1]==d[1] && dirs[i][2]==d[2]) { d_idx=i; break; }

	const auto &tbl = EndMoves3D[d_idx];
	//~ std::uniform_int_distribution<int> dist_move(0, (int)tbl.count-1);
	const EndMove3D &tpl = tbl.moves[dist_move(this->rng)];

	Position A_raw, B_raw;
	for (size_t i=0;i<3;++i){
		A_raw[i] = end_pos[i] + tpl.A_offset[i];
		B_raw[i] = A_raw[i] + tpl.B_offset[i];
	}


    // apply BC for occupancy checks
    const Position A_wrapped = config.applyBoundary(A_raw);
    const Position B_wrapped = config.applyBoundary(B_raw);

    // Feasibility checks (only chosen move)
    if (config.isSiteOccupied(A_wrapped) || config.isSiteOccupied(B_wrapped)) {
        moveArgs.status = MoveStatus::NotPossible;
        return;
    }

    // initialize first two moved monomers (end and neighbor)
    moveArgs.num_monomer_moved = 2;
    moveArgs.old_positions[0] = end_pos;        // old end
    moveArgs.new_positions[0] = B_raw;          // new end (raw)
    moveArgs.old_positions[1] = neighbor_pos;   // old neighbor
    moveArgs.new_positions[1] = A_raw;          // new neighbor (raw)

    // propagate along backbone (same logic as 2D, but positions are 3D)
    int offset = 2;
    int i = forward ? 2 : b_length - 3;
    while (forward ? (i < b_length) : (i >= 0)) {
        const Position curr_pos = backbone->getPosition(i);

        // stop if current is adjacent to the last new position
        if (config.areAdjacent(curr_pos, moveArgs.new_positions[offset - 1])) break;

        moveArgs.old_positions[offset] = curr_pos;
        moveArgs.new_positions[offset] = forward ? backbone->getPosition(i - 2)
                                                 : backbone->getPosition(i + 2);
        moveArgs.num_monomer_moved++;
        if (forward) ++i; else --i;
        ++offset;
    }

    moveArgs.status = MoveStatus::Possible;
}



