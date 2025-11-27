#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**-------------------------- Pull move ---------------------------**/
template<bool Forward>
void PullSquare<Forward>::checkMove(Configuration<SquareLattice>& config) {
    if (moveArgs.status == MoveStatus::NotPossible) return;

    const int bb = moveArgs.backbone_id;
    auto backbone = config.getBackbone_ptr(bb);

    double dE = 0.0;

    // --- Precompute types and global IDs for moved monomers ---
    for (int k = 0; k < moveArgs.num_monomer_moved; ++k) {
        int idx = Forward ? moveArgs.start_index + k : moveArgs.start_index - k;
        moveArgs.types[k] = backbone->getType(idx);
        moveArgs.global_ids[k] = config.computeGlobalID(bb, idx);
    }

    // --- 1. Remove old contributions using block energy ---
    dE -= config.computeBlockEnergyAndContacts(
        bb,
        moveArgs.start_index,
        moveArgs.num_monomer_moved,
        moveArgs.old_positions,
        moveArgs.types,
        moveArgs.old_contacts,
        moveArgs.old_counts,
        Forward
    );

    // --- 2. Temporarily move monomers into lattice ---
    if constexpr (Forward) {
        config.moveBackboneMinimalPullForward(
            bb, moveArgs.old_positions, moveArgs.new_positions,
            moveArgs.start_index, moveArgs.num_monomer_moved,
            config.getEnergy()
        );
    } else {
        config.moveBackboneMinimalPullBackward(
            bb, moveArgs.old_positions, moveArgs.new_positions,
            moveArgs.start_index, moveArgs.num_monomer_moved,
            config.getEnergy()
        );
    }

    // --- 3. Add new contributions using block energy ---
    dE += config.computeBlockEnergyAndContacts(
        bb,
        moveArgs.start_index,
        moveArgs.num_monomer_moved,
        moveArgs.new_positions,
        moveArgs.types,
        moveArgs.new_contacts,
        moveArgs.new_counts,
        Forward
    );

    // --- 4. Restore lattice to old state ---
    if constexpr (Forward) {
        config.moveBackboneMinimalPullForward(
            bb, moveArgs.new_positions, moveArgs.old_positions,
            moveArgs.start_index, moveArgs.num_monomer_moved,
            config.getEnergy()
        );
    } else {
        config.moveBackboneMinimalPullBackward(
            bb, moveArgs.new_positions, moveArgs.old_positions,
            moveArgs.start_index, moveArgs.num_monomer_moved,
            config.getEnergy()
        );
    }

    // --- 5. Save results ---
    moveArgs.dE = dE;
    moveArgs.status = MoveStatus::Possible;
}

template<bool Forward>
bool PullSquare<Forward>::apply(Configuration<SquareLattice>& config) {
    if (moveArgs.status != MoveStatus::Possible) return false;

    const int bb = moveArgs.backbone_id;

    // --- 1. Remove old contacts ---
    for (int k = 0; k < moveArgs.num_monomer_moved; ++k)
        for (int c = 0; c < moveArgs.old_counts[k]; ++c)
            config.contact_matrix.dec(bb, moveArgs.old_contacts[k][c]);

    // --- 2. Add new contacts ---
    for (int k = 0; k < moveArgs.num_monomer_moved; ++k)
        for (int c = 0; c < moveArgs.new_counts[k]; ++c)
            config.contact_matrix.inc(bb, moveArgs.new_contacts[k][c]);

    // --- 3. Apply actual backbone move ---
    if constexpr (Forward) {
        config.moveBackboneMinimalPullForward(
            bb, moveArgs.old_positions, moveArgs.new_positions,
            moveArgs.start_index, moveArgs.num_monomer_moved,
            config.getEnergy() + moveArgs.dE
        );
    } else {
        config.moveBackboneMinimalPullBackward(
            bb, moveArgs.old_positions, moveArgs.new_positions,
            moveArgs.start_index, moveArgs.num_monomer_moved,
            config.getEnergy() + moveArgs.dE
        );
    }

    return true;
}

template<bool Forward>
void PullSquare<Forward>::generateMoveArgs(Configuration<SquareLattice>& config) {
    const size_t bb = dist_backbone(this->rng);
    moveArgs.backbone_id = bb;

    auto backbone = config.getBackbone_ptr(bb);
    moveArgs.start_index = Forward ? dist_start(this->rng) : dist_start_backward(this->rng);

    // Neighbor index and positions
    int neighbor_idx = Forward ? moveArgs.start_index - 1 : moveArgs.start_index + 1;
    Position start_pos = backbone->getPosition(moveArgs.start_index);
    Position neighbor_pos = backbone->getPosition(neighbor_idx);

    // Relative displacement
    int dx = neighbor_pos[0] - start_pos[0];
    int dy = neighbor_pos[1] - start_pos[1];


    // --- Always pick one perpendicular site at random ---
	Position perp1 = {neighbor_pos[0] - dy, neighbor_pos[1] + dx};
	Position perp2 = {neighbor_pos[0] + dy, neighbor_pos[1] - dx};
	std::uniform_int_distribution<int> dist_perp(0, 1);
	Position chosen = (dist_perp(this->rng) == 0) ? perp1 : perp2;
	//~ chosen = config.applyBoundary(chosen);

	// --- Check occupancy of chosen site ---
	if (config.isSiteOccupied(config.applyBoundary(chosen))) {
		moveArgs.status = MoveStatus::NotPossible;
		return;
}

    // --- Compute L (propagation target) ---
    Position L = {chosen[0] - dx, chosen[1] - dy};
    //~ Position wrapped_L = config.applyBoundary(L);
    int g_id_perp_L = config.checkSite(config.applyBoundary(L));
    int g_id_for_corner = config.computeGlobalID(bb, Forward ? moveArgs.start_index + 1 : moveArgs.start_index - 1);

	// --- Always enter loop ---
	if (g_id_perp_L == g_id_for_corner) {
		// Corner move
		moveArgs.num_monomer_moved = 1;
		moveArgs.old_positions[0] = start_pos;
		moveArgs.new_positions[0] = chosen;
		moveArgs.status = MoveStatus::Possible;
		return;
	}
	else if (g_id_perp_L != -1) {
		// Occupied by something else → not possible
		moveArgs.status = MoveStatus::NotPossible;
		return;
	}

	// --- Standard propagation ---
	moveArgs.num_monomer_moved = 2;
	moveArgs.old_positions[0] = start_pos;
	moveArgs.new_positions[0] = chosen;

	Position second_old = backbone->getPosition(Forward ? moveArgs.start_index + 1 : moveArgs.start_index - 1);
	moveArgs.old_positions[1] = second_old;
	moveArgs.new_positions[1] = L;

	// Continue propagation loop as before
	int offset = 2;
	int i = Forward ? moveArgs.start_index + 2 : moveArgs.start_index - 2;

	while ((Forward ? i < backbone->getLength() : i >= 0)) {
		Position curr_pos = backbone->getPosition(i);
		if (config.areAdjacent(curr_pos, moveArgs.new_positions[offset - 1])) break;

		moveArgs.old_positions[offset] = curr_pos;
		moveArgs.new_positions[offset] = Forward ? backbone->getPosition(i - 2) : backbone->getPosition(i + 2);

		moveArgs.num_monomer_moved++;
		if (Forward) ++i; else --i;
		++offset;
	}

	moveArgs.status = MoveStatus::Possible;
}
