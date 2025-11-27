#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**-------------------------- Pull End move ---------------------------**/
void PullEndSquare::checkMove(Configuration<SquareLattice>& config) {
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


bool PullEndSquare::apply(Configuration<SquareLattice>& config) {
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



void PullEndSquare::generateMoveArgs(Configuration<SquareLattice>& config) {
    const size_t bb = dist_backbone(this->rng);
    moveArgs.backbone_id = bb;
    auto backbone = config.getBackbone_ptr(bb);
    b_length = backbone->getLength();

    // need at least two beads to do an end pull
    if (b_length < 2) { moveArgs.status = MoveStatus::NotPossible; return; }

    forward = (dist_end(this->rng) == 0);
    moveArgs.start_index = forward ? 0 : b_length - 1;
    int neighbor_idx = forward ? 1 : b_length - 2;

    Position end_pos = backbone->getPosition(moveArgs.start_index);
    Position neighbor_pos = backbone->getPosition(neighbor_idx);

    // direction from neighbor -> end (unit vector on square lattice)
    int dx = end_pos[0] - neighbor_pos[0];
    int dy = end_pos[1] - neighbor_pos[1];

    // perpendiculars
    Position perp1 = { -dy,  dx };
    Position perp2 = {  dy, -dx };

    // Pick one of the 7 allowed end-pull moves uniformly.
    // Mapping (move_id):
    // 0..2 : A = E + d; B = A + { d, perp1, perp2 }  (3 options)
    // 3..4 : A = E + perp1; B = A + { d, perp1 }     (2 options)
    // 5..6 : A = E + perp2; B = A + { d, perp2 }     (2 options)
    std::uniform_int_distribution<int> dist_move(0, 6);
    int move_id = dist_move(this->rng);

    Position A_raw, B_raw;
    switch (move_id) {
        case 0: // A = E + d; B = A + d
            A_raw = { end_pos[0] + dx,       end_pos[1] + dy };
            B_raw = { A_raw[0]   + dx,       A_raw[1]   + dy };
            break;
        case 1: // A = E + d; B = A + perp1
            A_raw = { end_pos[0] + dx,       end_pos[1] + dy };
            B_raw = { A_raw[0]   + perp1[0], A_raw[1]   + perp1[1] };
            break;
        case 2: // A = E + d; B = A + perp2
            A_raw = { end_pos[0] + dx,       end_pos[1] + dy };
            B_raw = { A_raw[0]   + perp2[0], A_raw[1]   + perp2[1] };
            break;
        case 3: // A = E + perp1; B = A + d
            A_raw = { end_pos[0] + perp1[0], end_pos[1] + perp1[1] };
            B_raw = { A_raw[0]   + dx,       A_raw[1]   + dy };
            break;
        case 4: // A = E + perp1; B = A + perp1
            A_raw = { end_pos[0] + perp1[0], end_pos[1] + perp1[1] };
            B_raw = { A_raw[0]   + perp1[0], A_raw[1]   + perp1[1] };
            break;
        case 5: // A = E + perp2; B = A + d
            A_raw = { end_pos[0] + perp2[0], end_pos[1] + perp2[1] };
            B_raw = { A_raw[0]   + dx,       A_raw[1]   + dy };
            break;
        case 6: // A = E + perp2; B = A + perp2
            A_raw = { end_pos[0] + perp2[0], end_pos[1] + perp2[1] };
            B_raw = { A_raw[0]   + perp2[0], A_raw[1]   + perp2[1] };
            break;
        default:
            moveArgs.status = MoveStatus::NotPossible;
            return;
    }

    // Apply boundary conditions (wrap)
    Position A = config.applyBoundary(A_raw);
    Position B = config.applyBoundary(B_raw);

    // --- Feasibility checks (only for the chosen move) ---
    // 1) A and B must be empty
    if (config.isSiteOccupied(A) || config.isSiteOccupied(B)) {
        moveArgs.status = MoveStatus::NotPossible;
        return;
    }

    // --- Initialize first two monomers (end and its neighbor) ---
    moveArgs.num_monomer_moved = 2;
    moveArgs.old_positions[0] = end_pos;        // end monomer
    moveArgs.new_positions[0] = B_raw;              // new end pos
    moveArgs.old_positions[1] = neighbor_pos;   // neighbor monomer
    moveArgs.new_positions[1] = A_raw;              // new neighbor pos

    // --- Propagate efficiently along backbone ---
    int offset = 2;
    int i = forward ? 2 : b_length - 3;
    while (forward ? (i < b_length) : (i >= 0)) {
        Position curr_pos = backbone->getPosition(i);
        // stop if current is adjacent to the last new position
        if (config.areAdjacent(curr_pos, moveArgs.new_positions[offset - 1])) break;

        moveArgs.old_positions[offset] = curr_pos;
        moveArgs.new_positions[offset] = forward ? backbone->getPosition(i - 2) : backbone->getPosition(i + 2);
        moveArgs.num_monomer_moved++;
        if (forward) ++i; else --i;
        ++offset;
    }

    moveArgs.status = MoveStatus::Possible;
}
