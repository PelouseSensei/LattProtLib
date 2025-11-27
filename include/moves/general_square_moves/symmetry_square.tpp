#include "../base_move.hpp"
#include "../custom_structs.hpp"

/**------------------------- Symmetry Square --------------------------**/

void SymmetrySquare::checkMove(Configuration<SquareLattice>& config) {
    auto backbone = config.getBackbone_ptr(moveArgs.backbone_id);

    double dE = 0.0;

    // --- Loop over all monomers in the backbone ---
    for (int i = 0; i < moveArgs.b_length; i++) {
        const auto& old_pos = moveArgs.old_positions[i];
        const auto& new_pos = moveArgs.new_positions[i];

        if (config.isSiteOccupiedIgnoringBackbonei(config.applyBoundary(new_pos),
                                                   moveArgs.backbone_id)) {
            moveArgs.status = MoveStatus::NotPossible;
            return;
        }

        int type = backbone->getType(i);

        moveArgs.old_counts[i] = 0;
        moveArgs.new_counts[i] = 0;

        dE -= config.computeLocalEnergyAndContactsIgnoringBackbonei(
            old_pos, type, moveArgs.old_contacts[i], moveArgs.old_counts[i], moveArgs.backbone_id);

        dE += config.computeLocalEnergyAndContactsIgnoringBackbonei(
            new_pos, type, moveArgs.new_contacts[i], moveArgs.new_counts[i], moveArgs.backbone_id);
    }

    moveArgs.dE = dE;
    moveArgs.status = MoveStatus::Possible;
}


bool SymmetrySquare::apply(Configuration<SquareLattice>& config) {
    int bb_id = moveArgs.backbone_id;

    for (int i = 0; i < moveArgs.b_length; i++) {
        for (int c = 0; c < moveArgs.old_counts[i]; c++) {
            int neighbor_bb = moveArgs.old_contacts[i][c];
            config.contact_matrix.dec(bb_id, neighbor_bb);
        }
    }

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


void SymmetrySquare::generateMoveArgs(Configuration<SquareLattice>& config) {
    const size_t bb = dist_backbone(this->rng);
    moveArgs.backbone_id = bb;

    auto backbone = config.getBackbone_ptr(bb);
    moveArgs.b_length = backbone->getLength();

    // Pick symmetry axis
    moveArgs.symmetry_axis = static_cast<SymmetryAxis>(dist_axis(this->rng));

    for (int i = 0; i < moveArgs.b_length; ++i) {
        const auto& curr = backbone->getPosition(i);
        moveArgs.old_positions[i] = curr;

        Position r = curr;
        switch (moveArgs.symmetry_axis) {
            case SymmetryAxis::X:     r[1] = -r[1]; break;
            case SymmetryAxis::Y:     r[0] = -r[0]; break;
            case SymmetryAxis::DiagXY: std::swap(r[0], r[1]); break;
            case SymmetryAxis::DiagYX: { int tmp = r[0]; r[0] = -r[1]; r[1] = -tmp; } break;
            default: break;
        }

        moveArgs.new_positions[i] = r;
    }
}
