#ifndef MOVESLIST_H
#define MOVESLIST_H

#include "base_move.hpp"
#include "custom_structs.hpp"

//-------------------------- Square moves ------------------------------
// General square moves
#include "general_square_moves/translation_square.hpp"
#include "general_square_moves/symmetry_square.hpp"
#include "general_square_moves/rotation_square.hpp"

// VSDH square moves
#include "VSHD_square/corner_square.hpp"
#include "VSHD_square/crankshaft_square.hpp"
#include "VSHD_square/tail_square.hpp"

// Pull square moves
#include "PULL_square/pull_square.hpp"
#include "PULL_square/pull_end_square.hpp"


//-------------------------- Cubic moves ------------------------------
//~ // General cubic moves
#include "general_cubic_moves/translation_cubic.hpp"
#include "general_cubic_moves/symmetry_cubic.hpp"
#include "general_cubic_moves/rotation_cubic.hpp"

//~ // Pull cubic moves
#include "PULL_cubic/pull_cubic.hpp"
#include "PULL_cubic/pull_end_cubic.hpp"


//-------------------------- Sequence moves ----------------------------
#include "sequence_moves/one_point_mutation.hpp"
#include "sequence_moves/two_point_switch.hpp"

#endif // MOVELIST_H
