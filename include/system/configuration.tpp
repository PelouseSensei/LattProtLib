#ifndef CONFIGURATION_T
#define CONFIGURATION_T

/**------------------------ Backbone Class --------------------------**/

template<typename LatticeType>
void Backbone<LatticeType>::print() const {
	//~ std::cout << "Backbone details:\n";
	std::cout << "  Dimension: " << Dim << "\n";
	std::cout << "  Length:    " << length_ << "\n";
	std::cout << "  Sequence:  ";
	for (const auto& s : sequence_) {
		std::cout << s << " ";
	}
	std::cout << "\n";

	std::cout << "  Positions:\n";
	for (size_t i = 0; i < positions_.size(); ++i) {
		std::cout << "    [" << i << "]: (";
		for (size_t j = 0; j < positions_[i].size(); ++j) {
			std::cout << positions_[i][j];
			if (j < positions_[i].size() - 1) std::cout << ", ";
		}
		std::cout << ")\n";
	}
}

/**--------------------- Interaction Matrix Class -------------------**/

InterMatrix InterMatrix::fromInterFile(const std::string& filename) {
	std::ifstream file(filename);
	if (!file) {
		throw std::runtime_error("Cannot open file: " + filename);
	}

	std::string line;
	size_t num_types = 0;

	// First line: Number of types = <integer>
	if (!std::getline(file, line)) {
		throw std::runtime_error("Empty file: " + filename);
	}
	{
		std::istringstream iss(line);
		std::string word;
		bool found_eq = false;
		while (iss >> word) {
			if (word == "=") {
				found_eq = true;
				break;
			}
		}
		if (!found_eq || !(iss >> num_types) || num_types == 0) {
			throw std::runtime_error("Invalid header format in .inter file");
		}
	}
	
	// Read the square matrix
	std::vector<double> values;
	values.reserve(num_types * num_types);

	double val;
	while (file >> val) {
		//~ std::cout << val << std::endl;
		values.push_back(val);
	}

	if (values.size() != num_types * num_types) {
		throw std::runtime_error("Matrix size does not match declared number of types");
	}

	// Optional: check symmetry
	for (size_t i = 0; i < num_types; ++i) {
		for (size_t j = 0; j < num_types; ++j) {
			if (values[i * num_types + j] != values[j * num_types + i]) {
				throw std::runtime_error("Matrix is not symmetric at (" +
										  std::to_string(i) + ", " + std::to_string(j) + ")");
			}
		}
	}

	// Construct InterMatrix
	return InterMatrix(num_types, std::move(values));
}


/**------------------------ Configuration Class ---------------------**/

//-------------------------- Private Methods ---------------------------

template<typename LatticeType>
int Configuration<LatticeType>::checkSite2D(const int x, const int y) const {
	return occupancy(x, y);
}

template<typename LatticeType>
int Configuration<LatticeType>::checkSite3D(const int x, const int y, const int z) const {
	return occupancy(x, y, z);
}


template<typename LatticeType>
bool Configuration<LatticeType>::isSiteOccupied2D(const int x, const int y) const {
	return occupancy(x, y) != -1;
}

template<typename LatticeType>
bool Configuration<LatticeType>::isSiteOccupied3D(const int x, const int y, const int z) const {
	return occupancy(x, y, z) != -1;
}


template<typename LatticeType>
bool Configuration<LatticeType>::isSiteOccupiedIgnoringBackbonei2D(const int x, const int y, const size_t ignored_b) const {
	int global_id = occupancy(x, y);
	if (global_id == -1) { return false; }
	else {
		return ignored_b != getBackboneId(global_id);
	}
}

template<typename LatticeType>
bool Configuration<LatticeType>::isSiteOccupiedIgnoringBackbonei3D(const int x, const int y, const int z, const size_t ignored_b) const {
	int global_id = occupancy(x, y, z);
	if (global_id == -1) { return false; }
	else {
		return ignored_b!= getBackboneId(global_id);
	}
}


template<typename LatticeType>
void Configuration<LatticeType>::setOccupancy2D(const int x, const int y, const int global_id) {
	occupancy(x, y) = global_id;
}

template<typename LatticeType>
void Configuration<LatticeType>::setOccupancy3D(const int x, const int y, const int z, const int global_id) {
	occupancy(x, y, z) = global_id;
}


//-------------------------------- Misc --------------------------------

template<typename LatticeType>
double Configuration<LatticeType>::computeLocalEnergy(const Position& position, const int g_id, const int pos_type) {
	double loc_energy = 0.0;
	
	for (Position& displacement : displacements) {
		Position neighbour_pos = position;
		if constexpr (Dim == 2) {
			neighbour_pos[0] += displacement[0];
			neighbour_pos[1] += displacement[1];
		} else if constexpr (Dim == 3) {
			neighbour_pos[0] += displacement[0];
			neighbour_pos[1] += displacement[1];
			neighbour_pos[2] += displacement[2];
		}
		auto global_id = checkSite(applyBoundary(neighbour_pos));
		if (global_id == -1) continue;
		
		if (areConsecutive(global_id, g_id)) continue;
		
		auto [backbone_id, local_id] = getBackboneAndLocalId(global_id);
		int neighbour_type = backbone_list[backbone_id]->getType(local_id);
		loc_energy += matrix(pos_type, neighbour_type);
		
	} 
	return loc_energy;
}

template<typename LatticeType>
double Configuration<LatticeType>::computeLocalEnergyIgnoringBackbonei(const Position& position, const int g_id, const int pos_type, const size_t ignored_b) {
	double loc_energy = 0.0;
	
	for (Position& displacement : displacements) {
		Position neighbour_pos = position;
		if constexpr (Dim == 2) {
			neighbour_pos[0] += displacement[0];
			neighbour_pos[1] += displacement[1];
		} else if constexpr (Dim == 3) {
			neighbour_pos[0] += displacement[0];
			neighbour_pos[1] += displacement[1];
			neighbour_pos[2] += displacement[2];
		}
		auto global_id = checkSite(applyBoundary(neighbour_pos));
		if (global_id == -1 || ignored_b == getBackboneId(global_id)) continue;
		
		//~ if (areConsecutive(global_id, g_id)) continue;
		
		auto [backbone_id, local_id] = getBackboneAndLocalId(global_id);
		int neighbour_type = backbone_list[backbone_id]->getType(local_id);
		loc_energy += matrix(pos_type, neighbour_type);
		
	} 
	return loc_energy;
}

template<typename LatticeType>
double Configuration<LatticeType>::computeLocalEnergySingleCount(const Position& position, const int g_id, const int pos_type) {
	
	double loc_energy = 0.0;

    auto [b_id, local_id] = getBackboneAndLocalId(g_id);

    for (const Position& displacement : displacements) {
        Position neighbour_pos = position;
        if constexpr (Dim == 2) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
        } else if constexpr (Dim == 3) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
            neighbour_pos[2] += displacement[2];
        }

        auto neighbour_g_id = checkSite(applyBoundary(neighbour_pos));
        if (neighbour_g_id != -1) {
            if (areConsecutive(g_id, neighbour_g_id)) {
                continue;
            }

            auto [nb_backbone, nb_local] = getBackboneAndLocalId(neighbour_g_id);
            // Only count intra-chain if neighbour's local index > current
            if (nb_backbone == b_id && nb_local < local_id) {
                continue;
            }

            int neighbour_type = backbone_list[nb_backbone]->getType(nb_local);
            loc_energy += matrix(pos_type, neighbour_type);
        }
    }
    return loc_energy;
}

template<typename LatticeType>
double Configuration<LatticeType>::computeLocalEnergyAndContactsSingleCount(
    const Position& pos,
    const int g_id,
    const int pos_type,
    //~ int backbone_id,
    std::array<int, csts::MAX_TOTAL_CONTACTS>& contacts,
    int& n_contacts) {
	
	double loc_energy = 0.0;
	n_contacts = 0;
	
	auto [b_id, local_id] = getBackboneAndLocalId(g_id);

	for (Position& displacement : displacements) {
        Position neighbour_pos = pos;
        if constexpr (Dim == 2) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
        } else if constexpr (Dim == 3) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
            neighbour_pos[2] += displacement[2];
        }

        auto global_id = checkSite(applyBoundary(neighbour_pos));
        if (global_id == -1) continue;

        if (areConsecutive(global_id, g_id)) continue;

        auto [nb_id, nb_local] = getBackboneAndLocalId(global_id);
        
        if (nb_id == b_id && nb_local < local_id) continue;
                
        int neighbour_type = backbone_list[nb_id]->getType(nb_local);
        loc_energy += matrix(pos_type, neighbour_type);

        // Record contact
        contacts[n_contacts++] = nb_id;
    }

    return loc_energy;
}

template<typename LatticeType>
double Configuration<LatticeType>::computeBlockEnergyAndContacts(
    const int backbone_id,
    const int start_index,
    const int num_monomer_moved,
    const std::array<Position, csts::MAX_BACKBONE_LENGTH>& positions,
    const std::array<int, csts::MAX_BACKBONE_LENGTH>& types,
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH>& contacts,
    std::array<int, csts::MAX_BACKBONE_LENGTH>& n_contacts,
    const bool forward)
{
    double total_energy = 0.0;

    // --- Precompute global IDs for the block and reset contact counts ---
    std::array<int, csts::MAX_BACKBONE_LENGTH> global_ids;
    for (int k = 0; k < num_monomer_moved; ++k) {
        int idx = forward ? start_index + k : start_index - k;
        global_ids[k] = computeGlobalID(backbone_id, idx);
        n_contacts[k] = 0;
    }

    // --- Iterate over the block ---
    for (int k = 0; k < num_monomer_moved; ++k) {
        const Position& pos = positions[k];
        const int pos_type = types[k];
        const int gid = global_ids[k];

        // --- Inline 2D displacement loop to avoid temporary objects ---
        for (const auto& d : displacements) {
            Position neighbour_pos;
            if constexpr (Dim == 2) {
                neighbour_pos[0] = pos[0] + d[0];
                neighbour_pos[1] = pos[1] + d[1];
            } else {
                neighbour_pos[0] = pos[0] + d[0];
                neighbour_pos[1] = pos[1] + d[1];
                neighbour_pos[2] = pos[2] + d[2];
            }

            // --- Apply boundary only once ---
            int nb_gid = checkSite(applyBoundary(neighbour_pos));
            if (nb_gid == -1 || areConsecutive(nb_gid, gid)) continue;

            // --- Move–move pair check ---
            auto [nb_bid, nb_local] = getBackboneAndLocalId(nb_gid);
            bool is_moved = (nb_bid == static_cast<unsigned int>(backbone_id) &&
                             ((forward && nb_local >= start_index && nb_local <= start_index + num_monomer_moved - 1) ||
                              (!forward && nb_local <= start_index && nb_local >= start_index - num_monomer_moved + 1)));

            if (is_moved && ((forward && gid >= nb_gid) || (!forward && gid <= nb_gid))) continue;

            // --- Add energy ---
            int neighbour_type = backbone_list[nb_bid]->getType(nb_local);
            total_energy += matrix(pos_type, neighbour_type);

            // --- Record contact ---
            contacts[k][n_contacts[k]++] = nb_bid;
        }
    }

    return total_energy;
}




template<typename LatticeType>
double Configuration<LatticeType>::computeLocalEnergyAndContacts(
    const Position& pos,
    const int g_id,
    const int pos_type,
    std::array<int, csts::MAX_TOTAL_CONTACTS>& contacts,
    int& n_contacts)
{
    double loc_energy = 0.0;
    n_contacts = 0;

    for (Position& displacement : displacements) {
        Position neighbour_pos = pos;
        if constexpr (Dim == 2) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
        } else if constexpr (Dim == 3) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
            neighbour_pos[2] += displacement[2];
        }

        auto global_id = checkSite(applyBoundary(neighbour_pos));
        if (global_id == -1) continue;

        if (areConsecutive(global_id, g_id)) continue;
        //~ std::cout << "Contact between " << global_id << " and " << g_id << std::endl;

        auto [bb_id, local_id] = getBackboneAndLocalId(global_id);
        int neighbour_type = backbone_list[bb_id]->getType(local_id);
        loc_energy += matrix(pos_type, neighbour_type);

        // Record contact
        contacts[n_contacts++] = bb_id;
    }

    return loc_energy;
}

template<typename LatticeType>
double Configuration<LatticeType>::computeLocalEnergyAndContactsIgnoringBackbonei(
    const Position& pos,
    const int pos_type,
    std::array<int, csts::MAX_TOTAL_CONTACTS>& contacts,
    int& n_contacts,
    const size_t ignored_b)
{
    double loc_energy = 0.0;
    n_contacts = 0;

    for (Position& displacement : displacements) {
        Position neighbour_pos = pos;
        if constexpr (Dim == 2) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
        } else if constexpr (Dim == 3) {
            neighbour_pos[0] += displacement[0];
            neighbour_pos[1] += displacement[1];
            neighbour_pos[2] += displacement[2];
        }

        auto global_id = checkSite(applyBoundary(neighbour_pos));
		if (global_id == -1 || ignored_b == getBackboneId(global_id)) continue;
		
        auto [bb_id, local_id] = getBackboneAndLocalId(global_id);
        int neighbour_type = backbone_list[bb_id]->getType(local_id);
        loc_energy += matrix(pos_type, neighbour_type);

        // Record contact
        contacts[n_contacts++] = bb_id;
    }

    return loc_energy;
}

template<typename LatticeType>
void Configuration<LatticeType>::printContactMatrix() const {
	std::cout << "Contact Matrix (upper-triangle + diagonal):\n";
    for (int i = 0; i < num_backbones; ++i) {
        for (int j = 0; j < num_backbones; ++j) {
            std::cout << contact_matrix.get(i, j) << "\t";  // uses your operator() to access
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

template<typename LatticeType>
void Configuration<LatticeType>::describe() const {
	std::cout << "---> Lattice: " << lattice.getLatticeType() << ", ";
	std::cout << "Dimension: " << Dim << ", ";
	std::cout << "Configuration energy: " << energy << ", ";
	std::cout << "Lattice size: " << L << ", ";
	std::cout << "Number of backbones: " << num_backbones << "\n\n";// ", ";
	//~ std::cout << "Total backbones length: " << total_backbones_length << ", ";
	//~ std::cout << "Cumulative lengths : ";
	//~ for (size_t i = 0; i<backbones_cumulative_lengths.size(); i++) {
	 //~ std::cout << backbones_cumulative_lengths[i] << ", ";
	//~ }
	//~ std::cout << "\n\n";
	for (int i = 0; i<num_backbones; i++ ) {
		std::cout << "Backbone num  " << i << ":\n";
		backbone_list[i]->print();
		std::cout << "\n";
	}
	printContactMatrix();
}

//~ template<typename LatticeType>
//~ bool Configuration<LatticeType>::areAdjacent(const Position& pos1, const Position pos2) const {
	//~ if constexpr (Dim == 2) {
		//~ if (std::abs(pos1[0] - pos2[0]) > 1) { return false; }
		//~ if (std::abs(pos1[1] - pos2[1]) > 1) { return false; }
		
		//~ return true;
	//~ } else if constexpr (Dim == 3) {
		//~ std::cout << "test" << std::endl;
		//~ if (std::abs(pos1[0] - pos2[0]) > 1) { return false; }
		//~ std::cout << "test2" << std::endl;
		//~ if (std::abs(pos1[1] - pos2[1]) > 1) { return false; }
		//~ std::cout << "test3" << std::endl;
		//~ if (std::abs(pos1[2] - pos2[2]) > 1) { return false; }
		//~ std::cout << "test4" << std::endl;
		
		//~ return true;
	//~ }
//~ }

template<typename LatticeType>
bool Configuration<LatticeType>::areAdjacent(const Position& pos1, const Position& pos2) const {
    int diff_count = 0;
    int diff_value = 0;

    for (size_t i = 0; i < Dim; ++i) {
        int d = pos1[i] - pos2[i];
        if (d != 0) {
            ++diff_count;
            diff_value = d;
        }
    }

    // Adjacent if exactly one coordinate differs by ±1
    return (diff_count == 1) && (diff_value == 1 || diff_value == -1);
}


template<typename LatticeType>
int Configuration<LatticeType>::computeGlobalID(const size_t backbone_id, const int monomer_id) const {
	return backbone_id * csts::MAX_BACKBONE_LENGTH + monomer_id;
}

template<typename LatticeType>
std::pair<size_t, int> Configuration<LatticeType>::getBackboneAndLocalId(int global_id) const { // No checks for -1 values, needs to be checked before
    return { global_id / csts::MAX_BACKBONE_LENGTH, global_id % csts::MAX_BACKBONE_LENGTH };
}

template<typename LatticeType>
size_t Configuration<LatticeType>::getBackboneId(const int global_id) const {
	return global_id / csts::MAX_BACKBONE_LENGTH;
}

template<typename LatticeType>
void Configuration<LatticeType>::loadFromCIFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open bb file: " + filename);

    auto readValueAfterEquals = [](std::istream& is, auto& value) {
        std::string token;
        while (is >> token) {
            if (token == "=") {
                if (!(is >> value)) {
                    throw std::runtime_error("Failed to read value after '='");
                }
                return;
            }
        }
        throw std::runtime_error("No '=' found in expected line");
    };

    std::string line;
    size_t num_backbones = 0;

    // --- Read number of backbones ---
    if (!std::getline(file, line)) throw std::runtime_error("Empty bb file");
    {
        std::istringstream iss(line);
        readValueAfterEquals(iss, num_backbones);
        if (num_backbones == 0) throw std::runtime_error("Invalid number of backbones");
    }

    for (size_t b = 0; b < num_backbones; ++b) {
        int b_length = 0;

        // --- Read backbone length ---
        if (!std::getline(file, line)) throw std::runtime_error("Unexpected end of file");
        {
            std::istringstream iss(line);
            readValueAfterEquals(iss, b_length);
            if (b_length <= 0) throw std::runtime_error("Invalid backbone length");
        }

        // --- Read positions ---
        if (!std::getline(file, line)) throw std::runtime_error("Unexpected end of file");
        std::vector<typename LatticeType::Position> positions(b_length);
        {
            std::istringstream iss(line);
            std::string tmp;
            while (iss >> tmp && tmp != "=") {}
            for (int i = 0; i < b_length; ++i) {
                for (size_t d = 0; d < LatticeType::dimension; ++d) {
                    if (!(iss >> positions[i][d])) {
                        throw std::runtime_error("Not enough coordinates for backbone positions");
                    }
                }
            }
        }

        // --- Read sequence ---
        if (!std::getline(file, line)) throw std::runtime_error("Unexpected end of file");
        std::vector<int> sequence(b_length);
        {
            std::istringstream iss(line);
            std::string tmp;
            while (iss >> tmp && tmp != "=") {}
            for (int i = 0; i < b_length; ++i) {
                if (!(iss >> sequence[i])) {
                    throw std::runtime_error("Backbone sequence length mismatch");
                }
            }
        }

        // --- Add backbone ---
        auto backbone = std::make_shared<Backbone<LatticeType>>(b_length, positions, sequence);
        if (!addBackboneToConfiguration(backbone)) {
            throw std::runtime_error("Failed to add backbone " + std::to_string(b));
        }
    }
}




//---------------------------- General Getters -------------------------



//---------------------------- General Setters -------------------------

//~ template<typename LatticeType>
//~ bool Configuration<LatticeType>::moveMonomer(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position ) {
	
	//~ if (isSiteOccupied(new_position)) { return false; }

	//~ // Do not check chain continuity before moving to new position
	//~ auto global_id = computeGlobalID(backbone_id, monomer_id);
	
	//~ int type = backbone_list[backbone_id]->getType(monomer_id);
	//~ energy -= computeLocalEnergy(old_position, global_id, type);
	
	//~ setOccupancy( new_position, global_id);
	//~ setOccupancy( old_position, -1);
	//~ backbone_list[backbone_id]->getPosition(monomer_id) = new_position;
	
	//~ energy += computeLocalEnergy(new_position, global_id, type);
	
	
		
	//~ return true;
//~ }

template<typename LatticeType>
bool Configuration<LatticeType>::moveMonomer(const size_t backbone_id,
                                            const int monomer_id,
                                            const Position& old_position,
                                            const Position& new_position)
{
    if (isSiteOccupied(new_position)) return false;

    auto global_id = computeGlobalID(backbone_id, monomer_id);
    int type = backbone_list[backbone_id]->getType(monomer_id);

    // --- 1. Remove old contacts and compute old energy ---
    std::array<int, csts::MAX_TOTAL_CONTACTS> old_contacts{};
    int n_old_contacts = 0;
    double old_energy = computeLocalEnergyAndContacts(old_position, global_id, type,
                                                      old_contacts, n_old_contacts);

    for (int c = 0; c < n_old_contacts; c++) {
        int neighbor_bb = old_contacts[c];
        contact_matrix.dec(backbone_id, neighbor_bb);
    }

    // --- 2. Update occupancy and position ---
    setOccupancy(new_position, global_id);
    setOccupancy(old_position, -1);
    backbone_list[backbone_id]->getPosition(monomer_id) = new_position;

    // --- 3. Add new contacts and compute new energy ---
    std::array<int, csts::MAX_TOTAL_CONTACTS> new_contacts{};
    int n_new_contacts = 0;
    double new_energy = computeLocalEnergyAndContacts(new_position, global_id, type,
                                                      new_contacts, n_new_contacts);

    for (int c = 0; c < n_new_contacts; c++) {
        int neighbor_bb = new_contacts[c];
        contact_matrix.inc(backbone_id, neighbor_bb);
    }

    // --- 4. Update energy ---
    energy += new_energy - old_energy;

    return true;
}


template<typename LatticeType>
bool Configuration<LatticeType>::moveMonomerIgnoringBackbonei(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position, const size_t ignored_b ) {
	
	if (isSiteOccupiedIgnoringBackbonei(new_position, ignored_b)) { return false; }
	
	// Do not check chain continuity before moving to new position
	auto global_id = computeGlobalID(backbone_id, monomer_id);
	
	int type = backbone_list[backbone_id]->sequence[monomer_id];
	energy -= computeLocalEnergy(old_position, global_id, type);
	
	setOccupancy( new_position, global_id);
	setOccupancy( old_position, -1);
	backbone_list[backbone_id]->positions[monomer_id] = new_position;
	
	energy += computeLocalEnergy(new_position, global_id, type);
		
	return true;
}

template<typename LatticeType>
void Configuration<LatticeType>::moveMonomerMinimal(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position, const double new_E ) {
	
	// Do not check chain continuity before moving to new position
	auto global_id = computeGlobalID(backbone_id, monomer_id);
	
	// Move monomer
	backbone_list[backbone_id]->getPosition(monomer_id) = new_position;
	// Update lattice occupancy
	setOccupancy(applyBoundary(new_position), global_id);
	setOccupancy(applyBoundary(old_position), -1);
	// Update energy
	energy = new_E;
}

template<typename LatticeType>
void Configuration<LatticeType>::moveMonomerOnly(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position ) {
	
	// Do not check chain continuity before moving to new position
	auto global_id = computeGlobalID(backbone_id, monomer_id);
	
	// Move monomer
	backbone_list[backbone_id]->getPosition(monomer_id) = new_position;
	// Update lattice occupancy
	setOccupancy(applyBoundary(new_position), global_id);
	setOccupancy(applyBoundary(old_position), -1);
}

template<typename LatticeType>
void Configuration<LatticeType>::moveBackboneMinimal(const size_t backbone_id, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& old_positions, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& new_positions, const int b_length, const double new_E) {
    //~ int length = old_positions.size();
    auto backbone = backbone_list[backbone_id];

    // 1. Clear old occupancies
    for (int i = 0; i < b_length; i++) {
        //~ auto gid = computeGlobalID(backbone_id, i);
        setOccupancy(applyBoundary(old_positions[i]), -1);
    }

    // 2. Set new occupancies + update positions
    for (int i = 0; i < b_length; i++) {
        auto gid = computeGlobalID(backbone_id, i);
        backbone->getPosition(i) = new_positions[i];
        setOccupancy(applyBoundary(new_positions[i]), gid);
    }

    // 3. Update energy
    energy = new_E;
}

template<typename LatticeType>
void Configuration<LatticeType>::moveBackboneMinimalPullBackward(const size_t backbone_id, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& old_positions, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& new_positions, const int start_index, const int b_length, const double new_E) {
    //~ int length = old_positions.size();
    auto backbone = backbone_list[backbone_id];

    // 1. Clear old occupancies
    for (int i = 0; i < b_length; i++) {
        //~ auto gid = computeGlobalID(backbone_id, i);
        setOccupancy(applyBoundary(old_positions[i]), -1);
    }

    // 2. Set new occupancies + update positions
    for (int i = 0; i < b_length; i++) {
        auto gid = computeGlobalID(backbone_id, start_index - i);
        //~ std::cout << "gid : " << gid << std::endl;
        backbone->getPosition(start_index - i) = new_positions[i];
        setOccupancy(applyBoundary(new_positions[i]), gid);
    }

    // 3. Update energy
    energy = new_E;
}


template<typename LatticeType>
void Configuration<LatticeType>::moveBackboneMinimalPullForward(const size_t backbone_id, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& old_positions, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& new_positions, const int start_index, const int b_length, const double new_E) {
    //~ int length = old_positions.size();
    auto backbone = backbone_list[backbone_id];

    // 1. Clear old occupancies
    for (int i = 0; i < b_length; i++) {
        //~ auto gid = computeGlobalID(backbone_id, i);
        setOccupancy(applyBoundary(old_positions[i]), -1);
    }

    // 2. Set new occupancies + update positions
    for (int i = 0; i < b_length; i++) {
        auto gid = computeGlobalID(backbone_id, start_index + i);
        //~ std::cout << "gid : " << gid << std::endl;
        backbone->getPosition(start_index + i) = new_positions[i];
        setOccupancy(applyBoundary(new_positions[i]), gid);
    }

    // 3. Update energy
    energy = new_E;
}




template<typename LatticeType>
void Configuration<LatticeType>::mutateMonomer(const size_t backbone_id, const int monomer_id, const int new_type, const double new_E ) {
	auto backbone = getBackbone_ptr(backbone_id);
	backbone->setType(monomer_id, new_type);
	energy = new_E;
}

//-------------------------- Occupancy Tensors -------------------------



//-------------------------- Backbones gestion -------------------------

template<typename LatticeType>
bool Configuration<LatticeType>::addBackboneToConfiguration( std::shared_ptr<Backbone<LatticeType>> backbone) {
	
	int b_id     = num_backbones;
	int b_length = backbone->getLength();
	
	for (int i = 0; i<b_length; i++) { // Try to add the backbone
		Position b_pos = applyBoundary(backbone->getPosition(i));
		if (!isSiteOccupied(b_pos)) {
			setOccupancy(b_pos, computeGlobalID(b_id, i));
		}
		else {
			for (int j = i - 1; j >= 0; --j) { // If site not empty, addition not possible -> revert all modifications of the occupancy tensor
				setOccupancy(applyBoundary(backbone->getPosition(j)), -1);
			}
			return false;
		}
	}
	
	backbone_list.push_back(backbone);
	total_backbones_length += b_length;
	//~ backbones_cumulative_lengths.push_back(total_backbones_length);
	num_backbones++;
	
	int n_contacts = 0;
	std::array<int, csts::MAX_TOTAL_CONTACTS> contacts{};
	
	// Update energy
	//~ std::cout << "Adding backbone " << num_backbones << std::end;
	for (int i = 0; i<b_length; i++) { // Try to add the backbone
		energy += computeLocalEnergyAndContactsSingleCount(backbone->getPosition(i), computeGlobalID(b_id, i), backbone->getType(i), contacts, n_contacts);
		
		// increment new contacts
		for (int i = 0; i < n_contacts; ++i) {
			int nb = contacts[i];
			contact_matrix.inc(b_id, nb);
		}
	}
	
	
	return true;
}

template<typename LatticeType>    
bool Configuration<LatticeType>::removeBackboneFromConfiguration( const int backbone_id ) {
	
	auto backbone = backbone_list[backbone_id];
	size_t b_length = backbone->getLength();
	
	int n_contacts = 0;
	std::array<int, csts::MAX_TOTAL_CONTACTS> contacts{};
	
	// Update energy
	for (size_t i = 0; i<b_length; i++) {
		auto pos_to_remove = applyBoundary(backbone->getPosition(i));
		energy -= computeLocalEnergyAndContactsSingleCount(backbone->getPosition(i), computeGlobalID(num_backbones-1, i), backbone->getType(i), contacts, n_contacts);
		
		// increment new contacts
		for (int i = 0; i < n_contacts; ++i) {
			int nb = contacts[i];
			contact_matrix.dec(num_backbones-1, nb);
		}
		setOccupancy(pos_to_remove, -1);
	}
	
	backbone_list.erase(backbone_list.begin() + backbone_id);
	total_backbones_length -= b_length;
	num_backbones--;
	
	return true;
}

#endif // CONFIGURATION_T
