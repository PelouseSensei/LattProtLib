#ifndef CUSTOM_STRUCTS_H
#define CUSTOM_STRUCTS_H

/**--------------------------- Custom Structs -----------------------**/

enum class MoveStatus {
	Possible,
	NotPossible
};

enum class SymmetryAxis {
    None,
    X,       // reflection across X-axis
    Y,       // reflection across Y-axis
    DiagXY,  // reflection across y = x
    DiagYX   // reflection across y = -x
};

enum class SymmetryAxis3D {
    None,
    X,    // reflection across YZ-plane (flip x -> -x)
    Y,    // reflection across XZ-plane (flip y -> -y)
    Z,    // reflection across XY-plane (flip z -> -z)
    XY,   // reflection across plane x=y (swap x,y)
    XZ,   // reflection across plane x=z (swap x,z)
    YZ    // reflection across plane y=z (swap y,z)
    // You could add diagonals like x=y=z if you want later
};

enum class RotationAxis {
    X,  // rotate around X-axis
    Y,  // rotate around Y-axis
    Z   // rotate around Z-axis
};



/**---------------------- Move arguments structs --------------------**/

template<typename LatticeType>
struct OneMonoMoveArgs {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;
    Position old_position{};
    Position new_position{};

    size_t backbone_id;
    int monomer_id;
    int type;
    double dE;

    std::array<int, csts::MAX_TOTAL_CONTACTS> old_contacts{};
    std::array<int, csts::MAX_TOTAL_CONTACTS> new_contacts{};
    int old_count = 0;
    int new_count = 0;
};




template<typename LatticeType>
struct TwoMonoMoveArgs {
	
	using Position = typename LatticeType::Position;
	static constexpr size_t Dim = LatticeType::dimension;
	
	MoveStatus status;
    Position old_position1{};
    Position new_position1{};
    Position old_position2{};
    Position new_position2{};
    
    Position pos_mono_0{};
    Position pos_mono_3{};
	
	size_t backbone_id1;
	size_t backbone_id2;
	int monomer_id1;
	int monomer_id2;
	int type1;
	int type2;
	double dE;
	
	std::array<int, csts::MAX_TOTAL_CONTACTS> old_contacts1{};
	std::array<int, csts::MAX_TOTAL_CONTACTS> old_contacts2{};
    std::array<int, csts::MAX_TOTAL_CONTACTS> new_contacts1{};
    std::array<int, csts::MAX_TOTAL_CONTACTS> new_contacts2{};
    int old_count1 = 0;
    int old_count2 = 0;
    int new_count1 = 0;
    int new_count2 = 0;
};



template<typename LatticeType>
struct TranslationMoveArgs {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;
    Position translation_vec{};
	std::array<Position, csts::MAX_BACKBONE_LENGTH> old_positions{};
	std::array<Position, csts::MAX_BACKBONE_LENGTH> new_positions{};

    size_t backbone_id;
    int b_length;
    double dE;

    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> old_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> old_counts{};
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> new_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> new_counts{};
};

template<typename LatticeType>
struct RotMoveArgs {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;

    // Rotation info
    int rotation_angle{};      // rotation angle (e.g., 90, 180, 270 degrees)
    int pivot_index;
    Position rotation_center{}; // pivot point for the rotation

    // Geometry before/after move
    std::array<Position, csts::MAX_BACKBONE_LENGTH> old_positions{};
    std::array<Position, csts::MAX_BACKBONE_LENGTH> new_positions{};

    size_t backbone_id{};
    int b_length{};
    double dE{};

    // Contacts before/after move
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> old_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> old_counts{};
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> new_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> new_counts{};
};

template<typename LatticeType>
struct RotMoveArgs3D {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;

    // Rotation info
    int rotation_angle{};      // rotation angle (e.g., 90, 180, 270 degrees)
    int pivot_index;
    Position rotation_center{}; // pivot point for the rotation
    RotationAxis rotation_axis{}; // Axis of the rotation

    // Geometry before/after move
    std::array<Position, csts::MAX_BACKBONE_LENGTH> old_positions{};
    std::array<Position, csts::MAX_BACKBONE_LENGTH> new_positions{};

    size_t backbone_id{};
    int b_length{};
    double dE{};

    // Contacts before/after move
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> old_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> old_counts{};
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> new_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> new_counts{};
};


template<typename LatticeType>
struct SymMoveArgs {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;

    // Symmetry info
    SymmetryAxis symmetry_axis{SymmetryAxis::None};

    // Geometry before/after move
    std::array<Position, csts::MAX_BACKBONE_LENGTH> old_positions{};
    std::array<Position, csts::MAX_BACKBONE_LENGTH> new_positions{};

    size_t backbone_id{};
    int b_length{};
    double dE{};

    // Contacts before/after move
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> old_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> old_counts{};
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> new_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> new_counts{};
};


template<typename LatticeType>
struct SymMoveArgs3D {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;

    // Symmetry info
    SymmetryAxis3D symmetry_axis{SymmetryAxis3D::None};

    // Geometry before/after move
    std::array<Position, csts::MAX_BACKBONE_LENGTH> old_positions{};
    std::array<Position, csts::MAX_BACKBONE_LENGTH> new_positions{};

    size_t backbone_id{};
    int b_length{};
    double dE{};

    // Contacts before/after move
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> old_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> old_counts{};
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> new_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> new_counts{};
};

template<typename LatticeType>
struct PullMoveArgs {
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    MoveStatus status;

    size_t backbone_id{};
    //~ int b_length{};
    int num_monomer_moved;
    double dE{};

    int start_index{}; // monomer chosen to pull
    std::array<Position, csts::MAX_BACKBONE_LENGTH> old_positions{};
    std::array<Position, csts::MAX_BACKBONE_LENGTH> new_positions{};
    
    std::array<int, csts::MAX_BACKBONE_LENGTH> types{};
	std::array<int, csts::MAX_BACKBONE_LENGTH> global_ids{};


    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> old_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> old_counts{};
    std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH> new_contacts{};
    std::array<int, csts::MAX_BACKBONE_LENGTH> new_counts{};
};



template<typename LatticeType>
struct OneTypeMoveArgs {
	using Position = typename LatticeType::Position;
	static constexpr size_t Dim = LatticeType::dimension;
	
	MoveStatus status;
	Position position{};
	
	size_t backbone_id;
	int monomer_id;
	double dE;
	
	int old_type;
	int new_type;
};

template<typename LatticeType>
struct TwoTypeMoveArgs {
	using Position = typename LatticeType::Position;
	static constexpr size_t Dim = LatticeType::dimension;
	
	MoveStatus status;
	Position position1{};
	Position position2{};
	
	size_t backbone_id1;
	size_t backbone_id2;
	int monomer_id1;
	int monomer_id2;
	double dE;
	
	int old_type1;
	int old_type2;
};

#endif // CUSTOM_STRUCTS_H
