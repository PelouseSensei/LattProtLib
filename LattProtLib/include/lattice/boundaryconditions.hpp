#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "lattice.hpp"

// ----------------ENUM-----------------------------

enum class BoundaryType {
	Periodic
};




// ------------------------- BASE CLASS ------------------------- 

template<typename LatticeType>
class BoundaryConditionsBase {
public:

	using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    virtual ~BoundaryConditionsBase() = default;
    virtual Position applyBoundaryConditions(const int lattice_size, const Position& position) const = 0;
    virtual BoundaryType getBoundaryCondition() const = 0;
    virtual int getDimension() const = 0;
};




// ------------------------- DERIVED CLASSES -------------------------

template<typename LatticeType>
class PeriodicBoundary : public BoundaryConditionsBase<LatticeType> {
public:

	using Position = typename BoundaryConditionsBase<LatticeType>::Position;
	
    Position applyBoundaryConditions(const int lattice_size, const Position& position) const override {
		auto pos_ = position;

		if constexpr (this->Dim == 2) {
			pos_[0] = ((pos_[0] % lattice_size) + lattice_size) % lattice_size;
			pos_[1] = ((pos_[1] % lattice_size) + lattice_size) % lattice_size;
		} else if constexpr (this->Dim == 3) {
			pos_[0] = ((pos_[0] % lattice_size) + lattice_size) % lattice_size;
			pos_[1] = ((pos_[1] % lattice_size) + lattice_size) % lattice_size;
			pos_[2] = ((pos_[2] % lattice_size) + lattice_size) % lattice_size;
		}

		return pos_;
	}

    
    BoundaryType getBoundaryCondition() const override { return BoundaryType::Periodic; }
    
    int getDimension() const override {
		if constexpr (this->Dim == 2) { return 2; }
		if constexpr (this->Dim == 3) { return 3; }
	}
};



// ------------------------- FACTORY METHOD -------------------------

template<typename LatticeType>
std::unique_ptr<BoundaryConditionsBase<LatticeType>> create_BC(BoundaryType type) {
    switch (type) {
        case BoundaryType::Periodic:
			return std::make_unique<PeriodicBoundary<LatticeType>>();
        default:
            throw std::invalid_argument("Unknown boundary type");
    }
}



#endif // BOUNDARYCONDITIONS_H
