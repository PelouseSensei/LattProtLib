#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../lattice/lattice.hpp"
#include "../lattice/boundaryconditions.hpp"
#include "occupancygrid.hpp"



//--------------------------- Constants --------------------------------

namespace csts {
	const int MAX_BACKBONE_LENGTH = 1000;
	const size_t MAX_TOTAL_MONOMERS = 100;
	const int MAX_NUM_BACKBONE = 10;
	const int MAX_TOTAL_CONTACTS = 20;
}


//------------------------ Backbone Class ------------------------------


template <typename LatticeType>
class Backbone {
public:
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    Backbone(int length,
             std::vector<Position> positions,
             std::vector<int> sequence)
        : length_(length),
          sequence_(std::move(sequence)),
          positions_(std::move(positions))
    {
        if (sequence_.size() != static_cast<size_t>(length_) ||
            positions_.size() != static_cast<size_t>(length_)) {
            throw std::invalid_argument("Sequence or position length mismatch with backbone length.");
        }
    }

    int getLength() const { return length_; }

    const Position& getPosition(int i) const { return positions_[i]; }
    Position& getPosition(int i) { return positions_[i]; }
    
    int getType(int i) const { return sequence_[i]; }
    
    void setType(const int monomer_id, const int new_type) { sequence_[monomer_id] = new_type; }

    void print() const;

    void fillFromFile(const std::string& filename);

private:
    int length_;
    std::vector<int> sequence_;
    std::vector<Position> positions_;
};





//---------------------- Interaction Matrix Class ----------------------

class InterMatrix {
public:
    InterMatrix(size_t num_types, std::initializer_list<double> values)
        : num_types(num_types), matrix(values) {
        if (values.size() != num_types * num_types) {
            throw std::invalid_argument("Initializer size does not match matrix dimensions.");
        }
    }

    InterMatrix(size_t num_types, std::vector<double> values)
        : num_types(num_types), matrix(std::move(values)) {
        if (matrix.size() != num_types * num_types) {
            throw std::invalid_argument("Vector size does not match matrix dimensions.");
        }
    }

    double& operator()(size_t i, size_t j) {
        return matrix[i * num_types + j];
    }

    const double& operator()(size_t i, size_t j) const {
        return matrix[i * num_types + j];
    }

    void print() const {
        for (size_t i = 0; i < num_types; ++i) {
            for (size_t j = 0; j < num_types; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    size_t size() const { return num_types; }

    static InterMatrix fromInterFile(const std::string& filename);

private:
    size_t num_types;
    std::vector<double> matrix;
};






//-------------------------- Contact Matrix ----------------------------

struct ContactMatrix {
    std::array<uint16_t, (csts::MAX_NUM_BACKBONE*(csts::MAX_NUM_BACKBONE+1))/2> data{};

    static constexpr int index(int i, int j) noexcept {
        if (i > j) std::swap(i,j);
        return (i*(2*csts::MAX_NUM_BACKBONE - i + 1))/2 + (j-i);
    }

    void inc(int i, int j) noexcept { data[index(i,j)]++; }
    void dec(int i, int j) noexcept { data[index(i,j)]--; }
    uint16_t get(int i, int j) const noexcept { return data[index(i,j)]; }
};







//------------------------ Configuration Class -------------------------
template<typename LatticeType>
class Configuration {
public:
	//--------------------- Public Attributes --------------------------
	ContactMatrix contact_matrix;
	
	
	//------------------------- Constructor ----------------------------
    
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    Configuration(int L, const BoundaryType BC_type, InterMatrix matrix) : 
		L(L), 
		lattice(), 
		occupancy(L), 
		displacements(lattice.getDisplacements()), 
		
		boundary_conditions(create_BC<LatticeType>(BC_type)),
		
		energy(0.0),
		
		total_backbones_length(0),
		num_backbones(0),
		
		matrix(matrix)
	{}
	
	//------------------------------ Misc ------------------------------
	
    void describe() const;
    
    void printContactMatrix() const;
    
    double computeLocalEnergy(const Position& position, const int g_id, const int pos_type);
    
    double computeLocalEnergyIgnoringBackbonei(const Position& position, const int g_id, const int pos_type, const size_t ignored_b);
    
    double computeLocalEnergySingleCount(const Position& position, const int g_id, const int pos_type);
    
    double computeBlockEnergyAndContacts(const int backbone_id, const int start_index, const int num_monomer_moved, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& positions, const std::array<int, csts::MAX_BACKBONE_LENGTH>& types, std::array<std::array<int, csts::MAX_TOTAL_CONTACTS>, csts::MAX_BACKBONE_LENGTH>& contacts, std::array<int, csts::MAX_BACKBONE_LENGTH>& n_contacts, const bool forward)   ;
    
    double computeLocalEnergyAndContacts(const Position& pos, const int g_id, const int pos_type, std::array<int, csts::MAX_TOTAL_CONTACTS>& contacts, int& n_contacts);
    
    double computeLocalEnergyAndContactsIgnoringBackbonei(const Position& pos, const int pos_type, std::array<int, csts::MAX_TOTAL_CONTACTS>& contacts, int& n_contacts, const size_t ignored_b);
    
    double computeLocalEnergyAndContactsSingleCount(const Position& pos, const int g_id, const int pos_type, std::array<int, csts::MAX_TOTAL_CONTACTS>& contacts, int& n_contacts);
    
    void printOccupancy() const {
		if constexpr (Dim == 2) {
			occupancy.print2D();
		} else if constexpr (Dim == 3) {
			occupancy.print3D();
		} else {
			throw std::domain_error("Unsupported dimension");
		}
	}

    
    bool areConsecutive(int g_id1, int g_id2) const {
		int b1 = g_id1 / csts::MAX_BACKBONE_LENGTH;
		int b2 = g_id2 / csts::MAX_BACKBONE_LENGTH;
		if (b1 != b2) return false;
		return std::abs((g_id1 % csts::MAX_BACKBONE_LENGTH) - (g_id2 % csts::MAX_BACKBONE_LENGTH)) == 1;
	}
    
	std::pair<size_t, int> getBackboneAndLocalId(int global_id) const;
	
	size_t getBackboneId(const int global_id) const;
	
	int computeGlobalID(const size_t backbone_id, const int monomer_id) const;
	
	bool areAdjacent(const Position& pos1, const Position& pos2) const;
	
    void loadFromCIFile(const std::string& filename);
    
    int getTotalLength() { return total_backbones_length; }
    



    //-------------------------- General Getters -----------------------
    
    int getDimension() const { return Dim; }
    
    std::string getLatticeType() const { return lattice->getLatticeType(); }
    
    int getLatticeSize() const { return L;}
    
    BoundaryType getBoundaryType() const { return boundary_conditions->getBoundaryCondition(); }

    const std::vector<Position>& getDisplacements() const {
        return displacements;
    }
    
    double getEnergy() const { return energy; }
    
    std::vector<std::shared_ptr<Backbone<LatticeType>>>& getBackbones() {
		return backbone_list;
	}
	
	int getNumBackbones() const { return num_backbones; }
    



	//-------------------------- General Setters -----------------------
	
    Position applyBoundary(const Position& position) const {
        return boundary_conditions->applyBoundaryConditions(L, position);
    }
    
    bool moveMonomer(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position );
    
    void moveMonomerMinimal(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position, const double new_E );
    
    void moveMonomerOnly(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position );
    
    bool moveMonomerIgnoringBackbonei(const size_t backbone_id, const int monomer_id, const Position& old_position, const Position& new_position, const size_t ignored_b );
    
    void moveBackboneMinimal(const size_t backbone_id, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& old_positions, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& new_positions, const int b_length, const double new_E);
    
    void moveBackboneMinimalPullBackward(const size_t backbone_id, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& old_positions, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& new_positions, const int start_index, const int b_length, const double new_E);
    
    void moveBackboneMinimalPullForward(const size_t backbone_id, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& old_positions, const std::array<Position, csts::MAX_BACKBONE_LENGTH>& new_positions, const int start_index, const int b_length, const double new_E);
    
    void mutateMonomer(const size_t backbone_id, const int monomer_id, const int new_type, const double new_E );




    //------------------ Occupancy Tensors -----------------------------
    
    int checkSite(const Position& position) const {
		if constexpr (Dim ==2) {
			return checkSite2D(position[0], position[1]);
		} else if constexpr (Dim == 3) {
			return checkSite3D(position[0], position[1], position[2]);
		}
	}
    
    bool isSiteOccupied(const Position& position) const {
		if constexpr (Dim ==2) {
			return isSiteOccupied2D(position[0], position[1]);
		} else if constexpr (Dim == 3) {
			return isSiteOccupied3D(position[0], position[1], position[2]);
		}
	}
    
    bool isSiteOccupiedIgnoringBackbonei(const Position& position, const size_t backbone_id) const {
		if constexpr (Dim ==2) {
			return isSiteOccupiedIgnoringBackbonei2D(position[0], position[1], backbone_id);
		} else if constexpr (Dim == 3) {
			return isSiteOccupiedIgnoringBackbonei3D(position[0], position[1], position[2], backbone_id);
		}
	}
    
    void setOccupancy(const Position& pos, const int global_id) {
		if constexpr (Dim ==2) {
			return setOccupancy2D(pos[0], pos[1], global_id);
		} else if constexpr (Dim == 3) {
			return setOccupancy3D(pos[0], pos[1], pos[2], global_id);
		}
	}
    
    
    

    //------------------------ Backbones gestion -----------------------
    
    std::shared_ptr<Backbone<LatticeType>> getBackbone_ptr(const int backbone_id) { return backbone_list[backbone_id]; }
    
    bool addBackboneToConfiguration( std::shared_ptr<Backbone<LatticeType>> backbone);
    
    bool removeBackboneFromConfiguration( const int backbone_id );
    



private:
	
	//--------------------------- Attributes ---------------------------	
	int L;
    LatticeType lattice;
    OccupancyGrid<LatticeType> occupancy;
    std::vector<Position> displacements;

    std::unique_ptr<BoundaryConditionsBase<LatticeType>> boundary_conditions;
    
    double energy;
    
    std::vector<std::shared_ptr<Backbone<LatticeType>>> backbone_list;
    int total_backbones_length;
    int num_backbones;
    //~ std::vector<int> backbones_cumulative_lengths;
	
	InterMatrix matrix;
	
	
	//---------------------------- Methods -----------------------------
	
	int checkSite2D(const int x, const int y) const;
	int checkSite3D(const int x, const int y, const int z) const;
	
	bool isSiteOccupied2D(const int x, const int y) const;
	bool isSiteOccupied3D(const int x, const int y, const int z) const;
	
	bool isSiteOccupiedIgnoringBackbonei2D(const int x, const int y, const size_t ignored_b) const;
	bool isSiteOccupiedIgnoringBackbonei3D(const int x, const int y, const int z, const size_t ignored_b) const;
	
	void setOccupancy2D(const int x, const int y, const int global_id);
	void setOccupancy3D(const int x, const int y, const int z, const int global_id);
	
};





 
#include "configuration.tpp"


#endif // CONFIGURATION_H
