#ifndef BASE_MOVE_H
#define BASE_MOVE_H

#include <memory>
#include <any>
#include <random>
#include "../system/configuration.hpp"
#include "../lattice/lattice.hpp"
#include "custom_structs.hpp"

/**-------------------------- Base Move Class -----------------------**/

template<typename LatticeType>
class Move {
public:

	using Position = typename LatticeType::Position;
	static constexpr size_t Dim = LatticeType::dimension;

	Move(int seed) : rng(seed) {}
    virtual ~Move() {}
    
    virtual bool apply(Configuration<LatticeType>& config) = 0;
    
    virtual void generateMoveArgs(Configuration<LatticeType>& config) = 0;
    
    virtual void checkMove(Configuration<LatticeType>& config) = 0;
    
    virtual double getdE() = 0;
    
    virtual MoveStatus getStatus() = 0;

protected:
	std::mt19937 rng;

};

#endif // BASE_MOVE_H

