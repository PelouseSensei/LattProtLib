#ifndef LATTICE_H
#define LATTICE_H

#include <memory>
#include <string>
#include <vector>
#include <array>

//----------------------------- Base Class -----------------------------

template <size_t Dim>
class LatticeBase {
public:
    static constexpr size_t dimension = Dim;

    using Position = std::array<int, Dim>;

    virtual std::vector<Position> getDisplacements() const = 0;
    //~ virtual Position applyPeriodic(const Position& pos, int L) const = 0;
    virtual std::string getLatticeType() const = 0;

    virtual ~LatticeBase() = default;
};


//------------------------- Square Lattice -----------------------------

class SquareLattice : public LatticeBase<2> {
public:
    std::vector<Position> getDisplacements() const override {
        return {{ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }};
    }
    
    std::string getLatticeType() const override { return "Square"; }

    //~ Position applyPeriodic(const Position& pos, int L) const override {
        //~ return { (pos[0] + L) % L, (pos[1] + L) % L };
    //~ }
};


//-------------------------- Cubic Lattice -----------------------------

class CubicLattice : public LatticeBase<3> {
public:
    std::vector<Position> getDisplacements() const override {
        return {{ {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1} }};
    }
    
    std::string getLatticeType() const override { return "Cubic"; }

    //~ Position applyPeriodic(const Position& pos, int L) const override {
        //~ return { (pos[0] + L) % L, (pos[1] + L) % L, (pos[2] + L) % L };
    //~ }
};


#endif // LATTICE_H
