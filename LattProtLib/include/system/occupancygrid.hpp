#ifndef OCCUPANCYGRID_H
#define OCCUPANCYGRID_H

#include <stdexcept>
#include <iostream>
//~ #include <type_traits>
#include "../lattice/lattice.hpp"



template <typename LatticeType, typename T = int>
class OccupancyGrid {
public:
    static constexpr size_t Dim = LatticeType::dimension;
    static_assert(Dim == 2 || Dim == 3, "Only 2D or 3D supported.");

    OccupancyGrid(int L)
        : L(L), L2(Dim == 3 ? L * L : 0),
          data(static_cast<size_t>(std::pow(L, Dim)), -1) {}

    // 2D accessor
    T& operator()(int x, int y) requires (Dim == 2) {
        return data[index(x, y)];
    }

    const T& operator()(int x, int y) const requires (Dim == 2) {
        return data[index(x, y)];
    }

    // 3D accessor
    T& operator()(int x, int y, int z) requires (Dim == 3) {
        return data[index(x, y, z)];
    }

    const T& operator()(int x, int y, int z) const requires (Dim == 3) {
        return data[index(x, y, z)];
    }
	
	// Only defined if dimension == 2
    void print2D() const requires (Dim == 2) {
        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; ++x) {
                std::cout << (*this)(x, y) << " ";
            }
            std::cout << '\n';
        }
    }
    
    // Only defined if dimension == 3
	void print3D() const requires (Dim == 3) {
		for (int z = 0; z < L; ++z) {
			std::cout << "Slice z = " << z << ":\n";
			for (int y = 0; y < L; ++y) {
				for (int x = 0; x < L; ++x) {
					std::cout << (*this)(x, y, z) << " ";
				}
				std::cout << '\n';
			}
			std::cout << '\n';
		}
	}


private:
    int L;
    int L2;  // L * L for 3D
    std::vector<T> data;

    // 2D index
    inline size_t index(int x, int y) const {
        return static_cast<size_t>(x * L + y);
    }

    // 3D index
    inline size_t index(int x, int y, int z) const {
        return static_cast<size_t>(x * L2 + y * L + z);
    }
};


#endif //OCCUPANCYGRID_H
