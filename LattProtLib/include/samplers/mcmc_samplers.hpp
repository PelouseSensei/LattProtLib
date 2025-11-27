#ifndef MCMC_SAMPLERS_H
#define MCMC_SAMPLERS_H

#include "../lattice/lattice.hpp"
#include "../moves/moves.hpp"
#include "../system/configuration.hpp"
#include <array>

//------------------------ Acceptance Policies -------------------------

struct MetropolisAcceptance {
    template<typename RNG>
    bool operator()(double dE, double beta, RNG& rng) const {
        if (dE <= 0) return true;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(rng) < std::exp(-beta * dE);
    }
};






//---------------------------- Base sampler ----------------------------

template<typename LatticeType, typename AcceptancePolicy, typename RNG = std::mt19937>
class BaseSampler {
public:
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;

    BaseSampler(int seed, double beta)
        : seed_(seed),
          rng_(seed),
          beta_(beta),
          acceptance_policy_()
    {}

    virtual ~BaseSampler() = default;

    virtual void setupMoves() = 0;

    virtual void sample(Configuration<LatticeType>& config) {
		int move_index = move_selector_(rng_);
		auto& move = *moves_[move_index];

		move.generateMoveArgs(config);
		move.checkMove(config);

		++attempted_moves_[move_index];

		if (move.getStatus() == MoveStatus::Possible) {
			++possible_moves_[move_index];

			if (acceptance_policy_(move.getdE(), beta_, rng_)) {
				move.apply(config);
				++accepted_moves_[move_index];
			}
		}
	}

    void run(Configuration<LatticeType>& config, int nsteps) {
        for (int i = 0; i < nsteps; ++i) {
            sample(config);
        }
    }

	void printStatistics() const {
		for (size_t i = 0; i < moves_.size(); ++i) {
			double hit_rate = attempted_moves_[i] ?
							  static_cast<double>(possible_moves_[i]) / attempted_moves_[i] : 0.0;

			double acceptance_rate = attempted_moves_[i] ?
							  static_cast<double>(accepted_moves_[i]) / attempted_moves_[i] : 0.0;

			double global_acceptance = possible_moves_[i] ?
							  static_cast<double>(accepted_moves_[i]) / possible_moves_[i] : 0.0;

			std::cout << "Move " << i << ": "
					  << "Attempts = " << attempted_moves_[i]
					  << ", Possible = " << possible_moves_[i]
					  << ", Accepted = " << accepted_moves_[i]
					  << ", Hit Rate = " << hit_rate
					  << ", Acceptance Rate = " << acceptance_rate
					  << ", Global Acceptance = " << global_acceptance
					  << "\n";
		}
	}
    
    void resetStatistics() {
        std::fill(attempted_moves_.begin(), attempted_moves_.end(), 0);
        std::fill(accepted_moves_.begin(), accepted_moves_.end(), 0);
        std::fill(possible_moves_.begin(), possible_moves_.end(), 0);
    }
    
    std::vector<double> getAcceptanceRates() const {
        std::vector<double> rates(moves_.size(), 0.0);
        for (size_t i = 0; i < moves_.size(); ++i) {
            if (attempted_moves_[i] > 0) {
                rates[i] = static_cast<double>(accepted_moves_[i]) / attempted_moves_[i];
            }
        }
        return rates;
    }
    
    std::vector<double> getHitRates() const {
		std::vector<double> rates(moves_.size(), 0.0);
		for (size_t i = 0; i < moves_.size(); ++i) {
			if (attempted_moves_[i] > 0) {
				rates[i] = static_cast<double>(possible_moves_[i]) / attempted_moves_[i];
			}
		}
		return rates;
	}
	
	std::vector<double> getGlobalAcceptanceRates() const {
		std::vector<double> rates(moves_.size(), 0.0);
		for (size_t i = 0; i < moves_.size(); ++i) {
			if (possible_moves_[i] > 0) {
				rates[i] = static_cast<double>(accepted_moves_[i]) / possible_moves_[i];
			}
		}
		return rates;
	}
    
    void printWeights() {
		std::cout << "Move weigths : " << std::endl;
		for (size_t i = 0; i<move_weights_.size(); i++) {
			std::cout << "Move " << i << ": " << move_weights_[i] << " - ";
		}
		std::cout << std::endl;
	}
	
	void setWeightMovei(const int i, const double new_weight) {
		move_weights_[i] = new_weight;
		updateSelector();
	}
	
	template<typename ConcreteMove>
	ConcreteMove* getMovei(const int i) {
		return dynamic_cast<ConcreteMove*>(moves_[i].get());
	}
	
	void setBeta(const double new_beta) { beta_ = new_beta; };

protected:
    void addMove(std::unique_ptr<Move<LatticeType>> move, double weight) {
        moves_.push_back(std::move(move));
        move_weights_.push_back(weight);
        updateSelector();
        attempted_moves_.push_back(0);
        accepted_moves_.push_back(0);
        possible_moves_.push_back(0);
    }

    void updateSelector() {
        move_selector_ = std::discrete_distribution<int>(move_weights_.begin(), move_weights_.end());
    }

protected:
    int seed_;
    RNG rng_;
    double beta_;

    AcceptancePolicy acceptance_policy_;

    std::vector<std::unique_ptr<Move<LatticeType>>> moves_;
    std::vector<double> move_weights_;
    std::discrete_distribution<int> move_selector_;

    std::vector<int> attempted_moves_;
    std::vector<int> accepted_moves_;
    std::vector<int> possible_moves_;

};




//---------------------- Square Metropolis-Hasting VSHD -----------------------

class SquareMetropolisSampler : public BaseSampler<SquareLattice, MetropolisAcceptance> {
public:
    SquareMetropolisSampler(int seed, double beta)
        : BaseSampler<SquareLattice, MetropolisAcceptance>(seed, beta) 
    {
        this->setupMoves();
    }

    void setupMoves() override {
        this->addMove(std::make_unique<CornerSquare>(this->seed_), 1.0);
        this->addMove(std::make_unique<TailSquare>(this->seed_),   1.0);
        this->addMove(std::make_unique<CrankshaftSquare>(this->seed_),   1.0);
        this->addMove(std::make_unique<TranslationSquare>(this->seed_),   1.0);
        this->addMove(std::make_unique<RotationSquare>(this->seed_), 1.0);
        this->addMove(std::make_unique<SymmetrySquare>(this->seed_), 1.0);
    }
};


//---------------------- Square Metropolis-Hasting Pull -----------------------

class SquareMetropolisSampler_pull : public BaseSampler<SquareLattice, MetropolisAcceptance> {
public:
    SquareMetropolisSampler_pull(int seed, double beta)
        : BaseSampler<SquareLattice, MetropolisAcceptance>(seed, beta) 
    {
        this->setupMoves();
    }

    void setupMoves() override {
        this->addMove(std::make_unique<PullBackwardSquare>(this->seed_), 1.0);
        this->addMove(std::make_unique<PullForwardSquare>(this->seed_), 1.0);
        this->addMove(std::make_unique<PullEndSquare>(this->seed_), 1.0);
        this->addMove(std::make_unique<TranslationSquare>(this->seed_),   1.0);
        this->addMove(std::make_unique<RotationSquare>(this->seed_), 1.0);
        this->addMove(std::make_unique<SymmetrySquare>(this->seed_), 1.0);
    }
};


//---------------------- Cubic Metropolis-Hasting Pull -----------------------

class CubicMetropolisSampler_pull : public BaseSampler<CubicLattice, MetropolisAcceptance> {
public:
    CubicMetropolisSampler_pull(int seed, double beta)
        : BaseSampler<CubicLattice, MetropolisAcceptance>(seed, beta) 
    {
        this->setupMoves();
    }

    void setupMoves() override {
        this->addMove(std::make_unique<PullBackwardCubic>(this->seed_), 1.0);
        this->addMove(std::make_unique<PullForwardCubic>(this->seed_), 1.0);
        this->addMove(std::make_unique<PullEndCubic>(this->seed_), 1.0);
        this->addMove(std::make_unique<TranslationCubic>(this->seed_),   1.0);
        this->addMove(std::make_unique<RotationCubic>(this->seed_), 1.0);
        this->addMove(std::make_unique<SymmetryCubic>(this->seed_), 1.0);
    }
};


//-------------------- Arbitrary lattice Sequence Sampler ----------------------

template<typename LatticeType>
class MHSequenceSampler : public BaseSampler<LatticeType, MetropolisAcceptance> {
public:
	MHSequenceSampler(int seed, double beta) : BaseSampler<LatticeType, MetropolisAcceptance>(seed, beta) {
		this->setupMoves();
	}
	
	void setupMoves() override {
		this->addMove(std::make_unique<OnePointMutation<LatticeType>>(this->seed_), 1.0);
		this->addMove(std::make_unique<TwoPointSwitch<LatticeType>>(this->seed_), 1.0);
	}
};

#endif // MCMC_SAMPLERS_H
