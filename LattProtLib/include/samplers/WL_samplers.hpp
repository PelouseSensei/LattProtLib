#ifndef WL_SAMPLERS_H
#define WL_SAMPLERS_H

#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <highfive/H5File.hpp>
#include <string>
#include "../lattice/lattice.hpp"
#include "../moves/moves.hpp"
#include "../system/configuration.hpp"
#include "../samplers/mcmc_samplers.hpp"
#include "../schedulers/task_scheduler.hpp"

struct WangLandauAcceptance {
    template<typename RNG>
    bool operator()(double g_ratio, double /*unsused beta param asked by the base sampler*/, RNG& rng) const {
        if (g_ratio > 1) return true;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(rng) < g_ratio;
    }
};


template<typename LatticeType, typename RNG = std::mt19937>
class WangLandauSampler
    : public BaseSampler<LatticeType, WangLandauAcceptance, RNG>
{
public:
    using ConfigType = Configuration<LatticeType>;
    using Base = BaseSampler<LatticeType, WangLandauAcceptance, RNG>;
    using EnergyFn = std::function<int(double)>; // map energy -> bin index

    // Constructor using uniform linear binning
    WangLandauSampler(int seed,
                      double initial_ln_f = 1.0,    // ln(f) initial
                      double emin = -100.0,
                      double emax = 100.0,
                      int n_bins = 1001,
                      TaskScheduler<ConfigType>* scheduler = nullptr,
                      RNG rng = RNG())
        : Base(seed, /*beta not used in WL but keep base ctor param*/ 1.0),
          ln_f_(initial_ln_f),
          emin_(emin),
          emax_(emax),
          n_bins_(n_bins),
          scheduler_(scheduler),
          flatness_threshold_(0.8),
          min_ln_f_(1e-8),
          check_interval_(1000)
    {
        if (n_bins_ <= 0) throw std::invalid_argument("n_bins must be > 0");
        bin_width_ = (emax_ - emin_) / static_cast<double>(n_bins_);
        if (bin_width_ <= 0.0) throw std::invalid_argument("Invalid energy range");
        
        bins_inf_bound_ = 0;
        bins_sup_bound_ = n_bins_-1;

        logg_.assign(n_bins_, 0.0);
        hist_.assign(n_bins_, 0);

        // default energy->bin mapping uses uniform linear grid
        energy_to_bin_ = [this](double E) -> int {
            int b = static_cast<int>(std::floor((E - emin_) / bin_width_));
            return b;
        };

        // set RNG inside base (note Base has rng_ member, but we constructed Base with seed; override if provided)
        this->rng_ = rng;
    }

    // Alternate ctor with custom energy -> bin mapping
    WangLandauSampler(int seed,
                      double initial_ln_f,
                      double emin,
                      double emax,
                      int n_bins,
                      EnergyFn energy_to_bin,
                      TaskScheduler<ConfigType>* scheduler = nullptr,
                      RNG rng = RNG())
        : Base(seed, 1.0),
          ln_f_(initial_ln_f),
          emin_(emin),
          emax_(emax),
          n_bins_(n_bins),
          scheduler_(scheduler),
          flatness_threshold_(0.8),
          min_ln_f_(1e-8),
          check_interval_(1000)
    {
        if (n_bins_ <= 0) throw std::invalid_argument("n_bins must be > 0");
        bin_width_ = (emax_ - emin_) / static_cast<double>(n_bins_);
        std::cout << "bin width: " << bin_width_ << std::endl;
        if (bin_width_ <= 0.0) throw std::invalid_argument("Invalid energy range");
        
        bins_inf_bound_ = 0;
        bins_sup_bound_ = n_bins_-1;
        
        logg_.assign(n_bins_, 0.0);
        hist_.assign(n_bins_, 0);
        
        energy_to_bin_ = std::move(energy_to_bin);
        this->rng_ = rng;
    }

    // ----- Configuration helpers -----

    // Must be called before stepping to initialize the visited bin of the current configuration
    void initialize(ConfigType& config, bool mark_visit = true) {
        double E = config.getEnergy();
        int bin = energyToBinChecked(E);
        if (mark_visit) {
            logg_[bin] += ln_f_;
            hist_[bin] += 1;
        }
    }

    // Set scheduler after construction if you like
    void setScheduler(TaskScheduler<ConfigType>* scheduler) {
        scheduler_ = scheduler;
    }

    // customize flatness criterion (default 0.8)
    void setFlatnessThreshold(double t) {
        if (t <= 0.0 || t > 1.0) throw std::invalid_argument("flatness threshold must be in (0,1]");
        flatness_threshold_ = t;
    }

    void setCheckInterval(int interval) {
        if (interval <= 0) throw std::invalid_argument("interval must be > 0");
        check_interval_ = interval;
    }

    void setMinLnF(double min_ln_f) {
        min_ln_f_ = min_ln_f;
    }
    
    void setBinsInfBound(const int new_inf_bound) {
		if (new_inf_bound >= 0) {
			bins_inf_bound_ = new_inf_bound;
        } else {
			throw std::out_of_range("setBinsInfBound: new inf bin bound < 0");
		}
	}
    
    void setBinsSupBound(const int new_sup_bound) {
		if (new_sup_bound <= n_bins_) {
			bins_sup_bound_ = new_sup_bound;
        } else {
			throw std::out_of_range("setBinsInfBound: new inf bin bound > number of bins");
		}
	}
	
	
	// Accessors
	int getBinsInfBound() {
		return bins_inf_bound_;
	}
	
	int getBinsSupBound() {
		return bins_sup_bound_;
	}
    
    const std::vector<double>& getLogDOS() const { return logg_; } // ln g
    const std::vector<int>& getHistogram() const { return hist_; }
    double getLnF() const { return ln_f_; }
    int getNBins() const { return n_bins_; }
    
    int getFirstOccupiedLogGBin() {
		for (size_t i = 0; i<logg_.size(); i++) {
			if (logg_[i] != 0) { return i; }
		}
		return logg_.size() - 1;
	}
	
	int getLastOccupiedLogGBin() {
		for (int i = static_cast<int>(logg_.size()) - 1; i >= 0; --i) {
			if (logg_[i] != 0) { return i; }
		}
		return 0;
	}


    void sample(ConfigType& config) override {
		// choose move via base move_selector_
		int move_index = this->move_selector_(this->rng_);
		auto& move = *this->moves_[move_index];

		move.generateMoveArgs(config);
		move.checkMove(config);

		++this->attempted_moves_[move_index];
		
		double E = config.getEnergy();
		int binE = energyToBin(E);  // always get current bin
		//~ if (binE < bins_inf_bound_ || binE > bins_sup_bound_) { 
			//~ return;
		//~ }
		

		int visited_bin = binE;  // default: current energy bin

		if (move.getStatus() == MoveStatus::Possible) {
			++this->possible_moves_[move_index];

			double dE = move.getdE();
			double Enew = E + dE;

			int binEnew = energyToBin(Enew);

			if (binEnew >= bins_inf_bound_ && binEnew <= bins_sup_bound_) {
				double ln_gE = logg_[binE];
				double ln_gEnew = logg_[binEnew];

				// acceptance probability = min(1, exp(ln_gE - ln_gEnew))
				double g_ratio_exp = std::exp(ln_gE - ln_gEnew);

				// Use Wang–Landau acceptance policy
				if (this->acceptance_policy_(g_ratio_exp, 1.0, this->rng_)) {
					move.apply(config);
					++this->accepted_moves_[move_index];
					//~ visited_bin = energyToBinChecked(config.getEnergy());
					visited_bin = binEnew;
				}
			}
			// Update ln g(E) as long as the move is possible, even if the move was refused because not in the energy bounds
			logg_[visited_bin] += ln_f_;
			hist_[visited_bin] += 1;
		}
	}


    // Check whether histogram meets the flatness criterion
    // (returns true if flat enough)
   bool checkFlatness(int min_count = 1) const {
		// Count only bins that have been visited at least min_count times
		int considered_bins = 0;
		long long sum = 0;
		int minv = std::numeric_limits<int>::max();

		for (int h : hist_) {
			if (h >= min_count) {
				++considered_bins;
				sum += h;
				if (h < minv) minv = h;
			}
		}

		// If no bins meet the min_count threshold → not flat
		if (considered_bins == 0) return false;

		double avg = static_cast<double>(sum) / static_cast<double>(considered_bins);
		if (avg <= 0.0) return false;

		double ratio = static_cast<double>(minv) / avg;
		return ratio >= flatness_threshold_;
	}
	
	// Check whether histogram meets the flatness criterion
	// threshold: minimum value for all non-empty bins
	bool checkFlatnessThreshold(int threshold = 1) const {
		int considered_bins = 0;
		long long sum = 0;
		int minv = std::numeric_limits<int>::max();

		for (int h : hist_) {
			if (h > 0) {  // only non-empty bins are considered
				++considered_bins;
				sum += h;
				if (h < minv) minv = h;
			}
		}

		// If no bins are non-empty → not flat
		if (considered_bins == 0) return false;

		// Original flatness criterion: min / avg >= flatness_threshold_
		double avg = static_cast<double>(sum) / static_cast<double>(considered_bins);
		if (avg <= 0.0) return false;
		double ratio = static_cast<double>(minv) / avg;

		// Also enforce that all considered bins are above threshold
		return (ratio >= flatness_threshold_) && (minv >= threshold);
	}
	
	// Checks flatness on masked bins
	// Eg Checks flatness only on occupied bins found with getOccupiedEnergyBinIdx()
	bool checkFlatnessMask(const std::vector<int>& maskedBins) {
		int considered_bins = 0;
		long long sum = 0;
		int minv = std::numeric_limits<int>::max();
		
		for (int bin_idx : maskedBins) {
			++considered_bins;
			long h = hist_[bin_idx];
			sum += h;
			if (h < minv) { minv = h; }
		}
		
		// If no bins meet the min_count threshold → not flat
		if (considered_bins == 0) return false;

		double avg = static_cast<double>(sum) / static_cast<double>(considered_bins);
		if (avg <= 0.0) return false;

		double ratio = static_cast<double>(minv) / avg;
		return ratio >= flatness_threshold_;
	}
	
	// Returns the indices of energy bins occupied by at least threshold counts
	std::vector<int> getOccupiedEnergyBinIdx(int threshold = 0) {
		
		std::vector<int> occupied_bins; 
		
		for (size_t i = 0; i<logg_.size(); ++i) {
			if (logg_[i] > threshold) { occupied_bins.push_back(i) ; }
		}
		
		return occupied_bins;
	}
	
	

	
	double getMinLnF() {
		return min_ln_f_;
	}

    bool updateModificationFactor() {
        if (ln_f_ <= min_ln_f_) return false;
        ln_f_ /= 2.0;
        std::fill(hist_.begin(), hist_.end(), 0);
        return true;
    }
    
    void setLnF(const double new_ln_f) {
		ln_f_ = new_ln_f;
	}

    void step(ConfigType& config, int global_step_counter) {
        sample(config);

        if (scheduler_) scheduler_->run(config, global_step_counter);

        ++steps_since_check_;
        if (steps_since_check_ >= check_interval_) {
            steps_since_check_ = 0;
            if (checkFlatness()) {
                bool updated = updateModificationFactor();
                // user may want to be notified - run tasks once more with special step id -1 or provide callback
                // For now we'll call scheduler again with the same step so they can check ln_f
                if (scheduler_) scheduler_->run(config, global_step_counter);
                // Optionally print status:
                std::cout << "[WangLandau] Histogram flat at step " << global_step_counter
                          << ", ln f reduced to " << ln_f_ << "\n";
                (void)updated;
            }
        }
    }

    // Run N steps (convenience). It's a wrapper invoking step() N times.
    void run(ConfigType& config, int nSteps, int step_offset = 0) {
        for (int i = 0; i < nSteps; ++i) {
            step(config, step_offset + i);
            //~ if (ln_f_ <= min_ln_f_) {
				//~ std::cout << "DOS converged... stopping run" << std::endl;
				//~ this->printStatistics();
				//~ return;
			//~ }
        }
    }

    // Reset histogram and/or ln g
    void resetHistogram() { std::fill(hist_.begin(), hist_.end(), 0); }
    void resetLogDOS(double value = 0.0) { std::fill(logg_.begin(), logg_.end(), value); }

    // energy -> bin mapping with range check (returns -1 if out of range)
    int energyToBin(double E) const {
        if (energy_to_bin_) {
            int b = energy_to_bin_(E);
            return b;
        } else {
            // shouldn't happen since ctor sets energy_to_bin_
            int b = static_cast<int>(std::floor((E - emin_) / bin_width_));
            return b;
        }
    }

    // checked variant that throws if out of range
    int energyToBinChecked(double E) const {
        int b = energyToBin(E);
        if (b < 0 || b >= n_bins_) {
            throw std::out_of_range("energyToBinChecked: energy out of bin range");
        }
        return b;
    }
    
    void printHistogram(std::ostream& os = std::cout,
						bool skip_empty = false,
						bool show_bin_centers = true) const
	{
		os << "# Histogram H(E)\n";
		os << "# Columns: bin_index";
		if (show_bin_centers) os << "  E_center";
		os << "  H(E)\n";

		for (int i = 0; i < n_bins_; ++i) {
			int h = hist_[i];
			if (skip_empty && h == 0) continue;

			if (show_bin_centers) {
				double E_center = emin_ + (i + 0.5) * bin_width_;
				os << i << "\t" << E_center << "\t" << h << "\n";
			} else {
				os << i << "\t" << h << "\n";
			}
		}
		os << std::flush;
	}

	void printLogDOS(std::ostream& os = std::cout,
                 bool skip_unvisited = false,
                 bool show_bin_centers = true,
                 std::function<double(int)> bin_to_energy = nullptr) const
{
    os << "# Log density of states ln g(E)\n";
    os << "# Columns: bin_index";
    if (show_bin_centers) os << "  E_center";
    os << "  ln_g(E)\n";

    for (int i = 0; i < n_bins_; ++i) {
        double ln_g = logg_[i];
        if (skip_unvisited && hist_[i] == 0) continue;

        if (show_bin_centers) {
            double E_center;

            if (bin_to_energy) {
                // Use custom nonlinear energy mapping
                double E_low = bin_to_energy(i);
                double E_high;
                if (i < n_bins_ - 1)
                    E_high = bin_to_energy(i + 1);
                else
                    // Approximate last bin center
                    E_high = E_low + (E_low - bin_to_energy(i - 1));

                E_center = 0.5 * (E_low + E_high);
            } else {
                // Fall back to linear binning
                E_center = emin_ + (i + 0.5) * bin_width_;
            }

            os << i << "\t" << E_center << "\t" << ln_g << "\n";
        } else {
            os << i << "\t" << ln_g << "\n";
        }
    }

    os << std::flush;
}

	
	// Save Wang–Landau data to an HDF5 file
	void saveToHDF5(const std::string& filename, const std::string& group_name = "/WangLandau") const {
		using namespace HighFive;

		// Create (or open existing) file in read/write mode
		File file(filename, File::ReadWrite | File::Create | File::Truncate);

		// Create a group for this dataset (e.g. /WangLandau)
		Group group = file.createGroup(group_name);

		// Write datasets
		group.createDataSet("logg", logg_);
		group.createDataSet("histogram", hist_);

		// Write attributes (metadata)
		group.createAttribute("emin", emin_);
		group.createAttribute("emax", emax_);
		group.createAttribute("n_bins", n_bins_);
		group.createAttribute("bin_width", bin_width_);
		group.createAttribute("ln_f", ln_f_);

		//~ std::cout << "[WangLandauSampler] Saved HDF5 file: " << filename
				  //~ << " (group: " << group_name << ")\n";
	}
	
	void appendToHDF5(const std::string& filename,
					  int mc_step,
					  bool append = true) const {
		using namespace HighFive;

		// Open or create the file
		File file(filename,
				  append ? (File::ReadWrite | File::Create)
						 : (File::ReadWrite | File::Create | File::Truncate));

		// --- Ensure /metadata group exists ---
		Group meta;
		if (!file.exist("/metadata")) {
			meta = file.createGroup("/metadata");
			meta.createDataSet("emin", emin_);
			meta.createDataSet("emax", emax_);
			meta.createDataSet("bin_width", bin_width_);
			meta.createDataSet("n_bins", n_bins_);
		} else {
			meta = file.getGroup("/metadata");
		}

		// --- Ensure /WL_data exists ---
		Group wl_group;
		if (!file.exist("/WL_data")) {
			wl_group = file.createGroup("/WL_data");
		} else {
			wl_group = file.getGroup("/WL_data");
		}

		// --- Create a group for this snapshot ---
		std::string group_name = "step_" + std::to_string(mc_step);
		if (wl_group.exist(group_name)) {
			wl_group.unlink(group_name); // overwrite existing
		}
		Group g = wl_group.createGroup(group_name);

		// --- Write datasets ---
		g.createDataSet("logg", logg_);
		g.createDataSet("histogram", hist_);

		// --- Write attributes ---
		g.createAttribute("ln_f", ln_f_);
		g.createAttribute("MC_step", mc_step);

		// Optional logging (quiet by default)
		// std::cout << "[WL] Appended step " << mc_step
		//           << " to " << filename << " (/WL_data/" << group_name << ")\n";
	}
	
	// --- Dynamically modify energy bounds (emin, emax) ---
    void setEnergyRange(double new_emin, double new_emax) {
        if (new_emax <= new_emin)
            throw std::invalid_argument("setEnergyRange: emax must be > emin");
        emin_ = new_emin;
        emax_ = new_emax;
        bin_width_ = (emax_ - emin_) / static_cast<double>(n_bins_);
        std::cout << "[WangLandauSampler] Energy range updated: ["
                  << emin_ << ", " << emax_ << "], bin width = "
                  << bin_width_ << "\n";
    }

    // --- Modify the number of bins and reinitialize arrays ---
    void setNBins(int new_n_bins, bool keep_logg = false) {
        if (new_n_bins <= 0)
            throw std::invalid_argument("setNBins: number of bins must be > 0");

        n_bins_ = new_n_bins;
        bin_width_ = (emax_ - emin_) / static_cast<double>(n_bins_);

        std::vector<double> new_logg(n_bins_, 0.0);
        std::vector<int>    new_hist(n_bins_, 0);

        if (keep_logg && !logg_.empty()) {
            // crude resampling: copy overlapping part
            int copy_bins = std::min<int>(n_bins_, logg_.size());
            std::copy_n(logg_.begin(), copy_bins, new_logg.begin());
        }

        logg_.swap(new_logg);
        hist_.swap(new_hist);

        std::cout << "[WangLandauSampler] Number of bins set to "
                  << n_bins_ << ", bin width = " << bin_width_ << "\n";
    }

    // --- Initialize current logg_ from an external vector ---
    void initializeLogDOS(const std::vector<double>& new_logg) {
        if (static_cast<int>(new_logg.size()) != n_bins_)
            throw std::invalid_argument(
                "initializeLogDOS: size of new_logg does not match n_bins");
        logg_ = new_logg;
        std::cout << "[WangLandauSampler] logg_ initialized from external data\n";
    }



private:
    double ln_f_;                    // ln(f)
    std::vector<double> logg_;       // ln g(E) per bin
    std::vector<int> hist_;          // histogram H(E) per bin

    // energy binning
    double emin_, emax_;
    int n_bins_;
    double bin_width_;
    EnergyFn energy_to_bin_;
    int bins_inf_bound_;
    int bins_sup_bound_;

    // scheduler to run tasks after each step
    TaskScheduler<ConfigType>* scheduler_;

    // flatness / termination parameters
    double flatness_threshold_;
    double min_ln_f_;
    int check_interval_;
    int steps_since_check_ = 0;
};




//---------------------- Square Wang-Landau Pull -----------------------
template<typename RNG = std::mt19937>
class SquareWangLandauSampler_pull : public WangLandauSampler<SquareLattice> {
public:
    SquareWangLandauSampler_pull(int seed,
                      double initial_ln_f = 1.0,    // ln(f) initial
                      double emin = -100.0,
                      double emax = 100.0,
                      int n_bins = 1001,
                      TaskScheduler<ConfigType>* scheduler = nullptr,
                      RNG rng = RNG())
        : WangLandauSampler<SquareLattice>(seed, initial_ln_f, emin, emax, n_bins, scheduler, rng)
    {
        this->setupMoves();
    }
    // Alternate ctor with custom energy -> bin mapping
    SquareWangLandauSampler_pull(int seed,
                      double initial_ln_f,
                      double emin,
                      double emax,
                      int n_bins,
                      EnergyFn energy_to_bin,
                      TaskScheduler<ConfigType>* scheduler = nullptr,
                      RNG rng = RNG())
		: WangLandauSampler<SquareLattice>(seed, initial_ln_f, emin, emax, n_bins, energy_to_bin, scheduler, rng)
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


#endif // WL_SAMPLERS_H
