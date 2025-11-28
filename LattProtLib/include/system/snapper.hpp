#ifndef SNAPPER_H
#define SNAPPER_H

#include <highfive/H5File.hpp>
#include <functional>
#include <array>
#include "configuration.hpp"

//---------------------- System Snapshot structs  ----------------------

template<typename LatticeType>
struct energySnapshot {
    double energy;
    
    void takeSnap(Configuration<LatticeType>& config) { 
        energy = config.getEnergy();
    }

    void print(std::ostream& os = std::cout) const {
        os << energy << "\n";
    }

    void writeToFile(std::ofstream& ofs) const {
        ofs.write(reinterpret_cast<const char*>(&energy), sizeof(double));
    }
};



template<typename LatticeType>
struct posSnapshot {
    
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;
    
    size_t nParticles;
    int coords[csts::MAX_TOTAL_MONOMERS * Dim];
    
    void takeSnap(Configuration<LatticeType>& config) {
        nParticles = config.getTotalLength();
        auto backbones = config.getBackbones();

        int tot_i = 0;

        for (auto backbone : backbones) {
            auto b_length = backbone->getLength();
            for (int id = 0; id < b_length; id++) {
                Position pos = backbone->getPosition(id);

                if constexpr (Dim == 2) {
                    coords[tot_i++] = pos[0];
                    coords[tot_i++] = pos[1];
                }
                if constexpr (Dim == 3) {
                    coords[tot_i++] = pos[0];
                    coords[tot_i++] = pos[1];
                    coords[tot_i++] = pos[2];
                }
            }
        }
    }

    void print(std::ostream& os = std::cout) const {
        os << "Particles: " << nParticles << "\n";
        for (size_t i = 0; i < nParticles; i++) {
            os << "(";
            for (size_t d = 0; d < Dim; d++) {
                os << coords[i*Dim + d];
                if (d + 1 < Dim) os << ",";
            }
            os << ")\n";
        }
    }

    void writeToFile(std::ofstream& ofs) const {
        ofs.write(reinterpret_cast<const char*>(&nParticles), sizeof(size_t));
        ofs.write(reinterpret_cast<const char*>(coords), nParticles * Dim * sizeof(int));
    }
    
    
};

template<typename LatticeType>
struct contactSnapshot {
    size_t nBackbones{};
    std::vector<uint16_t> contacts; // flattened storage

    void takeSnap(const Configuration<LatticeType>& config) {
        nBackbones = config.getNumBackbones();
        contacts.resize(nBackbones * (nBackbones + 1) / 2); // only store needed entries

        const auto& cm = config.contact_matrix;
        for (size_t i = 0; i < nBackbones; i++) {
            for (size_t j = i; j < nBackbones; j++) {
                contacts[(i*(2*nBackbones - i + 1))/2 + (j-i)] = cm.get(i, j);
            }
        }
    }

    void print(std::ostream& os = std::cout) const {
        os << "Contact matrix (" << nBackbones << " backbones):\n";
        for (size_t i = 0; i < nBackbones; i++) {
            for (size_t j = 0; j < nBackbones; j++) {
                if (j < i) {
                    os << contacts[(i*(2*nBackbones - i + 1))/2 + (j-i)] << " ";
                } else {
                    os << contacts[(i*(2*nBackbones - i + 1))/2 + (j-i)] << " ";
                }
            }
            os << "\n";
        }
    }

    void writeToFile(std::ofstream& ofs) const {
        ofs.write(reinterpret_cast<const char*>(&nBackbones), sizeof(size_t));
        ofs.write(reinterpret_cast<const char*>(contacts.data()), contacts.size() * sizeof(uint16_t));
    }
};


template<typename LatticeType>
struct seqSnapshot {
    size_t nBackbones{};
    std::vector<std::vector<int>> sequences; // store each backbone's sequence

    void takeSnap(Configuration<LatticeType>& config) {
        auto backbones = config.getBackbones();
        nBackbones = backbones.size();
        sequences.clear();
        sequences.reserve(nBackbones);

        for (const auto& backbone : backbones) {
            size_t length = backbone->getLength();
            std::vector<int> seq(length);
            for (size_t i = 0; i < length; ++i) {
                seq[i] = backbone->getType(i); 
            }
            sequences.push_back(std::move(seq));
        }
    }

    void print(std::ostream& os = std::cout) const {
        os << "Sequences (" << nBackbones << " backbones):\n";
        for (size_t b = 0; b < sequences.size(); ++b) {
            os << "Backbone " << b << ": ";
            for (size_t i = 0; i < sequences[b].size(); ++i) {
                os << sequences[b][i];
                if (i + 1 < sequences[b].size()) os << ",";
            }
            os << "\n";
        }
    }

    void writeToFile(std::ofstream& ofs) const {
        ofs.write(reinterpret_cast<const char*>(&nBackbones), sizeof(size_t));
        for (const auto& seq : sequences) {
            size_t length = seq.size();
            ofs.write(reinterpret_cast<const char*>(&length), sizeof(size_t));
            ofs.write(reinterpret_cast<const char*>(seq.data()), length * sizeof(int));
        }
    }
};


//----------------------------- Classes --------------------------------

//~ template<typename LatticeType, int nSteps>
//~ class EnergySnapContainer {
//~ public:
    //~ using Position = typename LatticeType::Position;
    //~ static constexpr size_t Dim = LatticeType::dimension;
    
    //~ EnergySnapContainer() : head_(0) { container_.resize(nSteps); }
    
    //~ void takeAndSaveSnapshot(Configuration<LatticeType>& config) {
        //~ if (head_ >= nSteps) return; // avoid overflow
        //~ snap_.takeSnap(config);
        //~ container_[head_] = snap_;
        //~ head_++;
    //~ }
    
    //~ void print() const {
        //~ for (size_t i = 0; i < head_; i++) {
            //~ container_[i].print();
        //~ }
    //~ }

    //~ void writeHDF5(const std::string& filename) const {
        //~ using namespace HighFive;
        //~ File file(filename, File::Overwrite); // create new file
        //~ std::vector<double> energies;
        //~ energies.reserve(head_);
        //~ for (size_t i = 0; i < head_; i++) {
            //~ energies.push_back(container_[i].energy);
        //~ }
        //~ DataSet dataset = file.createDataSet<double>("energies",
                                //~ DataSpace::From(energies));
        //~ dataset.write(energies);
    //~ }
    
    //~ void writeTxt(const std::string& filename) const {
        //~ std::ofstream ofs(filename);
        //~ if (!ofs) {
            //~ throw std::runtime_error("Cannot open file " + filename);
        //~ }

        //~ for (size_t i = 0; i < head_; i++) {
            //~ ofs << container_[i].energy << "\n";
        //~ }
        //~ ofs.close();
    //~ }
    
    //~ void resetHead() { head_ = 0; }
    
//~ private:
    //~ energySnapshot<LatticeType> snap_;
    //~ size_t head_;
    //~ std::vector<energySnapshot<LatticeType>> container_;
//~ };

template<typename DataType, int nSteps>
class GenericSnapContainer {
public:
    GenericSnapContainer(const std::string& datasetName = "data")
        : head_(0), datasetName_(datasetName)
    {
        container_.resize(nSteps);
    }

    // Record a single data value
    void takeAndSaveSnapshot(const DataType& value) {
        if (head_ >= nSteps) return;
        container_[head_] = value;
        head_++;
    }

    void print(std::ostream& os = std::cout) const {
        for (size_t i = 0; i < head_; i++) {
            os << container_[i] << "\n";
        }
    }

    void writeTxt(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) throw std::runtime_error("Cannot open file " + filename);
        for (size_t i = 0; i < head_; i++) ofs << container_[i] << "\n";
        ofs.close();
    }

    void writeHDF5(const std::string& filename) const {
		using namespace HighFive;
		File file(filename, File::Overwrite);

		std::vector<DataType> data(container_.begin(), container_.begin() + head_);
		DataSet dataset = file.createDataSet<DataType>(
			datasetName_, DataSpace::From(data));
		dataset.write(data);
	}

    void resetHead() { head_ = 0; }

    size_t size() const { return head_; }

private:
    size_t head_;
    std::vector<DataType> container_;
protected:
    std::string datasetName_;
};




template<typename LatticeType, int nSteps>
class EnergySnapContainer : public GenericSnapContainer<double, nSteps> {
    using Base = GenericSnapContainer<double, nSteps>;
public:
    using Base::Base; // inherit constructor
    
    EnergySnapContainer() {
		this->datasetName_ = "energies";
	}
    
    void takeAndSaveSnapshot(Configuration<LatticeType>& config) {
        double e = config.getEnergy();
        Base::takeAndSaveSnapshot(e);
    }
};




template<typename LatticeType, int nSteps>
class PosSnapContainer {
public:
    using Position = typename LatticeType::Position;
    static constexpr size_t Dim = LatticeType::dimension;
    
    PosSnapContainer() : head_(0) { container_.resize(nSteps); }
    
    void takeAndSaveSnapshot(Configuration<LatticeType>& config) {
        if (head_ >= nSteps) return; // avoid overflow
        snap_.takeSnap(config);
        container_[head_] = snap_;
        head_++;
    }
    
    void print() const {
        for (size_t i = 0; i < head_; i++) {
            container_[i].print();
        }
    }

    void writeHDF5(const std::string& filename) const {
		using namespace HighFive;
		File file(filename, File::Overwrite);

		// Number of columns actually used
		size_t nCols = container_[0].nParticles * Dim;

		// Pre-allocate outer vector
		std::vector<std::vector<int>> all_coords;
		all_coords.reserve(head_);

		// Single reusable buffer
		std::vector<int> row_buffer;
		row_buffer.resize(nCols);

		for (size_t s = 0; s < head_; s++) {
			// Copy only the valid coordinates into the buffer
			std::copy(container_[s].coords,
					  container_[s].coords + nCols,
					  row_buffer.begin());

			// Move a copy of the buffer into the outer vector
			all_coords.push_back(row_buffer);
		}

		std::vector<size_t> dims = {head_, nCols};
		DataSet dataset = file.createDataSet<int>("coords", DataSpace(dims));
		dataset.write(all_coords);
	}
	
	void writeTxt(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) {
            throw std::runtime_error("Cannot open file " + filename);
        }

        for (size_t s = 0; s < head_; s++) {
            // optional: write nParticles first
            ofs << container_[s].nParticles << " ";

            for (size_t i = 0; i < container_[s].nParticles * Dim; i++) {
                ofs << container_[s].coords[i] << " ";
            }
            ofs << "\n";
        }
        ofs.close();
    }
    
    void resetHead() { head_ = 0; }
    
private:
    posSnapshot<LatticeType> snap_;
    size_t head_;
    std::vector<posSnapshot<LatticeType>> container_;
};

template<typename LatticeType, int nSteps>
class ContactSnapContainer {
public:
    ContactSnapContainer() : head_(0) { container_.resize(nSteps); }

    void takeAndSaveSnapshot(const Configuration<LatticeType>& config) {
        if (head_ >= nSteps) return; // avoid overflow
        snap_.takeSnap(config);
		
        container_[head_] = snap_;
        
        head_++;
    }

    void print() const {
        for (size_t i = 0; i < head_; i++) {
            container_[i].print();
        }
    }

    void writeTxt(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) throw std::runtime_error("Cannot open file " + filename);

        for (size_t s = 0; s < head_; s++) {
            const auto& snap = container_[s];
            ofs << snap.nBackbones << " ";
            for (auto val : snap.contacts) {
                ofs << val << " ";
            }
            ofs << "\n";
        }
        ofs.close();
    }

    void writeHDF5(const std::string& filename) const {
        using namespace HighFive;
        File file(filename, File::Overwrite);

        // Precompute max entries among all snapshots for rectangular dataset
        size_t maxEntries = 0;
        for (size_t s = 0; s < head_; s++) {
            maxEntries = std::max(maxEntries, snapLength(container_[s]));
        }

        std::vector<std::vector<uint16_t>> all_contacts(head_, std::vector<uint16_t>(maxEntries, 0));

        for (size_t s = 0; s < head_; s++) {
            const auto& snap = container_[s];
            //~ container_[s].print();
            std::copy(snap.contacts.begin(), snap.contacts.end(), all_contacts[s].begin());
        }

        std::vector<size_t> dims = {head_, maxEntries};
        DataSet dataset = file.createDataSet<uint16_t>("contacts", DataSpace(dims));
        dataset.write(all_contacts);
    }

    void resetHead() { head_ = 0; }

private:
    contactSnapshot<LatticeType> snap_;
    size_t head_;
    std::vector<contactSnapshot<LatticeType>> container_;

    static size_t snapLength(const contactSnapshot<LatticeType>& snap) {
        return snap.contacts.size();
    }
};


template<typename LatticeType, int nSteps>
class SeqSnapContainer {
public:
    SeqSnapContainer() : head_(0) { container_.resize(nSteps); }

    void takeAndSaveSnapshot(Configuration<LatticeType>& config) {
        if (head_ >= nSteps) return; // avoid overflow
        snap_.takeSnap(config);
        container_[head_] = snap_;
        head_++;
    }

    void print() const {
        for (size_t i = 0; i < head_; i++) {
            container_[i].print();
        }
    }

    void writeTxt(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) throw std::runtime_error("Cannot open file " + filename);

        for (size_t s = 0; s < head_; s++) {
            const auto& snap = container_[s];
            ofs << snap.nBackbones << "\n";
            for (const auto& seq : snap.sequences) {
                ofs << seq.size() << " ";
                for (auto val : seq) ofs << val << " ";
                ofs << "\n";
            }
            ofs << "\n";
        }
        ofs.close();
    }

    void writeHDF5(const std::string& filename) const {
        using namespace HighFive;
        File file(filename, File::Overwrite);

        // Find the maximum sequence length across all backbones and snapshots
        size_t maxLen = 0;
        size_t maxBackbones = 0;
        for (size_t s = 0; s < head_; s++) {
            maxBackbones = std::max(maxBackbones, container_[s].nBackbones);
            for (const auto& seq : container_[s].sequences) {
                maxLen = std::max(maxLen, seq.size());
            }
        }

        // Build rectangular 3D array: [snapshots, backbones, seq_length]
        std::vector<std::vector<std::vector<int>>> all_data(
            head_, std::vector<std::vector<int>>(maxBackbones, std::vector<int>(maxLen, -1))
        );

        for (size_t s = 0; s < head_; s++) {
            const auto& snap = container_[s];
            for (size_t b = 0; b < snap.sequences.size(); b++) {
                const auto& seq = snap.sequences[b];
                for (size_t i = 0; i < seq.size(); i++) {
                    all_data[s][b][i] = seq[i];
                }
            }
        }

        // Flatten into a 2D matrix (snapshots × (backbones*seq_length))
        std::vector<std::vector<int>> flat_data;
        flat_data.reserve(head_);
        for (size_t s = 0; s < head_; s++) {
            std::vector<int> row;
            row.reserve(maxBackbones * maxLen);
            for (size_t b = 0; b < maxBackbones; b++) {
                row.insert(row.end(), all_data[s][b].begin(), all_data[s][b].end());
            }
            flat_data.push_back(std::move(row));
        }

        std::vector<size_t> dims = {head_, maxBackbones * maxLen};
        DataSet dataset = file.createDataSet<int>("sequences", DataSpace(dims));
        dataset.write(flat_data);
    }

    void resetHead() { head_ = 0; }

private:
    seqSnapshot<LatticeType> snap_;
    size_t head_;
    std::vector<seqSnapshot<LatticeType>> container_;
};

    
#endif // SNAPPER_H
