#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <cmath>
#include <filesystem>
#include "lattice/lattice.hpp"
#include "lattice/boundaryconditions.hpp"
#include "system/configuration.hpp"
#include "system/snapper.hpp"
#include "moves/moves.hpp"
#include "samplers/mcmc_samplers.hpp"
#include "samplers/WL_samplers.hpp"
#include "schedulers/task_scheduler.hpp"
//~ #include "config_WL_PULL_square.hpp"

using Position = typename SquareLattice::Position;

int main() {
	
	/**######################### Parameters #########################**/
    const int LATTICE_SIZE = 10;
    const int REPEAT = 10;
    const int STEPS = 1000000;
    const double T = 1.0;
    const int SEED = 1245;
    const std::string INTER_PATH = "config_files/HP.inter";
    const std::string BB_PATH = "config_files/backbones_CI.bb";


	auto discrete_energy_to_bin = [](double E){ 
		int b = static_cast<int>(E+5); 
		b = std::abs(b);
		return b;
	};


    /**####################### Initialization #######################**/
    
    // Creating and initiliasing interaction matrix
    InterMatrix matrix = InterMatrix::fromInterFile(INTER_PATH);
	
	// Creating and initialising configuration
    Configuration<SquareLattice> config(LATTICE_SIZE, BoundaryType::Periodic, matrix);
    config.loadFromCIFile(BB_PATH);
    config.describe();
    config.printOccupancy();
    
    //~ // Creating the task scheduler
    //~ TaskScheduler<Configuration<SquareLattice>> scheduler;
	//~ // add a logging task
	//~ scheduler.addTask([](Configuration<SquareLattice>& cfg, int step){
		//~ if ((step % 100000) == 0) {
			//~ std::cout << "step " << step << ", energy = " << cfg.getEnergy() << "\n";
		//~ }
	//~ });
	
	// Creation of a square metropolis sampler
    //~ SquareWangLandauSampler_pull Sampler(SEED, 1.0, -5.5, 0.5, 6, discrete_energy_to_bin);//, &scheduler);
    SquareWangLandauSampler_pull Sampler(SEED, 1.0, -10.5, 0.5, 11);//, &scheduler);
    
    // Configuration on the moves and their weights
    Sampler.setWeightMovei(0, 1.);
    auto pullb_move = Sampler.getMovei<PullBackwardSquare>(0);
    pullb_move->setBackboneDist(1,1);
    pullb_move->setMonomerDist(1, 4);
    
    Sampler.setWeightMovei(1, 1.);
    auto pullf_move = Sampler.getMovei<PullForwardSquare>(1);
    pullf_move->setBackboneDist(1,1);
    pullf_move->setMonomerDist(1, 4);
    
    Sampler.setWeightMovei(2, 7./4);
    auto pulle_move = Sampler.getMovei<PullEndSquare>(2);
    pulle_move->setBackboneDist(1,1);
    
    Sampler.setWeightMovei(3, 1.);
    auto translation_move = Sampler.getMovei<TranslationSquare>(3);
    translation_move->setBackboneDist(0, 0);
    
    Sampler.setWeightMovei(4, 1.);
    auto rotation_move = Sampler.getMovei<RotationSquare>(4);
    rotation_move->setBackboneDist(1,1);
    rotation_move->setDistMono(0,2);
    
    Sampler.setWeightMovei(5, 0.);
    auto symmetry_move = Sampler.getMovei<SymmetrySquare>(5);
    symmetry_move->setBackboneDist(0,0);
    
    Sampler.printWeights();
	
	
	Sampler.initialize(config); // mark starting energy

	Sampler.run(config, 1000000); // runs steps (will call scheduler after each step)
	Sampler.printHistogram();
	Sampler.printLogDOS();
	
	Sampler.saveToHDF5("data/WL_test.h5");

}
