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
#include "schedulers/task_scheduler.hpp"
#include "config_SA_TM_square.hpp"

using Position = typename SquareLattice::Position;

int main() {
	
    /**######################### Parameters #########################**/
    //~ const int LATTICE_SIZE = 10;
    //~ const int REPEAT = 10;
    //~ const int K = 20;
    //~ const int TK_STEPS = 50000;
    double T = 2.0;
    //~ const double T_START = 2.0;
    //~ const double ALPHA = 0.9;
    //~ const int SEED = 1245;





    /**####################### Initialization #######################**/
    
    // Creating and initiliasing interaction matrix
    InterMatrix matrix = InterMatrix::fromInterFile(INTER_PATH);
	
	// Creating and initialising configuration
    Configuration<SquareLattice> config(LATTICE_SIZE, BoundaryType::Periodic, matrix);
    config.loadFromCIFile(BB_PATH);
    config.describe();
    config.printOccupancy();
	
	// Creation of a square metropolis sampler
    MHSequenceSampler<SquareLattice> MSS(SEED, T_START);
    
    // Configuration on the moves and their weights
    MSS.setWeightMovei(0, 1.);
    auto mutation_move = MSS.getMovei<OnePointMutation<SquareLattice>>(0);
    mutation_move->setBackboneDist(1,1);
    mutation_move->setMonomerDist(0,5);
    mutation_move->setTypeDist(0,5);
    
    MSS.setWeightMovei(1, 1.);
    auto type_switch_move = MSS.getMovei<TwoPointSwitch<SquareLattice>>(1);
    type_switch_move->setBackboneDist(1, 1);
    type_switch_move->setMonomerDist(0,5);
    
    MSS.printWeights();
	
	
	
	
	// Creation of the containers for energies positions and contacts
    EnergySnapContainer<SquareLattice, TK_STEPS*K> energy_container;
    SeqSnapContainer<SquareLattice, TK_STEPS*K> sequence_container;
    
    
    
    
    
    //---------------- Set up scheduler for this run -----------------//
    TaskScheduler<Configuration<SquareLattice>> scheduler;
	
    // MCMC sampling
    scheduler.addTask([&](Configuration<SquareLattice>& cfg, int) {
        MSS.sample(cfg);
    });
	
    // Energy snapshot every step
    scheduler.addTask([&](auto& cfg, int step) {
        if (step % 1 == 0) energy_container.takeAndSaveSnapshot(cfg);
    });
    
    // Sequence snapshot every step
    scheduler.addTask([&](auto& cfg, int step) {
        if (step % 1 == 0) sequence_container.takeAndSaveSnapshot(cfg);
    });
    
    //~ // Debug task
    //~ int a;
    //~ scheduler.addTask([&](auto& cfg, int step) {
		//~ double e = cfg.getEnergy();
		//~ if (e < -1.0) {
			//~ throw std::runtime_error("Energy below 5, stopping simulation");
		//~ }
		//~ std::cin >> a;
	//~ });






    /**######################## Simulation ##########################**/
    double avg_cpu_time = 0;

	// Main loop
    for (int j = 0; j < REPEAT; j++) {
        //~ c = 0;
		
		std::clock_t start = std::clock();
		
		for (int k =0; k<K; k++) {
			T = T_START * std::pow(ALPHA, k); // Geometric schedule for temperature
			std::cout << T << std::endl;
			MSS.setBeta(T);
			scheduler.run_nSteps(config, TK_STEPS);
			
		}

		std::clock_t end = std::clock();
        double cpu_time = 1000.0 * (end - start) / CLOCKS_PER_SEC; // ms
        avg_cpu_time += cpu_time;

        MSS.printStatistics();
        
        std::string energies_filename  = "SEQ_energies_"  + std::to_string(j)+".h5";
		std::string sequences_filename = "SEQ_sequences_" + std::to_string(j)+".h5";
        std::filesystem::path energies_path = std::filesystem::path(DATA_DIR) / energies_filename;
		std::filesystem::path sequences_path = std::filesystem::path(DATA_DIR) / sequences_filename;
        
        energy_container.writeHDF5(energies_path.string());
        sequence_container.writeHDF5(sequences_path.string());

        std::cout << "CPU time: " << cpu_time << " ms\n";

        energy_container.resetHead();
        sequence_container.resetHead();
        
        MSS.setBeta(T_START);
        //~ k = 0;
        
        MSS.resetStatistics();
    }

    std::cout << "Average CPU time: " << avg_cpu_time / REPEAT << " ms\n";
}
