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
#include "config_MH_VSHD_square.hpp"

using Position = typename SquareLattice::Position;

int main() {
	
    /**######################### Parameters #########################**/
    //~ const int LATTICE_SIZE = 10;
    //~ const int REPEAT = 10;
    //~ const int STEPS = 1000000;
    //~ const double T = 1.0;
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
    SquareMetropolisSampler SMCS(SEED, 1./T);
    
    // Configuration on the moves and their weights
    SMCS.setWeightMovei(0, 5.);
    auto corner_move = SMCS.getMovei<CornerSquare>(0);
    corner_move->setBackboneDist(1,1);
    corner_move->setMonomerDist(1,4);
    
    SMCS.setWeightMovei(1, 3.);
    auto tail_move = SMCS.getMovei<TailSquare>(1);
    tail_move->setBackboneDist(1, 1);
    
    SMCS.setWeightMovei(2, 1);
    auto crankshaft_move = SMCS.getMovei<CrankshaftSquare>(2);
    crankshaft_move->setBackboneDist(1,1);
    crankshaft_move->setMonomerDist(0,2);
    
    SMCS.setWeightMovei(3, 1.);
    auto translation_move = SMCS.getMovei<TranslationSquare>(3);
    translation_move->setBackboneDist(0,0);
    
    SMCS.printWeights();
	
	
	
	
	// Creation of the containers for energies positions and contacts
    EnergySnapContainer<SquareLattice, STEPS> energy_container;
    PosSnapContainer<SquareLattice, STEPS> pos_container;
    ContactSnapContainer<SquareLattice, STEPS> contact_container;
    
    
    
    
    
    //---------------- Set up scheduler for this run -----------------//
    TaskScheduler<Configuration<SquareLattice>> scheduler;
	
    // MCMC sampling
    scheduler.addTask([&](Configuration<SquareLattice>& cfg, int) {
        SMCS.sample(cfg);
    });
	
    // Energy snapshot every step
    scheduler.addTask([&](auto& cfg, int step) {
        if (step % 1 == 0) energy_container.takeAndSaveSnapshot(cfg);
    });

    // Position snapshot every step
    scheduler.addTask([&](auto& cfg, int step) {
        if (step % 1 == 0) pos_container.takeAndSaveSnapshot(cfg);
    });
	
    // Contact snapshot every step
    scheduler.addTask([&](auto& cfg, int step) {
        if (step % 1 == 0) contact_container.takeAndSaveSnapshot(cfg);
    });
	
    // cumulative contact count for contact free energy computation
    int c = 0;
    scheduler.addTask([&](auto& cfg, int /**step**/) {
        int step_c = cfg.contact_matrix.get(0,1);
        if (step_c != 0) { c++; }
    });






    /**######################## Simulation ##########################**/
    double avg_cpu_time = 0;
    double avg_dF = 0;

    // Creation and initialisation of the dF file
    std::string results_filename  = "results_VSHD_square.txt";
    std::filesystem::path results_path = std::filesystem::path(DATA_DIR) / results_filename;
    std::ofstream result_file(results_path);
    result_file << "Run REPEATed " << REPEAT << " times\n";
    

	// Main loop
    for (int j = 0; j < REPEAT; j++) {
        c = 0;
		
		std::clock_t start = std::clock();
		
		// Execute scheduler tasks nSteps times
        scheduler.run_nSteps(config, STEPS);

		std::clock_t end = std::clock();
        double cpu_time = 1000.0 * (end - start) / CLOCKS_PER_SEC; // ms
        avg_cpu_time += cpu_time;

        SMCS.printStatistics();
        
        std::string energies_filename  = "VSHD_energies_"  + std::to_string(j)+".h5";
		std::string positions_filename = "VSHD_positions_" + std::to_string(j)+".h5";
		std::string contacts_filename = "VSHD_contacts_" + std::to_string(j)+".h5";
        std::filesystem::path energies_path = std::filesystem::path(DATA_DIR) / energies_filename;
		std::filesystem::path positions_path = std::filesystem::path(DATA_DIR) / positions_filename;
		std::filesystem::path contacts_path = std::filesystem::path(DATA_DIR) / contacts_filename;
        
        energy_container.writeHDF5(energies_path.string());
        pos_container.writeHDF5(positions_path.string());
        contact_container.writeHDF5(contacts_path.string());

        double dF = -std::log(1/(-1+1/(static_cast<double>(c)/STEPS)));
        avg_dF += dF;

        std::cout << "CPU time: " << cpu_time << " ms\n";
        std::cout << "dF: " << dF << std::endl;

        result_file << "Run " << j+1 << " dF: " << dF << "\n";

        energy_container.resetHead();
        pos_container.resetHead();
        
        SMCS.resetStatistics();
    }

    std::cout << "Average CPU time: " << avg_cpu_time / REPEAT << " ms\n";
    std::cout << "Average dF: " << avg_dF / REPEAT << std::endl;

    // append average dF to result file
    result_file << avg_dF / REPEAT << "\n";

    result_file.close();
}
