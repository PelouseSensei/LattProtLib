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
#include "config_MH_PULL_cubic.hpp"

using Position = typename CubicLattice::Position;

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
    Configuration<CubicLattice> config(LATTICE_SIZE, BoundaryType::Periodic, matrix);
    config.loadFromCIFile(BB_PATH);
    config.describe();
    config.printOccupancy();
	
	// Creation of a cubic metropolis sampler
    CubicMetropolisSampler_pull SMCS(SEED, 1./T);
    
    // Configuration on the moves and their weights
    SMCS.setWeightMovei(0, 1.);
    auto pullb_move = SMCS.getMovei<PullBackwardCubic>(0);
    pullb_move->setBackboneDist(1,1);
    pullb_move->setMonomerDist(1, 4);
    
    SMCS.setWeightMovei(1, 1.);
    auto pullf_move = SMCS.getMovei<PullForwardCubic>(1);
    pullf_move->setBackboneDist(1,1);
    pullf_move->setMonomerDist(1, 4);
    
    SMCS.setWeightMovei(2, 21./8.);
    auto pulle_move = SMCS.getMovei<PullEndCubic>(2);
    pulle_move->setBackboneDist(1,1);
    
    SMCS.setWeightMovei(3, 1.);
    auto translation_move = SMCS.getMovei<TranslationCubic>(3);
    translation_move->setBackboneDist(0, 0);
    
    SMCS.setWeightMovei(4, 1.);
    auto rotation_move = SMCS.getMovei<RotationCubic>(4);
    rotation_move->setBackboneDist(1,1);
    rotation_move->setDistMono(0,2);
    
    SMCS.setWeightMovei(5, 1.);
    auto symmetry_move = SMCS.getMovei<SymmetryCubic>(5);
    symmetry_move->setBackboneDist(0,0);
    
    SMCS.printWeights();
	
	
	
	
	// Creation of the containers for energies positions and contacts
    EnergySnapContainer<CubicLattice, STEPS> energy_container;
    PosSnapContainer<CubicLattice, STEPS> pos_container;
    ContactSnapContainer<CubicLattice, STEPS> contact_container;
    
    
    
    
    
    //---------------- Set up scheduler for this run -----------------//
    TaskScheduler<Configuration<CubicLattice>> scheduler;
	
	//~ int a;
    // MCMC sampling
    scheduler.addTask([&](Configuration<CubicLattice>& cfg, int /*step*/) {
        SMCS.sample(cfg);
        //~ config.describe();
        //~ config.printOccupancy();
        //~ std::cin >> a;
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
    std::string results_filename  = "results_PULL_cubic.txt";
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
        
        std::string energies_filename  = "PULL_energies_"  + std::to_string(j)+".h5";
		std::string positions_filename = "PULL_positions_" + std::to_string(j)+".h5";
		std::string contacts_filename = "PULL_contacts_" + std::to_string(j)+".h5";
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
