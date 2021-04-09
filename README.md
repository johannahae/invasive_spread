# invasive_spread 

** Code and data accompanying the paper "Invasive spread in meta-food-webs depends on landscape structure, fertilization and species characteristics" by Johanna Häussler, Remo Ryser and Ulrich Brose, 2021, Oikos, DOI: 10.1111/oik.07503 ** 

The directory "code/" contains all the source code in C++, R and Bash needed to replicate the results in theh manuscript. In detail, this means to generate the food webs (adjacency matrices and body mass vectors), landscapes (XY-coordinates) and to run the simulations. 
Successfully running these files requires:

* Installation of C++ and the libraries -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial. 
* An installion of R and several additional R packages. For the manuscript R versions 3.5.1 and 3.6.3. were used for simulations respecitively analysis. Package requirements can be checked with the Rscript "code/package_requirements.R" from the R package "requirements"). 

# CODE INSTRUCTIONS (run from code/)
Step 1: Create executables with ./cmake (creates "webs" and "simulation")
Step 2: Make food webs with ./build_webs.sh (saved in data/webs/)
Step 3: Make invasion webs with remove_spp.R (saved in data/webs)
Step 4: Make landscapes with make_landscapes.R (saved in data/landscapes/)
Step 5: Make bash files to simulate invasion processes with make_bash_files.R
Step 6: Simulate invasion processes by sourcing each bash-file, e.g. source 88516_001.sh etc. 
=> This runs the file invasion_process.sh:
I) executes ./simulation for the specified settings WITHOUT the invasive species (data/output/*_0.csv); 
II) use simulation output to make invasion input (data/invasion/) with make_invasion_input.R; 
III) executes ./simulation for the specified settings WITH the invasive species (data/output/*_1.csv). 
Step 7: Summarize the simulation output with make_summary.R (data/results/summary.csv)
Step 8: Prepare the summarized data for further analysis with prepare_data_plot.R (data/results/dat_plots.rds)
Step 9: Create figures and run statistical tests with the corresponding Rscripts (plot_figure1.R, plots_abiotic.R, plots_invaderchars.R, stats_invaderchars.R, plots_SI.R). 

# DATA STRUCTURE
data/webs/ --- 5 food webs generated with code/build_webs.sh and modified with code/remove_spp.R
data/landscapes/ --- 20 landscapes generated with code/make_landscapes.R
data/output/ --- simulation output with and without invasion 
data/invasion/ --- output from code/make_invasion_input.R (= input for simulations with invasion)
data/results/ ---- output summarized for analysis and plotting

All included code is covered under the GPL, version 3, available here: (http://www.gnu.org/licenses/gpl-3.0.en.html).

Any question about the code can be sent to Johanna Häussler (johanna.haeussler@idiv.de). 

