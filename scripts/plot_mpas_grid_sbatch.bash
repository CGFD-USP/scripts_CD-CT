#!/bin/bash
#-----------------------------------------------------------------------------#
# !SCRIPT: run_post_on_mpas_grid
#
# !DESCRIPTION:
#     Script to run the postprocessing of MONAN model using the original MPAS grid.
#
#-----------------------------------------------------------------------------#


# Set environment variables exports:
echo ""
echo -e "\033[1;32m==>\033[0m Moduling environment for MONAN model...\n"
. setenv.bash

# Standart directories variables:---------------------------------------
DIRHOMES=${DIR_SCRIPTS}/scripts_CD-CT; mkdir -p ${DIRHOMES}  
DIRHOMED=${DIR_DADOS}/scripts_CD-CT;   mkdir -p ${DIRHOMED}  
export SCRIPTS=${DIRHOMES}/scripts;    mkdir -p ${SCRIPTS}
DATAIN=${DIRHOMED}/datain;             mkdir -p ${DATAIN}
DATAOUT=${DIRHOMED}/dataout;           mkdir -p ${DATAOUT}
SOURCES=${DIRHOMES}/sources;           mkdir -p ${SOURCES}
EXECS=${DIRHOMED}/execs;               mkdir -p ${EXECS}
#----------------------------------------------------------------------


# Local variables------------------------------------------------------
## Grid file
GFILEPATH=${DATAIN}/fixed/lat_-15_lon_-65_ellipse_a_4000_b_6000_ang_30_iradius_2000_margin_800_hres_5_lres_15.region.grid.nc
## Output directory and filename to save plot
POSTFILEPATH=${DATAIN}/fixed/lat_-15_lon_-65_ellipse_a_4000_b_6000_ang_30_iradius_2000_margin_800_hres_5_lres_15.region.grid_with_optm.parallelized.png
#---------------------------------------------------------------------

cat << EOF0 > plot_grid.bash
#!/bin/bash -x
#SBATCH --job-name=${GRID_jobname}
#SBATCH --nodes=${GRID_nnodes}
#SBATCH --ntasks=${GRID_ntasks}
#SBATCH --cpus-per-task=${GRID_ncpt}
#SBATCH --partition=${GRID_QUEUE}
#SBATCH --time=${GRID_walltime}
#SBATCH --output=${DATAIN}/fixed/plot_grid.bash.o%j    # File name for standard output
#SBATCH --error=${DATAIN}/fixed/plot_grid.bash.e%j     # File name for standard error output
#SBATCH --exclusive
#SBATCH --mem=${GRID_memory}G

echo "Loading anaconda..."
module load anaconda3-2022.05-gcc-11.2.0-q74p53i
echo "Making sure we have access to conda env..."
conda init
source ~/.bashrc
conda config --add envs_dirs /home/guilherme.mendonca/.conda/envs
echo "Activating conda env..."
conda activate vtx_env
echo "Running mpas_plot_grid.py..."
time python3 ${SOURCES}/CGFD-USP-Post-Proc/mpas_plot_grid.py -g ${GFILEPATH} -o ${POSTFILEPATH} -nc ${GRID_ntasks} -vmin 4 -vmax 19

EOF0
chmod a+x plot_grid.bash

echo -e  "${GREEN}==>${NC} Submitting script to plot grid and waiting for finishing before exiting... \n"
echo -e  "${GREEN}==>${NC} Logs being generated at ${DATAIN}/fixed... \n"
echo -e  "sbatch ${SCRIPTS}/plot_grid.bash"
sbatch --wait plot_grid.bash
mv plot_grid.bash ${DATAIN}/fixed
