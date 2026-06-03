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
#GFILEPATH=../sources/CGFD-USP-Create-Mesh/vtx-mpas-meshes/mesh.grid.nc
GFILEPATH=${DATAIN}/fixed/lat_-21_lon_-55_ellipse_a_2500_b_3688_ang_-45_iradius_5000_margin_2000_hres_3_lres_1000.region.grid.nc
#GFILEPATH=/mnt/beegfs/guilherme.mendonca/scripts_CD-CT/circle_ellipse.grid.nc
## Output directory and filename to save plot
POSTFILEPATH=${DATAIN}/fixed/ellipse_large_global.grid.png
#---------------------------------------------------------------------

source ~/.bashrc
conda activate vtx_env
python3 ${SOURCES}/CGFD-USP-Post-Proc/mpas_plot_grid.py -g $GFILEPATH -o $POSTFILEPATH
#python3 ${SOURCES}/CGFD-USP-Post-Proc/mpas_plot_grid.py -g $GFILEPATH -o $POSTFILEPATH -vmax 15 -vmin 3 -lat_min -60 -lat_max 25 -lon_min -110 -lon_max -10
