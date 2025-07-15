###### SCRIPT DESARROLLADO POR M. SACCO - M DINAPOLI - M DE OTO ####

export LC_ALL="C"

####################################################################
############### CORROBORAR PATH DE COMPILADOR ######################
####################################################################
# export COMPILERVARS_ARCHITECTUR=intel64
# export COMPILERVARS_PLATFORM=linux
# export INTEL_BASE=/home/opt/intel/
# export CVER="I19"
# export INTEL_COMPILER_TOPDIR=$INTEL_BASE/parallel_studio_xe_2019
# export IMPI_BASE=$INTEL_BASE/impi/2019.1.144/
# source $INTEL_COMPILER_TOPDIR/psxevars.sh $COMPILERVARS_ARCHITECTUR
# source $IMPI_BASE/intel64/bin/mpivars.sh $COMPILERVARS_ARCHITECTUR
####################################################################

############################################
########## DEFINICIONES DE PATHS ###########
############################################
# export RAIZ='/data/share/bin/'
# export DIRLIB="$RAIZ/WRF-4.0_${CVER}/LIBRARIES"
# export NETCDF=$DIRLIB/netcdf-4.4.1_${CVER}
export NETCDF="/usr/local/netcdf/ /usr/local/netcdff/"
####################################################

export PATH=$NETCDF/bin:$PATH

export LD_LIBRARY_PATH=$NETCDF/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$LD_LIBRARY_PATH:$LD_RUN_PATH
export INCLUDE=$NETCDF/include
#############################################
