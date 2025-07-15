#!/bin/bash        
# MSOT TdF (D. Badagnani, dani.en.villa.rica@gmail.com, abril 2024)
export FECHA_INIT=20190101_00 # Fecha inical y final: AnnioMesDia_Hora
export FECHA_FIN=20190101_06
#   RANK=0: Modelo de más baja resolución, forzado por marea
#   RANK=1,2,...: Modelo de mayor resolución forzado por la solución del de rango inferior inmediato
export RANK=0
export DT=12 # Paso temporal EN SEGUNDOS
export RUNDIR="." # Directorio de ejecucion
export DIR_OUTPUT="./" # Directorio de salida
export OMP_NUM_THREADS=8 # Cantidad de nodos, tiene que se consistente con lo compilado (?????)
export DOMINIO="$RUNDIR/croco_grd.nc" # Archivo con la batimetria
export tidemodel=FES2014  # Posibles: FES2014

#
# COMPILAR editando el número de componentes de marea!!!!!!!
#
if [[ $tidemodel == 'FES2014' ]]
then
  #componentes=${componentes}\ 2n2
  #componentes=${componentes}\ eps2
  #componentes=${componentes}\ j1
  componentes=${componentes}\ k1
  componentes=${componentes}\ k2
  #componentes=${componentes}\ l2
  #componentes=${componentes}\ la2
  componentes=${componentes}\ m2
  #componentes=${componentes}\ m3
  #componentes=${componentes}\ m4
  #componentes=${componentes}\ m6
  #componentes=${componentes}\ m8
  #componentes=${componentes}\ mf
  #componentes=${componentes}\ mks2
  #componentes=${componentes}\ mm
  #componentes=${componentes}\ mn4
  #componentes=${componentes}\ ms4
  #componentes=${componentes}\ msf
  #componentes=${componentes}\ msqm
  #componentes=${componentes}\ mtm
  #componentes=${componentes}\ mu2
  componentes=${componentes}\ n2
  #componentes=${componentes}\ n4
  #componentes=${componentes}\ nu2
  componentes=${componentes}\ o1
  componentes=${componentes}\ p1
  componentes=${componentes}\ q1
  #componentes=${componentes}\ r2
  #componentes=${componentes}\ s1
  componentes=${componentes}\ s2
  #componentes=${componentes}\ s4
  #componentes=${componentes}\ sa
  #componentes=${componentes}\ ssa
  #componentes=${componentes}\ t2
fi

export MIEM=01 # 

echo ########################################
echo MSOT TdF Corrientes
echo Universidad Nacional de Tierra del Fuego
echo
echo Ver output en log file
echo ########################################

echo MSOT TdF Corrientes                        >./log_run_${FECHA_INIT}_${MIEM}.out
echo Universidad Nacional de Tierra del Fuego  >>./log_run_${FECHA_INIT}_${MIEM}.out
echo                                           >>./log_run_${FECHA_INIT}_${MIEM}.out
echo                                           >>./log_run_${FECHA_INIT}_${MIEM}.out

source ./config.ini # Path a librerias necesarias
if [[ $RANK -eq 0 ]] # si es simulacion padre
then
   echo MODELO DE MAREA: $tidemodel
   echo COMPONENTES: $componentes
   echo MODELO DE MAREA: $tidemodel   >./log_run_${FECHA_INIT}_${MIEM}.out
   echo COMPONENTES: $componentes    >>./log_run_${FECHA_INIT}_${MIEM}.out
   if [[ $tidemodel == 'TPXO9' ]]
   then
      ./crear_marea_tpxo9.py $DOMINIO $componentes
      export MAREA="$RUNDIR/marea_tpxo9v5a_msot.nc"
   else
      if [[ $tidemodel == 'FES2014' ]]
      then
         ./crear_marea_fes2014.py $DOMINIO $componentes 
         export MAREA="$RUNDIR/marea_fes2014_msot.nc"
      fi
   fi
else # si es anidado es forzado por solucion padre
   export MAREA=False
fi
echo INICIO: $(date -d now)
echo INICIO: $(date -d now) >>./log_run_${FECHA_INIT}_${MIEM}.out
python  ./crear_archivos_forzantes.py # Creacion de archivos forzantes
# Ejecucion de la simulacion
$RUNDIR/pronomar-ndm-cima ${RUNDIR}/${FECHA_INIT}_${MIEM}.in  2>./log_run_${FECHA_INIT}_${MIEM}.err >>./log_run_${FECHA_INIT}_${MIEM}.out
echo FIN: $(date -d now)
echo FIN: $(date -d now) >>./log_run_${FECHA_INIT}_${MIEM}.out
