#!/usr/bin/env python
#####################################################################################
# Uso:
# > ./crear_marea_para_msot.py dominio comp1 comp2 ... compN
#####################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os

MSOT_HOME = os.environ["MSOT_HOME"]
carpeta_marea               = MSOT_HOME + '/ocean_tide_FES2014'        # Carpeta con componentes FES2014 nivel del mar
#carpeta_corrientemeridional = MSOT_HOME + '/fes2014_northward_current' # Carpeta con componentes FES2014 velocidades norte
#carpeta_corrientezonal      = MSOT_HOME + '/fes2014_eastward_current'  # Carpeta con componentes FES2014 velocidades este

nargs = len(sys.argv)
#####################################################################################
# Leo archivo batimetria pasado como parametro
#####################################################################################
if nargs > 2:
   bat_filename = sys.argv[1]
   if os.path.isfile(bat_filename):
      print(f"Procesando batimetria '{bat_filename}'.")
      grd = xr.open_dataset(bat_filename)
   else:
      print(f"Error: El archivo '{bat_filename}' no existe. Saliendo.")
      sys.exit(1)  # Termina el script con un código de salida 1
else:
   print("Error: batimetria no especificada. Saliendo.")
   sys.exit(1)  # Termina el script con un código de salida 1
print("Iniciando...")

ncompmareas        = nargs-2
compmareas         = []
mapas_interpolados = []
for i in range(ncompmareas):
    compo = sys.argv[i+2]
    compmareas.append(compo)
    mapacompo   = xr.open_dataset(carpeta_marea               + '/'+compo+'.nc')
    #mapacorri_n = xr.open_dataset(carpeta_corrientemeridional + '/'+compo+'.nc')
    #mapacorri_e = xr.open_dataset(carpeta_corrientezonal      + '/'+compo+'.nc')
    tide_file2 = mapacompo.assign_coords(lat = mapacompo.lat.data, 
                                         lon = mapacompo.lon.data - 360)
    rawret  = tide_file2['amplitude'] * np.cos(tide_file2['phase']*np.pi/180.0) / 100.0
    rawimt  =-tide_file2['amplitude'] * np.sin(tide_file2['phase']*np.pi/180.0) / 100.0
    ret     = rawret.where(rawret < 100, other=0)
    imt     = rawimt.where(rawimt < 100, other=0)
    sshR = ret.interp(lon = grd.lon_rho, lat = grd.lat_rho, method="linear").drop(["lon", "lat"])
    sshI = imt.interp(lon = grd.lon_rho, lat = grd.lat_rho, method="linear").drop(["lon", "lat"])
    marea_compleja_interp = sshR - 1j * sshI
    Eamp = np.abs(marea_compleja_interp)
    Ephase = np.nan * marea_compleja_interp
    Ephase[:] = np.mod(np.angle(marea_compleja_interp), 2 * np.pi)
    Ephase = Ephase.astype(float) * 180 / np.pi
    Eamp.name = "ha"
    Ephase.name = "hp"
    mapas_interpolados.append(xr.merge([Eamp, Ephase]))
marea = xr.concat(mapas_interpolados, pd.Index(compmareas, name='nc'))
marea.to_netcdf("marea_fes2014_msot.nc")
