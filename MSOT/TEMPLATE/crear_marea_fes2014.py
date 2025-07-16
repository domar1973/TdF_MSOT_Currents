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
carpeta_corrientemeridional = MSOT_HOME + '/fes2014_northward_current' # Carpeta con componentes FES2014 velocidades norte
carpeta_corrientezonal      = MSOT_HOME + '/fes2014_eastward_current'  # Carpeta con componentes FES2014 velocidades este

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
    mapacorri_n = xr.open_dataset(carpeta_corrientemeridional + '/'+compo+'.nc')
    mapacorri_e = xr.open_dataset(carpeta_corrientezonal      + '/'+compo+'.nc')
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
    # Corrientes: cargo amplitud y fase de FES2014
    u_file = mapacorri_e.assign_coords(lat = mapacorri_e.lat.data,
                                       lon = mapacorri_e.lon.data - 360)
    v_file = mapacorri_n.assign_coords(lat = mapacorri_n.lat.data,
                                       lon = mapacorri_n.lon.data - 360)
    Ureal = (u_file['Ua'] * np.cos(u_file['Ug']*np.pi/180.0) / 100.0)
    Uimag = (-u_file['Ua'] * np.sin(u_file['Ug']*np.pi/180.0) / 100.0)
    Vreal = (v_file['Va'] * np.cos(v_file['Vg']*np.pi/180.0) / 100.0)
    Vimag = (-v_file['Va'] * np.sin(v_file['Vg']*np.pi/180.0) / 100.0)
    Ureal_i = Ureal.interp(lon = grd.lon_rho, lat = grd.lat_rho, method="linear").drop(["lon", "lat"])
    Uimag_i = Uimag.interp(lon = grd.lon_rho, lat = grd.lat_rho, method="linear").drop(["lon", "lat"])
    Vreal_i = Vreal.interp(lon = grd.lon_rho, lat = grd.lat_rho, method="linear").drop(["lon", "lat"])
    Vimag_i = Vimag.interp(lon = grd.lon_rho, lat = grd.lat_rho, method="linear").drop(["lon", "lat"])
    U_complex = Ureal_i - 1j * Uimag_i
    V_complex = Vreal_i - 1j * Vimag_i
    wp = (U_complex + 1j*V_complex)/2
    wm = (U_complex - 1j*V_complex)/2
    Cmax = np.abs(wp) + np.abs(wm)
    Cmin = np.abs(wp) - np.abs(wm)
    Cangle = (np.angle(wp) - np.angle(wm)) * 0.5 * 180/np.pi
    Cphase = (np.angle(wp) + np.angle(wm)) * 0.5 * 180/np.pi
    Cangle = np.mod(Cangle, 360)
    Cphase = np.mod(Cphase, 360)
    Cmax_da = Cmax.rename("tide_Cmax")
    Cmin_da = Cmin.rename("tide_Cmin")
    Cangle_da = Cangle.rename("tide_Cangle")
    Cphase_da = Cphase.rename("tide_Cphase")
    mapas_interpolados.append(xr.merge([Eamp, Ephase, Cmax_da, Cmin_da, Cangle_da, Cphase_da]))
marea = xr.concat(mapas_interpolados, pd.Index(compmareas, name='nc'))
marea.to_netcdf("marea_fes2014_msot.nc")
