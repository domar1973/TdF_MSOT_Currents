#!/usr/bin/env python
#####################################################################################
# Matias Dinapoli, modificado por Daniel Badagnani
# Uso:
# > ./crear_marea_para_msot.py ocean_grd_<version>.nc
#####################################################################################
import xarray as xr
import numpy as np
import sys
import os

MSOT_HOME = os.environ["MSOT_HOME"]
archivo_marea = MSOT_HOME + "/marea_tpxo9v5a_austral.nc"

#####################################################################################
# Leo archivo batimetria pasado como parametro
#####################################################################################
if len(sys.argv) > 1:
   bat_filename = sys.argv[1]
   if os.path.isfile(bat_filename):
      print(f"Procesando batimetria '{bat_filename}'.")
   else:
      print(f"Error: El archivo '{bat_filename}' no existe. Saliendo.")
      sys.exit(1)  # Termina el script con un código de salida 1
else:
   print("Error: batimetria no especificada. Saliendo.")
   sys.exit(1)  # Termina el script con un código de salida 1
print("Iniciando...")


#####################################################################################
# Lectura de los campos del modelo TPXO y el dominio de interes
#####################################################################################
# Leo el archivo de marea de TPOX9v5a
componentes_de_marea = ["m2", "s2", "n2", "k2", "k1", "o1", "p1", "q1"]
tide_file = xr.open_dataset(archivo_marea)[["hRe", "hIm", "lon_z", "lat_z"]]
tide_file = tide_file.assign_coords(nc = componentes_de_marea, 
                                    nx = tide_file.nx.data,
                                    ny = tide_file.ny.data)
# Leo la grilla
#grd = xr.open_dataset("./ocean_grd.nc")
grd = xr.open_dataset(bat_filename)

#####################################################################################
# Interpolacion
#####################################################################################
# Interpolo
tide_file = tide_file.assign_coords(lat_z = tide_file.lat_z[:,0].data, 
                                    lon_z = tide_file.lon_z[0,:].data - 360)
tide_file["hRe"] = (["nc", "lat_z", "lon_z"], tide_file.hRe.data)
tide_file["hIm"] = (["nc", "lat_z", "lon_z"], tide_file.hIm.data)
sshR = tide_file.hRe.interp(lon_z = grd.lon_rho, lat_z = grd.lat_rho, method="linear").drop(["lon_z", "lat_z"])
sshI = tide_file.hIm.interp(lon_z = grd.lon_rho, lat_z = grd.lat_rho, method="linear").drop(["lon_z", "lat_z"])
marea_compleja_interp = sshR - 1j * sshI

#####################################################################################
# Creacion del campo para alimentar a CROCO
#####################################################################################
# Calculo la amplitud y fase
Eamp = np.abs(marea_compleja_interp)
Ephase = np.nan * marea_compleja_interp
Ephase[:] = np.mod(np.angle(marea_compleja_interp), 2 * np.pi)
Ephase = Ephase.astype(float) * 180 / np.pi
# Creo el netcdf y lo guardo
### Eamp.name = "Eamp"
### Ephase.name = "Ephase"
Eamp.name = "ha"
Ephase.name = "hp"
marea = xr.merge([Eamp, Ephase])
marea.to_netcdf("marea_tpxo9v5a_msot.nc")
