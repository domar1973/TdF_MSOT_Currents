#!/usr/bin/env python
#####################################################################################
# Uso:
# > ./crear_marea_tpxo9.py dominio comp1 comp2 ... compN
#####################################################################################
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os

MSOT_HOME     = os.environ["MSOT_HOME"]
carpeta_marea = MSOT_HOME + '/TPXO9/DATA/'
h_filename    = carpeta_marea + 'h_tpxo9.v5a.nc'

def loncorregida(lon):
    lonc = lon.copy()
    for k in range(len(lon)):
        if (lon[k] > 180):
            lonc[k] = lon[k] - 360
    return lonc

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

######################################################################################
# Armo lista con las componentes de marea
######################################################################################
ncompmareas        = nargs-2
compmareas         = []
for k in range(ncompmareas):
    compmareas.append(sys.argv[k+2])


######################################################################################
# Abro modelo de mareas (amplitudes y fases)
######################################################################################
tpxo9tidefile = xr.open_dataset(h_filename)[["con", "hRe", "hIm", "lon_z", "lat_z"]]
byte_array = np.array(tpxo9tidefile.con.values)
componentes_de_marea = [item.decode('utf-8').strip() for item in byte_array]
tpxo9tidefile = tpxo9tidefile.assign_coords(nc = componentes_de_marea, 
                                            nx = tpxo9tidefile.nx.data,
                                            ny = tpxo9tidefile.ny.data)

#####################################################################################
# Reasigno coordenadas a TPXO9 e interpolo en dominio grd
#####################################################################################
tpxo9tidefile = tpxo9tidefile.assign_coords(lat_z = tpxo9tidefile.lat_z[0,:].data, 
                                            lon_z = loncorregida(tpxo9tidefile.lon_z[:,0].data))
tpxo9tidefile["hRe"] = (["nc", "lon_z", "lat_z"], tpxo9tidefile.hRe.data)
tpxo9tidefile["hIm"] = (["nc", "lon_z", "lat_z"], tpxo9tidefile.hIm.data)
sshR = tpxo9tidefile.hRe.interp(lon_z = grd.lon_rho, lat_z = grd.lat_rho, method="linear").drop(["lon_z", "lat_z"])
sshI = tpxo9tidefile.hIm.interp(lon_z = grd.lon_rho, lat_z = grd.lat_rho, method="linear").drop(["lon_z", "lat_z"])
marea_compleja_interp = sshR - 1j * sshI

### Filtrar nc con compmareas
mask = marea_compleja_interp.nc.isin(compmareas)
marea_compleja_componentes_filtradas = marea_compleja_interp.where(mask, drop=True)

#####################################################################################
# Creo Xarray de amplitudes y fases y guardo a .nc
#####################################################################################
Eamp = np.abs(marea_compleja_componentes_filtradas)
Ephase = np.nan * marea_compleja_componentes_filtradas
Ephase[:] = np.mod(np.angle(marea_compleja_componentes_filtradas), 2 * np.pi)
Ephase = Ephase.astype(float) * 180 / np.pi
Eamp.name = "ha"
Ephase.name = "hp"
marea = xr.merge([Eamp, Ephase])
marea.to_netcdf("marea_tpxo9v5a_msot.nc")
