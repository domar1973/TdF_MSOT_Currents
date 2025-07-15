#!/usr/bin/env python3
################################################################################
# Liberias
################################################################################
import xarray as xr                     # Uso de campos
import pandas as pd                     # Uso de dataframes y fechas
import numpy as np                      # Operaciones matematicas
import os                               # Instrucciones a linux
import glob                             # Lectura de archivos en el SO
from pyproj import Proj, Geod           # Para la creacion de la grilla

################################################################################

################################################################################
# Definciones generales
################################################################################
formato_netcdf = "netcdf4"
dominio = "../croco_grd.nc" # COMENTAR CARACTERISTICA
MSOT_HOME = os.environ["MSOT_HOME"]
GEBCO = MSOT_HOME+"/gebco_2023_n-20.0_s-70.0_w-105.0_e-10.0.nc"

Lat_i = -60   # Límite Sur
Lat_f = -45   # Limite Norte

Lon_i = -70   # Limite Oeste
Lon_f = -50   # Limite Este

DivisorLat = 24 # 1/24: reso 5 Km
DivisorLon = 12 # 1/12: reso 5 Km en TdF

#################################################################################
# Bordes del dominio
#################################################################################
#longitudes = np.arange(Lon_i, Lon_f + 1/DivisorLon, 1/DivisorLon)
longitudes = np.linspace(Lon_i, Lon_f, num = int((Lon_f - Lon_i)*DivisorLon + 1))
#latitudes  = np.arange(Lat_i, Lat_f + 1/DivisorLat, 1/DivisorLat)
latitudes  = np.linspace(Lat_i, Lat_f, num = int((Lat_f - Lat_i)*DivisorLat + 1))

Nlat = int((Lat_f - Lat_i)*DivisorLat+1)
Nlon = int((Lon_f - Lon_i)*DivisorLon+1)

################################################################################
# Dominio
################################################################################
def make_grd(nombre = dominio, 
             longitudes = longitudes, latitudes = latitudes,
             formato_netcdf = formato_netcdf, GEBCO = GEBCO):
    #################################################################################
    # Arakawa-C
    #################################################################################
    lon_rho, lat_rho = np.meshgrid(longitudes, latitudes)
    # Coordenadas en puntos U
    lon_u, lat_u = 0.5 * (lon_rho[:,1:] + lon_rho[:,:-1]), 0.5 * (lat_rho[:,1:] + lat_rho[:,:-1])
    # Coordenadas en puntos V
    lon_v, lat_v = 0.5 * (lon_rho[1:,:] + lon_rho[:-1,:]), 0.5 * (lat_rho[1:,:] + lat_rho[:-1,:])
    # Coordenadas en puntos PSI
    lon_psi, lat_psi = 0.5 * (lon_u[1:,:] + lon_u[:-1,:]), 0.5 * (lat_u[1:,:] + lat_u[:-1,:])
    # Coriolis
    f = 2.0 * 7.29e-5 * np.sin(lat_rho * np.pi / 180.0)
    #################################################################################
    # Calculo de la metrica
    #################################################################################
    # Deltas X e Y para las derivadas
    dx, dy = np.zeros_like(lon_rho), np.zeros_like(lat_rho)
    geod = Geod(ellps = 'WGS84')
    azimutX, _, dx[:, 1:-1] = geod.inv(lon_u[:, 1:], lat_u[:, 1:], lon_u[:, :-1], lat_u[:, :-1])
    dx[:, 0], dx[:, -1] = dx[:, 1], dx[:, -2]
    azimutY, _, dy[1:-1, :] = geod.inv(lon_v[1:, :], lat_v[1:, :], lon_v[:-1, :], lat_v[:-1, :])
    dy[0, :], dy[-1, :] = dy[1, :], dy[-2, :]
    pm, pn = 1./dx, 1./dy
    # Derivadas en la metrixa (xi; eta)
    dndx, dmde = np.zeros_like(lon_rho), np.zeros_like(lat_rho)
    dndx[1:-1, 1:-1] = 0.5 * (dy[1:-1, 2:] - dy[1:-1, :-2]) 
    dmde[1:-1, 1:-1] = 0.5 * (dx[2:, 1:-1] - dx[:-2, 1:-1])
    #################################################################################
    # Angulos entre celdas
    #################################################################################
    angle = np.zeros_like(lon_rho)
    angle[:, 1:-1] = np.pi/2 - azimutX * np.pi/180
    angle[:, 0] = angle[:, 1]
    angle[:, -1] = angle[:, -2]
    #################################################################################
    # Batimetria
    #################################################################################
    h = -xr.open_dataarray(GEBCO)
    h = h.interp(lon = lon_rho[0, :], lat = lat_rho[:, 0])
    #################################################################################
    # Mascaras
    #################################################################################
    h = h.where(h > 1, 1)
    antimask = h.where(h == 1, 0).data
    mask = 1 - antimask
    mask[0, 165:185] = 0 ### Tapando el estrecho de Magallanes
    mask_u = mask[:, 1:] * mask[:, :-1]
    mask_v = mask[1:, :] * mask[:-1, :]
    mask_psi = mask_u[1:, :] * mask_u[:-1, :]
    #################################################################################
    # Generacion de netcdf
    #################################################################################
    Mp, Lp = lon_rho.shape
    M, L = lon_psi.shape
    xl = lon_rho.shape[0]
    el = lat_rho.shape[1]
    variables = {"xi_rho": ("one", [Lp]), "eta_rho": ("one", [Mp]),
                 "xi_u": ("one", [L]), "eta_u": ("one", [Mp]),
                 "xi_v": ("one", [Lp]), "eta_v": ("one", [M]),
                 "xi_psi": ("one", [L]), "eta_psi": ("one", [M]),
                 "xl": ("one", [xl]), "el": ("one", [el]),
                'pm': (('eta_rho', 'xi_rho'),  pm), 
                'pn': (('eta_rho', 'xi_rho'), pn), 
                'dmde': (('eta_rho', 'xi_rho'), dmde), 
                'dndx': (('eta_rho', 'xi_rho'), dndx), 
                'angle': (('eta_rho', 'xi_rho'), angle), 
                'f': (('eta_rho', 'xi_rho'), f), 
                'h': (('eta_rho', 'xi_rho'), h.data), 
                'lon_rho': (('eta_rho', 'xi_rho'), lon_rho),
                'lat_rho': (('eta_rho', 'xi_rho'), lat_rho),
                'lon_u': (('eta_u', 'xi_u'), lon_u),
                'lat_u': (('eta_u', 'xi_u'), lat_u),
                'lon_v': (('eta_v', 'xi_v'), lon_v),
                'lat_v': (('eta_v', 'xi_v'), lat_v),
                'lon_psi': (('eta_psi', 'xi_psi'), lon_psi),
                'lat_psi': (('eta_psi', 'xi_psi'), lat_psi),
                "mask_rho": (('eta_rho', 'xi_rho'), mask.data),
                "mask_u": (('eta_u', 'xi_u'), mask_u.data),
                "mask_v": (('eta_v', 'xi_v'), mask_v.data),
                "mask_psi": (('eta_psi', 'xi_psi'), mask_psi.data),}
    # Creo el archivo de la batimetria
    grilla_netcdf = xr.Dataset(variables)
    grilla_netcdf["spherical"] = np.array(b'T', dtype='|S1')
    grilla_netcdf.attrs["Descripcion"] = "Grilla para MSOT / CROCO"
    grilla_netcdf.attrs["Autor"] = "grillado_msot.py"
    grilla_netcdf.attrs["Fecha de cracion"] = pd.Timestamp.now().isoformat()
    #################################################################################
    # Guardar
    #################################################################################
    grilla_netcdf.close()
    grilla_netcdf.fillna(1).to_netcdf(nombre)
    return
################################################################################
# Ejecucion de todas las funciones
################################################################################
if __name__ == "__main__":
    make_grd()

print(" \n")
print("Lista la grilla de ", Nlon, " x ",Nlat, "\n")
print("# En Run/params.h editar (Argentina o Argentina_Anidado1)\n")
print("# MMm0 = ",Nlat-2,"; LLm0 = ", Nlon-2)
