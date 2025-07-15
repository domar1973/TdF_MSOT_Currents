#!/usr/bin/env python
### Versión sin viento

################################################################################
# Liberias
################################################################################
import xarray as xr                     # Uso de campos
import pandas as pd                     # Uso de dataframes y fechas
import numpy as np                      # Operaciones matematicas
import uptide as tides                  # Analisis de marea, correciones nodales
import os                               # Instrucciones a linux
import glob                             # Lectura de archivos en el SO
import re                               # Regular expressions
from tqdm import tqdm
from scipy.interpolate import griddata
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
################################################################################

################################################################################
# Definciones generales
################################################################################
formato_netcdf = "NETCDF4"
nombre = os.environ["FECHA_INIT"]
fecha_inicial = pd.to_datetime(nombre, format = "%Y%m%d_%H")
fecha_final = os.environ["FECHA_FIN"]
fecha_final = pd.to_datetime(fecha_final, format = "%Y%m%d_%H")
miembro = int(os.environ["MIEM"])
path_salida = os.environ["DIR_OUTPUT"]
delta_t = int(os.environ["DT"])
dominio = os.environ["DOMINIO"]
archivo_marea = os.environ["MAREA"]
rundir= os.environ["RUNDIR"]
rank = int(os.environ["RANK"])

################################################################################
# Archivo con la deficion de inputs
################################################################################
def make_input(miembro = miembro, nombre = nombre, fecha_inicial = fecha_inicial, fecha_final = fecha_final,
               dominio = dominio, path_salida = path_salida, rundir=rundir, delta_t=delta_t, rank=rank):
    # Nombre del archivo
    nombre_input = "{}/{}_{:02d}.in".format(rundir,nombre, miembro)
    # Habilito la escritura del archivo
    txt = open(nombre_input, "w")
    # Nombre de la corrida
    txt.write("title:\n")
    txt.write(nombre + "_{:02d} \n".format(miembro))
    # Fecha de la simulacion
    txt.write("start_date:\n")
    txt.write(fecha_inicial.strftime("%Y-%b-%d %H:%M:%S") + "\n")
    txt.write("time_stepping: NTIMES   dt[sec]  NDTFAST  NINFO\n")
    # Pasos necesarios para correr el model
    ntimes   = (fecha_final - fecha_inicial).total_seconds()//delta_t
    txt.write("{} {} 1 {} \n".format(int(ntimes), 
                                     int(delta_t), int(21600/delta_t)))
    # Dominio
    txt.write("grid:  filename\n")
    txt.write(dominio + "\n")
    # Archivo forzante
    txt.write("forcing: filename\n")
    txt.write("{}{}_m{:02d}_frc.nc \n".format(path_salida, nombre, miembro))
    # Condiciones de contorno
    if rank != 0:
        txt.write("boundary: filename\n")
        txt.write("{}{}_r{:02d}_m{:02d}_bry.nc\n".format(path_salida, nombre, rank, miembro))
    # Archivo con condicion inicial
    txt.write("initial: NRREC / filename\n")
    txt.write("0\n")
    txt.write("{}{}_m{:02d}_ini.nc \n".format(path_salida, nombre, miembro))
    # Archivo para el restart - NO SE USA PERO TIENE QUE ESTAR
    txt.write("restart:          NRST, NRPFRST / filename\n")
    txt.write("{} 0 \n".format(int(86400//delta_t)))
    txt.write(path_salida + "croco_rst.nc\n")
    # Archivo con la solucion del modelo
    txt.write("history: LDEFHIS, NWRT, NRPFHIS / filename\n")
    txt.write("T {} 0 \n".format(int(3600//delta_t)))
    if rank == 0:
        txt.write("{}{}_m{:02d}_his.nc \n".format(path_salida, nombre, miembro))
    else:
        txt.write("{}{}_r{:02d}_m{:02d}_his.nc \n".format(path_salida, nombre, rank, miembro))
    txt.write("primary_history_fields: zeta UBAR VBAR QBFC  U  V  wrtT(1:NT)\n")
    txt.write("T    T   T T  F  F    30*F\n")
    # Parametros estaticos de la simulacion
    txt.write("rho0:\n")
    txt.write("1025 \n")
    txt.write("lateral_visc:   VISC2,    VISC4    [m^2/sec for all] \n")
    txt.write("0.       0.\n")
    txt.write("bottom_drag:     RDRG [m/s],  RDRG2,  Zob [m],  Cdb_min, Cdb_max\n")
    txt.write("0.0 2.25d-3 1.d-3     1.d-4    1.d-1\n")
    txt.write("gamma2:\n")
    txt.write("1.\n")
    txt.write("nudg_cof:    TauT_in, TauT_out, TauM_in, TauM_out  [days for all]\n")
    txt.write("1.       360.      3.      360.\n")
    # Archivo con la descarga climatologica de los 8 rio que hay datos
    txt.write("psource_ncfile:   Nsrc  Isrc  Jsrc  Dsrc qbardir  Lsrc  Tsrc   runoff file name\n")
    txt.write("{}{}_m{:02d}_runoff.nc \n".format(path_salida, nombre, miembro))
    txt.write("8 \n")
    txt.write("8 429 0 1 30*T 15. 0.\n")
    txt.write("5 445 0 1 30*T 15. 0.\n")
    txt.write("4 446 0 1 30*T 15. 0.\n")
    txt.write("29 265 1 -1 30*T 15. 0.\n")
    txt.write("24 293 0 1 30*T 15. 0.\n")
    txt.write("40 405 0 1 30*T 15. 0.\n")
    txt.write("11 108 1 -1 30*T 15. 0.\n")
    txt.write("4 88 0 1 30*T 15. 0.\n")
    ### # DOB: AGREGO linea runoff vacía
    ### txt.write("psource_ncfile:\n")
    ### # \DOB
    txt.close()
    return

################################################################################
# Archivo con condicion inicial
################################################################################
def make_init(miembro = miembro, dominio = dominio, fecha_inicial = fecha_inicial,
              nombre = nombre, path_salida = path_salida, formato_netcdf = formato_netcdf, rank=rank):
    # Cargo el dominio como archivo de referencia
    grd = xr.open_dataset(dominio)
    # Cargo las variables generales que pide CROCO
    variables = {"spherical": ("one", [1]), "Vtransform": ("one", [1]),
                 "Vstretching": ("one", [1]), "theta_s": ("one", [6.]),
                 "theta_b": ("one", [0.]), "Tcline": ("one", [0.5]),
                 "hc": ("one", [0.5]), "sc_r": ("one", [-0.5]),
                 "Cs_r": ("one", [-0.0497]), "tstart": ("one", [0]),
                 "tend": ("one", [0]), "ocean_time": ("one", [0]),
                 "scrum_time": ("one", [0]),}
    # Cargo la ultima solucion, sino existe inicializo con ceros
    nombre_previo = (fecha_inicial - pd.DateOffset(days = 0)).strftime("%Y%m%d_%H") + "_prono.nc"
    if rank == 0:
        if len(glob.glob(nombre_previo)) != 0:
            print(" <<< Inicializando con las soluciones del modelo: {} >>> ".format(nombre_previo))       
            # abro archivo que contiene la fecha de spin up
            inicio    = xr.open_dataset(nombre_previo).sel(miembro = miembro)
            # fecha del inicio del spin-up
            finicio = os.path.basename(nombre_previo).split("_")[0] +"_"+os.path.basename(nombre_previo).split("_")[1]
            estado_inicial = inicio.sel(time = pd.to_datetime(finicio, format = "%Y%m%d_%H"))
            print("condicion inicial")
            print(estado_inicial.time)
            variables.update({"ubar": (["time", "eta_u", "xi_u"], estado_inicial.ubar.fillna(0).values[np.newaxis, :, :]),
                              "vbar": (["time", "eta_v", "xi_v"], estado_inicial.vbar.fillna(0).values[np.newaxis, :, :]),
                              "zeta": (["time", "eta_rho", "xi_rho"], estado_inicial.zeta.fillna(0).values[np.newaxis, :, :]),})
        else:
            print(" <<< Inicializando con ceros >>> ")
            variables.update({"ubar": (["time", "eta_u", "xi_u"], 0 * grd.mask_u.values[np.newaxis, :, :]),
                              "vbar": (["time", "eta_v", "xi_v"], 0 * grd.mask_v.values[np.newaxis, :, :]),
                              "zeta": (["time", "eta_rho", "xi_rho"], 0 * grd.mask_rho.values[np.newaxis, :, :]),})
    else:
        if rank == 1:
            SolucionPadre_nombrearch = "../"+"{}{}_m{:02d}_his.nc".format(path_salida, nombre, miembro)
            SolucionPadre = xr.open_dataset(SolucionPadre_nombrearch).isel(time = 0)
        else:
            SolucionPadre = xr.open_dataset("{}{}_r{:02d}_m{:02d}_his.nc".format(path_salida, nombre, rank - 1, miembro)).isel(time = 0)
        # ZETA
        Coordenadas_hijo = grd[["lat_rho", "lon_rho"]].to_dataframe()
        rangos = [range(grd.eta_rho.data[0]), range(grd.xi_rho.data[0])]
        zeta0 = SolucionPadre.zeta.bfill("xi_rho").to_dataframe()
        campo_zeta = griddata(zeta0[["lat_rho", "lon_rho"]].values, 
                              zeta0["zeta"].values, 
                              Coordenadas_hijo.values, 
                              method = 'linear')
        campo_zeta = pd.DataFrame(data = campo_zeta, columns = ["zeta"], index = Coordenadas_hijo.index).to_xarray().zeta.fillna(0).data
        # UBAR
        Coordenadas_hijo = grd[["lat_u", "lon_u"]].to_dataframe()
        rangos = [range(grd.eta_rho.data[0]), range(grd.xi_u.data[0])]
        ubar0 = SolucionPadre.ubar.bfill("xi_u").to_dataframe()
        campo_ubar = griddata(ubar0[["lat_u", "lon_u"]].values, 
                              ubar0["ubar"].values, 
                              Coordenadas_hijo.values, 
                              method = 'linear')
        campo_ubar = pd.DataFrame(data = campo_ubar, columns = ["ubar"], index = Coordenadas_hijo.index).to_xarray().ubar.fillna(0).data
        # VBAR
        Coordenadas_hijo = grd[["lat_v", "lon_v"]].to_dataframe()
        rangos = [range(grd.eta_v.data[0]), range(grd.xi_rho.data[0])]
        vbar0 = SolucionPadre.vbar.bfill("xi_rho").to_dataframe()
        campo_vbar = griddata(vbar0[["lat_v", "lon_v"]].values, 
                              vbar0["vbar"].values, 
                              Coordenadas_hijo.values, 
                              method = 'linear')
        campo_vbar = pd.DataFrame(data = campo_vbar, columns = ["vbar"], index = Coordenadas_hijo.index).to_xarray().vbar.fillna(0).data
        # Guardo
        variables.update({"ubar": (["time", "eta_u", "xi_u"], campo_ubar[np.newaxis, :, :]),
                          "vbar": (["time", "eta_v", "xi_v"], campo_vbar[np.newaxis, :, :]),
                          "zeta": (["time", "eta_rho", "xi_rho"], campo_zeta[np.newaxis, :, :]),})
    # Creo el netcdf y le agrego atributos
    ini = xr.Dataset(variables)
    ini["tstart"].attrs = {"long_name": "start processing day","units": "day"}
    ini["tend"].attrs = {"long_name": "end processing day","units": "day"}
    ini["ocean_time"].attrs = {"long_name": "time since initialization","units": "second"}
    ini["scrum_time"].attrs = {"long_name": "time since initialization","units": "second"}
    ini.to_netcdf("{}{}_m{:02d}_ini.nc".format(path_salida, nombre, miembro), format = formato_netcdf)
    ini.close()  
    return

################################################################################
# Archivo con descarga de rios
################################################################################
def make_runoff(miembro = miembro, dominio = dominio, fecha_inicial = fecha_inicial,
                nombre = nombre, path_salida = path_salida, formato_netcdf = formato_netcdf, rank=rank):
    # Descargas climatologicas del Rio de la Plata
    descarga_mensual_climatologica = {1 : 22668.0, 2 : 24295.0, 3 : 27650.0, 
                                      4 : 29190.0, 5 : 30327.0, 6 : 29718.0, 
                                      7 : 28614.0, 8 : 25428.0, 9 : 22380.0, 
                                      10: 25140.0, 11: 26903.0, 12: 24324.0,}
    # Creo el vector descargas para los 10 dias de simulacion
    Qbar = descarga_mensual_climatologica[fecha_inicial.month] * np.ones((10, 3))
    # Proporcion de los afluentes
    Qbar[:,0] *= 0.24
    Qbar[:,1] *= 0.56
    Qbar[:,2] *= 0.20
    # Rios fuera del Rio de la Plata
    climato_rios = {1 :  ( 699, 207, 28, 27, 11),
                    2 :  ( 461, 131, 55, 25, 11),
                    3 :  ( 398,  91, 34, 22, 13),
                    4 :  ( 370,  76, 45, 18, 14),
                    5 :  ( 548,  83, 61, 21, 16),
                    6 :  ( 983,  91, 74, 21, 21),
                    7 :  (1121,  88, 86, 21, 17),
                    8 :  (1159,  81, 97, 21, 32),
                    9 :  (1028,  84, 71, 22, 34),
                    10:  (1136, 130, 46, 24, 37),
                    11:  (1261, 235, 33, 24, 23),
                    12:  (1073, 268, 19, 30, 15),}
    descarga_mensual = np.repeat([climato_rios[fecha_inicial.month]], repeats = Qbar.shape[0], axis = 0)
    # Agrego la descarga de los otros rios
    Qbar = np.concatenate((Qbar, descarga_mensual), axis = 1)
    # Tiempo de referencia para comenzar el con valor inicial de la descarga, consideracion una varacioan ciclica - NO SE USA
    cycle_length = 10950.
    # Time vector
    qbar_time = np.arange(len(Qbar))
	# Creo netCDF
    posiciones_de_rios = np.array([[8,429,0,1], [5,445,0, 1], [4,446,0,1], [29,265, 1, -1], 
                                  [24, 293, 0, 1], [40, 405, 0, 1], [11, 108, 1, -1], [4, 88, 0, 1]])
    # Dependiendo del rank
    if rank != 0:
        Qbar *= 0
    variables = {"qbar_time": ("qbar_time", qbar_time),
                 "runoff_position": (["n_qbar", "two"], posiciones_de_rios[:,:2]),
                 "runoff_direction": (["n_qbar", "two"], posiciones_de_rios[:,2:]),
                 "Qbar": (["n_qbar", "qbar_time"], Qbar.T) }
    rio = xr.Dataset(variables)
    rio.qbar_time.attrs = {"long_name": "runoff time", "units":"days", "cycle_length": cycle_length}
    rio.runoff_position.attrs = {"long_name": "position of the runoff (by line) in the ROMS grid"}
    rio.runoff_direction.attrs = {"long_name": "direction/sense of the runoff (by line) in the ROMS grid"}
    rio.Qbar.attrs = {"long_name": "runoff discharge", "units": "m3 s-1"}
    rio.to_netcdf("{}{}_m{:02d}_runoff.nc".format(path_salida, nombre, miembro), format = formato_netcdf)
    rio.close()
    return

################################################################################
# Archivo con el forzante atmosferico y de mareas
################################################################################
def make_forcing(miembro = miembro, dominio = dominio, fecha_inicial = fecha_inicial, fecha_final = fecha_final,
                 nombre = nombre, path_salida = path_salida, formato_netcdf = formato_netcdf,
                 archivo_marea = archivo_marea, rank=rank):

    # Cargo el archivo forzante 
    # Todas las lineas con frc_atmo comentadas por DOB ("### ") para "apagar el viento"
    ### frc_atmo = xr.open_dataset("vientos_prueba_enero-febrero_2019.nc")#.sel(number = miembro)
    ### frc_atmo = frc_atmo.sel(time = slice(fecha_inicial, fecha_final))
    # Interpolacion al dominio
    grd = xr.open_dataset(dominio)
    # Esto depende si la longitud del forzante esta en grados Oestes o Estes
    # frc_atmo = frc_atmo.assign_coords(longitude = frc_atmo.longitude - 360)
    ### frc_atmo_interp = frc_atmo.interp(longitude = grd.lon_rho, latitude = grd.lat_rho, method = "linear")
    ### frc_atmo_interp = frc_atmo_interp.bfill("xi_rho")
    # Variables
    ### w10 = np.hypot(frc_atmo_interp.u10, frc_atmo_interp.v10)
    # Calibracion empirica de Dinapoli et al (2022b)
    ### cd = -2.0144e-3 + 4.0e-4 * w10
    ### cd = cd.where(7.45 < w10, 0.9656e-3)
    ### sustr = (cd * 1.20 * w10 * frc_atmo_interp.u10).values
    ### svstr = (cd * 1.20 * w10 * frc_atmo_interp.v10).values
    ### frc_atmo_dt = (frc_atmo.time.diff(dim = 'time')//(3600e9)).astype('float').values
    ### frc_atmo_dt = np.append(frc_atmo_dt, frc_atmo_dt[-1])
    ### vector_tiempo = frc_atmo_dt.cumsum() - frc_atmo_dt/2
    ### vector_tiempo[0] *= -1
    ### variables = {"sms_time": ("sms_time", vector_tiempo/24.),
    ###              "sustr": (["sms_time", "eta_u", "xi_u"], 0.5 * (sustr[:,:,0:-1] + sustr[:,:,1:])),
    ###              "svstr": (["sms_time", "eta_v", "xi_v"], 0.5 * (svstr[:,0:-1,:] + svstr[:,1:,:])),
    ###              "Pair": (["sms_time", "eta_rho", "xi_rho"], frc_atmo_interp.sp.values)}
    ### # Variables con viento apagado:
    variables = {"sms_time": ("sms_time", [0,1]),
                 "sustr": (["sms_time", "eta_u", "xi_u"],    [np.zeros_like(grd.mask_u)  ] * 2),
                 "svstr": (["sms_time", "eta_v", "xi_v"],    [np.zeros_like(grd.mask_v)  ] * 2),
                 "Pair": (["sms_time", "eta_rho", "xi_rho"], [np.zeros_like(grd.mask_rho)] * 2)}
    # Creo el archivo forzante
    frc = xr.Dataset(variables)
    # Attributes
    frc.sms_time.attrs = {"long_name": 'surface momentum stress time',
                          "units": 'days', "cycle_length": float(10950),}
    frc.sustr.attrs = {"long_name": 'surface u-momentum stress',
                       "units": 'Newton meter-2',}
    frc.svstr.attrs = {"long_name": 'surface v-momentum stress',
                       "units": 'Newton meter-2',}
    frc.Pair.attrs = {"long_name": 'sea level pressure',
                       "units": 'Pascal',}

    # Tides
    if rank == 0:
        tide_file     = xr.open_dataset(archivo_marea)
        if "mm" in tide_file.dims:
            tide_file     = tide_file.sel(mm = miembro % 4).squeeze()
        nodal_factors = tides.Tides(list(tide_file.nc.values))
        nodal_factors.set_initial_time(fecha_inicial)
        # el archivo de marea presenta la misma dimension que el dominio padre
        # No hace falta realizar interpolaciones
        Eamp = tide_file.ha
        Ephase = tide_file.hp * np.pi / 180

        for p, cmd in enumerate(Eamp["nc"]):
            Eamp.loc[cmd] *= nodal_factors.f[p]
            Ephase.loc[cmd] = np.mod((Ephase.loc[cmd] - nodal_factors.phi[p] - nodal_factors.u[p]) *  180/np.pi, 360)

        frc["tide_period"] = ('tide_period', 2*np.pi/nodal_factors.omega/3600)
        frc.tide_period.attrs = {"long_name": 'Tide angular period',
                                 "units": 'hours'}
        frc["tide_Ephase"] = (['tide_period','eta_rho','xi_rho'], Ephase.data)
        frc.tide_Ephase.attrs = {"lon_name": 'Tidal elevation phase angle',
                                 "units": 'degrees'}
        frc["tide_Eamp"] = (['tide_period','eta_rho','xi_rho'], Eamp.data)
        frc.tide_Eamp.attrs = {"lon_name": 'Tidal elevation amplitude',
                "units": 'meters'}    

    # Friccion de fondo
    frc["QBFC"] = (['eta_rho','xi_rho'], 2.25e-3 * np.ones_like(grd.h.values))
    frc.QBFC.attrs = {"lon_name": 'Quadratic bottom friction',
                             "units": ''}
    # Save netCDF
    frc.to_netcdf("{}{}_m{:02d}_frc.nc".format(path_salida, nombre, miembro), format = formato_netcdf)
    return

################################################################################
# Condiciones de contorno
################################################################################
def make_boundary(miembro = miembro, dominio = dominio, fecha_inicial = fecha_inicial,
                 nombre = nombre, path_salida = path_salida, formato_netcdf = formato_netcdf,
                 archivo_marea = archivo_marea, rank=rank):
    # Cargo solucion del dominio anterior
    SolucionPadre_nombrearch = "../"+"{}{}_m{:02d}_his.nc".format(path_salida, nombre, miembro)
    SolucionPadre = xr.open_dataset(SolucionPadre_nombrearch)
    # Diccionario con definiciones necesaria para crear el netCDF de boundary
    variables = {"spherical": ("one", [1]), "Vtransform": ("one", [1]),
                 "Vstretching": ("one", [1]), "theta_s": ("one", [6.]),
                 "theta_b": ("one", [0.]), "Tcline": ("one", [0.5]),
                 "hc": ("one", [0.5]), "sc_r": ("one", [-0.5]),
                 "sc_w": ("two", [-1, 0]), "Cs_r": ("two", [-1, 0]),
                 "Cs_r": ("one", [-0.0497]), "tstart": ("one", [0]),
                 "tend": ("one", [SolucionPadre.time.shape[0]/4]),
                 "bry_time": ("bry_time", np.arange(SolucionPadre.time.shape[0])/(24.)), # Los datos de tiempo van en formato días
                 "zeta_time": ("zeta_time", np.arange(SolucionPadre.time.shape[0])/(24.)),
                 "v2d_time": ("v2d_time", np.arange(SolucionPadre.time.shape[0])/(24.))}
    # Creo el archivo de contorno
    bry = xr.Dataset(variables)
    BordeAbierto = xr.open_dataset(dominio)
    # ZETA
    Coordenadas_hijo = BordeAbierto[["lat_rho", "lon_rho"]].to_dataframe()
    rangos = [range(BordeAbierto.eta_rho.data[0]), range(BordeAbierto.xi_rho.data[0])]
    zeta = [None] * SolucionPadre.time.size
    for it, t in enumerate(tqdm(SolucionPadre.time)):
        zeta0 = SolucionPadre.sel(time = t).zeta.bfill("xi_rho").to_dataframe()
        campo = griddata(zeta0[["lat_rho", "lon_rho"]].values, 
                         zeta0["zeta"].values, 
                         Coordenadas_hijo.values, 
                         method = 'linear')
        zeta[it] = pd.DataFrame(data = campo, columns = ["zeta"], index = Coordenadas_hijo.index).to_xarray().zeta
        
    # UBAR
    Coordenadas_hijo = BordeAbierto[["lat_u", "lon_u"]].to_dataframe()
    rangos = [range(BordeAbierto.eta_rho.data[0]), range(BordeAbierto.xi_u.data[0])]
    ubar = [None] * SolucionPadre.time.size
    for it, t in enumerate(tqdm(SolucionPadre.time)):
        ubar0 = SolucionPadre.sel(time = t).ubar.bfill("xi_u").to_dataframe()
        ubar0 = ubar0[["lon_u", "lat_u", "ubar"]].dropna()
        campo = griddata(ubar0[["lat_u", "lon_u"]].values, 
                         ubar0["ubar"].values, 
                         Coordenadas_hijo.values, 
                         method = 'linear')
        ubar[it] = pd.DataFrame(data = campo, columns = ["ubar"], index = Coordenadas_hijo.index).to_xarray().ubar

    # VBAR
    Coordenadas_hijo = BordeAbierto[["lat_v", "lon_v"]].to_dataframe()
    rangos = [range(BordeAbierto.eta_v.data[0]), range(BordeAbierto.xi_rho.data[0])]
    vbar = [None] * SolucionPadre.time.size
    for it, t in enumerate(tqdm(SolucionPadre.time)):
        vbar0 = SolucionPadre.sel(time = t).vbar.bfill("xi_rho").to_dataframe()
        vbar0 = vbar0[["lon_v", "lat_v", "vbar"]].dropna()
        campo = griddata(vbar0[["lat_v", "lon_v"]].values, 
                         vbar0["vbar"].values, 
                         Coordenadas_hijo.values, 
                         method = 'linear')
        vbar[it] = pd.DataFrame(data = campo, columns = ["vbar"], index = Coordenadas_hijo.index).to_xarray().vbar

    ######################################################################
    # Compilo
    ######################################################################
    ZETA = xr.concat(zeta, dim = "time").assign_coords(time = SolucionPadre.time).bfill("xi_rho")
    UBAR = xr.concat(ubar, dim = "time").assign_coords(time = SolucionPadre.time).bfill("xi_u")
    VBAR = xr.concat(vbar, dim = "time").assign_coords(time = SolucionPadre.time).bfill("xi_v")

    ######################################################################
    # Cargo los datos
    ######################################################################
    for var in ["zeta", "ubar", "vbar"]:   
        if "ubar" in var:
            bry[var + '_north'] = (["v2d_time", 'north' + "_" + var[0]], UBAR[:,-1].data)
            bry[var + '_east'] = (["v2d_time", 'east' + "_" + var[0]], UBAR[:,:,-1].data)
            bry[var + '_south'] = (["v2d_time", 'south' + "_" + var[0]], UBAR[:,0].data)
            bry[var + '_west'] = (["v2d_time", 'west' + "_" + var[0]], UBAR[:,:,0].data)
        elif "vbar" in var:
            bry[var + '_north'] = (["v2d_time", 'north' + "_" + var[0]], VBAR[:,-1].data)
            bry[var + '_east'] = (["v2d_time", 'east' + "_" + var[0]], VBAR[:,:,-1].data)
            bry[var + '_south'] = (["v2d_time", 'south' + "_" + var[0]], VBAR[:,0].data)
            bry[var + '_west'] = (["v2d_time", 'west' + "_" + var[0]], VBAR[:,:,0].data)
        else:
            bry[var + '_north'] = (["zeta_time", "north_rho"], ZETA[:,-1].data)
            bry[var + '_east'] = (["zeta_time", "east_rho"], ZETA[:,:,-1].data)
            bry[var + '_south'] = (["zeta_time", "south_rho"], ZETA[:,0].data)
            bry[var + '_west'] = (["zeta_time", "west_rho"], ZETA[:,:,0].data)

    ######################################################################
    # Guardo
    ######################################################################
    bry.to_netcdf("{}{}_r{:02d}_m{:02d}_bry.nc".format(path_salida, nombre, rank, miembro))
    bry.close()
    del bry

################################################################################
# Ejecucion de todas las funciones
################################################################################
if __name__ == "__main__":
    make_input()
    make_init()
    make_runoff()
    make_forcing()
    if rank != 0:
        make_boundary()
