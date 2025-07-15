      subroutine get_psource
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
	      parameter (LLm0=240,   MMm0=359,  N=1)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=8)
      parameter (NSUB_X=2, NSUB_E=4)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Ntides
      parameter (Ntides=8)
      real D_wetdry
      parameter (D_wetdry=0.2D0)
      integer*4 Msrc
      parameter (Msrc=100)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=LLm, Mm=MMm)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=1)
      real sustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real svstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real svstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd
      integer*4 lsvsgrd,sms_tid, susid, svsid
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_cycle, sms_scale
      common /smsdat3/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat4/ lsvsgrd,sms_tid, susid, svsid
      integer*4 lwgrd, wid
      common /smsdat5/ lwgrd, wid
	real  Pair(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
	real  pair2d(0:Lm+1+padd_X,0:Mm+1+padd_E)
	integer pairid
	common /forces_pair/ Pair, pair2d
	common /pairdat/ pairid
	real  QBFC(0:Lm+1+padd_X,0:Mm+1+padd_E)
	integer qbfcid
	common /forces_qbfc/ QBFC
	common /qbfcdat/ qbfcid
      real bustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real bvstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real bvstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
      logical ldefhis
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_64bit_data
      integer*4 nf_cdf5
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_format_64bit_offset
      integer*4 nf_format_64bit_data
      integer*4 nf_format_cdf5
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_int64
      external        nf_put_att_int64
      integer*4         nf_get_att_int64
      external        nf_get_att_int64
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_compact
      parameter (nf_compact = 2)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_noquantize
      parameter (nf_noquantize = 0)
      integer*4 nf_quantize_bitgroom
      parameter (nf_quantize_bitgroom = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_estorage
      parameter (nf_estorage = -126)
      integer*4 nf_ebadchunk
      parameter (nf_ebadchunk = -127)
      integer*4 nf_enotbuilt
      parameter (nf_enotbuilt = -128)
      integer*4 nf_ediskless
      parameter (nf_ediskless = -129)
      integer*4 nf_ecantextend
      parameter (nf_ecantextend = -130)
      integer*4 nf_empi
      parameter (nf_empi = -131)
      integer*4 nf_efilter
      parameter (nf_efilter = -132)
      integer*4 nf_ercfile
      parameter (nf_ercfile = -133)
      integer*4 nf_enullpad
      parameter (nf_enullpad = -134)
      integer*4 nf_einmemory
      parameter (nf_einmemory = -135)
      integer*4 nf_enofilter
      parameter (nf_enofilter = -136)
      integer*4 nf_enczarr
      parameter (nf_enczarr = -137)
      integer*4 nf_es3
      parameter (nf_es3 = -138)
      integer*4 nf_eempty
      parameter (nf_eempty = -139)
      integer*4 nf_eobject
      parameter (nf_eobject = -140)
      integer*4 nf_enoobject
      parameter (nf_enoobject = -141)
      integer*4 nf_eplugin
      parameter (nf_eplugin = -142)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_zstandard
      external nf_def_var_zstandard
      integer*4 nf_inq_var_zstandard
      external nf_inq_var_zstandard
      integer*4 nf_def_var_szip
      external nf_def_var_szip
      integer*4 nf_inq_var_szip
      external nf_inq_var_szip
      integer*4 nf_def_var_quantize
      external nf_def_var_quantize
      integer*4 nf_inq_var_quantize
      external nf_inq_var_quantize
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_def_var_filter
      external nf_def_var_filter
      integer*4 nf_inq_var_filter
      external nf_inq_var_filter
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      integer*4 nf_set_log_level
      external nf_set_log_level
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxSSH
      parameter (indxSSH=indxVb+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 indxQBAR
      parameter (indxQBAR=indxSUSTR+80)
        integer*4 indxPAIR
        parameter (indxPAIR = indxSSH+12)
        integer*4 indxQBFC
        parameter (indxQBFC = indxSSH+21)
      integer*4 indxBhflx
      parameter (indxBhflx=indxSUSTR+131)
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisQBFC
      logical wrthis(500)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &      , hisQBFC
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , wrthis
      character*80 date_str, title, start_date
      character*80 origin_date, start_date_run
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8)             :: origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,   bry_file
      character*75  vname(20, 90)
      common /cncscrum/   date_str,   title,  start_date
     &         ,   origin_date, start_date_run
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,   bry_file
     &                                ,   vname
      real Qbar0(Msrc)
      common /sources_Qbar0/ Qbar0
      real Qbar(Msrc)
      common /sources_Qbar/ Qbar
      real Qsrc(Msrc,N)
      common /source_Qsrc/ Qsrc
      real Qshape(Msrc,N)
      common /source_Qshape/ Qshape
      real Tsrc(Msrc,N,1)
      common /source_Tsrc/ Tsrc
      real Tsrc0(Msrc,1)
      common /source_Tsrc0/ Tsrc0
      logical Lsrc(Msrc,30)
      common /source_Lsrc/ Lsrc
      real lasrc(Msrc)
      common /source_lasrc/ lasrc
      real losrc(Msrc)
      common /source_losrc/ losrc
      integer*4 Nsrc
      common /source_Nsrc/ Nsrc
      integer*4 Dsrc(Msrc)
      common /source_Dsrc/ Dsrc
      integer*4 Isrc(Msrc)
      common /source_Isrc/ Isrc
      integer*4 Jsrc(Msrc)
      common /source_Jsrc/ Jsrc
      real qbarg(Msrc,2)
      common /qbardat_qbarg/qbarg
      real    qbar_time(2)
      real    qbar_cycle
      integer*4 itqbar, qbar_ncycle, qbar_rec,  qbar_tid,  qbar_id
      common /qbardat1/ qbar_time
      common /qbardat2/ qbar_cycle
      common /qbardat3/ itqbar, qbar_ncycle, qbar_rec, qbar_tid, qbar_id
      real qbardir(Msrc)
      common /source_qbardir/ qbardir
      real cff
      integer*4 is,i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle,
     &        itrc
      integer*4 s(2),c(2)
      if (may_day_flag.ne.0) return
      if (itqbar.eq.0 .or. iic.eq.0) then
        lstr=lenstr(qbarname)
        if (ncidqbar.eq.-1) then
          ierr=nf_open(qbarname(1:lstr), nf_nowrite, ncidqbar)
          if (ierr. ne. nf_noerr) goto 4
        endif
        ierr=nf_inq_varid (ncidqbar, 'qbar_time', qbar_tid)
        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'qbar_time', qbarname(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxQBAR))
        ierr=nf_inq_varid (ncidqbar, vname(1,indxQBAR)(1:lvar), qbar_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidqbar, qbar_id, i)
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxQBAR)(1:lvar), qbarname(1:lstr)
          goto 99
        endif
         call set_cycle (ncidqbar,   qbar_tid,    ntqbar,
     &                  qbar_cycle, qbar_ncycle, qbar_rec)
        if (may_day_flag.ne.0) return
        itqbar=2
        qbar_time(1)=-1.D+20
        qbar_time(2)=-1.D+20
      endif
 1    i=3-itqbar
      cff=time+0.5D0*dt
      if (qbar_time(i).le.cff .and. cff.lt.qbar_time(itqbar)) return
      ierr=advance_cycle (qbar_cycle, ntqbar, qbar_ncycle, qbar_rec)
      if (ierr .ne. 0) goto 7
      ierr=nf_get_var1_double (ncidqbar, qbar_tid, qbar_rec, cff)
      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'qbar_time', qbar_rec
        goto 99
      endif
      qbar_time(i)=cff*day2sec+qbar_cycle*qbar_ncycle
      if (qbar_time(itqbar).eq.-1.D+20) qbar_time(itqbar)=qbar_time(i)
      s=(/ qbar_rec,1 /)
      c=(/ 1,Nsrc /)
      ierr=nf_get_vara_double(ncidqbar,qbar_id,s,c,qbarg(1:Nsrc,i))
      if (ierr .ne. nf_noerr) then
        write(stdout,6) 'QBAR', qbar_rec
        goto 99
      endif
      itqbar=i
       write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)')
     &     'GET_PSOURCE --',
     &     'Read Run-off flow fields      for time =', cff
      if (ntqbar.gt.1) goto 1
      if (ntqbar.eq.1) return
  3   format(/,' GET_PSOURCE - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in runoff NetCDF file: ',a)
  4   write(stdout,5) qbarname(1:lstr)
  5   format(/,' GET_PSOURCE - ERROR: unable to open runoff NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_PSOURCE - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) qbar_rec, ntqbar, qbarname(1:lstr), tdays,
     &                                      qbar_time(itqbar)*sec2day
  8   format(/,' GET_PSOURCE - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available',  /,  11x,'record ',I4,
     &       1x,'in runoff NetCDF file: ', a, /,   11x,'TDAYS = ',
     &                g12.4,2x,'last available QBAR_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end
      subroutine set_psource_tile(Istr,Iend,Jstr,Jend)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend,i, it1,it2
      real cff, cff1,cff2
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
	      parameter (LLm0=240,   MMm0=359,  N=1)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=8)
      parameter (NSUB_X=2, NSUB_E=4)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Ntides
      parameter (Ntides=8)
      real D_wetdry
      parameter (D_wetdry=0.2D0)
      integer*4 Msrc
      parameter (Msrc=100)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=LLm, Mm=MMm)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=1)
      real sustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real svstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real svstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd
      integer*4 lsvsgrd,sms_tid, susid, svsid
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_cycle, sms_scale
      common /smsdat3/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat4/ lsvsgrd,sms_tid, susid, svsid
      integer*4 lwgrd, wid
      common /smsdat5/ lwgrd, wid
	real  Pair(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
	real  pair2d(0:Lm+1+padd_X,0:Mm+1+padd_E)
	integer pairid
	common /forces_pair/ Pair, pair2d
	common /pairdat/ pairid
	real  QBFC(0:Lm+1+padd_X,0:Mm+1+padd_E)
	integer qbfcid
	common /forces_qbfc/ QBFC
	common /qbfcdat/ qbfcid
      real bustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real bvstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real bvstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
      logical ldefhis
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      real Qbar0(Msrc)
      common /sources_Qbar0/ Qbar0
      real Qbar(Msrc)
      common /sources_Qbar/ Qbar
      real Qsrc(Msrc,N)
      common /source_Qsrc/ Qsrc
      real Qshape(Msrc,N)
      common /source_Qshape/ Qshape
      real Tsrc(Msrc,N,1)
      common /source_Tsrc/ Tsrc
      real Tsrc0(Msrc,1)
      common /source_Tsrc0/ Tsrc0
      logical Lsrc(Msrc,30)
      common /source_Lsrc/ Lsrc
      real lasrc(Msrc)
      common /source_lasrc/ lasrc
      real losrc(Msrc)
      common /source_losrc/ losrc
      integer*4 Nsrc
      common /source_Nsrc/ Nsrc
      integer*4 Dsrc(Msrc)
      common /source_Dsrc/ Dsrc
      integer*4 Isrc(Msrc)
      common /source_Isrc/ Isrc
      integer*4 Jsrc(Msrc)
      common /source_Jsrc/ Jsrc
      real qbarg(Msrc,2)
      common /qbardat_qbarg/qbarg
      real    qbar_time(2)
      real    qbar_cycle
      integer*4 itqbar, qbar_ncycle, qbar_rec,  qbar_tid,  qbar_id
      common /qbardat1/ qbar_time
      common /qbardat2/ qbar_cycle
      common /qbardat3/ itqbar, qbar_ncycle, qbar_rec, qbar_tid, qbar_id
      real qbardir(Msrc)
      common /source_qbardir/ qbardir
      integer*4 IstrR,IendR,JstrR,JendR
      if (Istr.eq.1) then
        IstrR=Istr-1
      else
        IstrR=Istr
      endif
      if (Iend.eq.Lm) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (Jstr.eq.1) then
        JstrR=Jstr-1
      else
        JstrR=Jstr
      endif
      if (Jend.eq.Mm) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      it1=3-itqbar
      it2=itqbar
      cff=time+0.5D0*dt
      cff1=qbar_time(it2)-cff
      cff2=cff-qbar_time(it1)
      if (qbar_cycle.lt.0.D0) then
        if (iic.eq.ntstart) then
          do i=1,Nsrc
            Qbar(i)=qbarg(i,itqbar)
          enddo
        endif
      elseif (cff1.ge.0.D0 .and. cff2.ge.0.D0) then
        if (Istr+Jstr.eq.2 .and. cff1.lt.dt) synchro_flag=.TRUE.
        cff=1.D0/(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff
        do i=1,Nsrc
          Qbar(i)=cff1*qbarg(i,it1)+cff2*qbarg(i,it2)
        enddo
      else
        if (Istr+Jstr.eq.2) then
          write(stdout,1) 'qbar_time', tdays, qbar_time(it2)*sec2day
  1       format(/,' SET_PSOURCE - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
