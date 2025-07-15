      subroutine init_arrays (tile)
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
      integer*4 tile, trd
C$    integer*4 omp_get_thread_num
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lm+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lm)/2
      chunk_size_E=(Mm+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mm)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lm)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mm)
      call init_arrays_tile (Istr,Iend,Jstr,Jend)
      return
      end
      subroutine init_arrays_tile (Istr,Iend,Jstr,Jend)
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
      real h(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real hinv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real f(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real fomn(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_angler/angler
      real latr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real lonr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real latu(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real lonu(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real latv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real lonv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_latr/latr /grid_lonr/lonr
      common /grid_latu/latu /grid_lonu/lonu
      common /grid_latv/latv /grid_lonv/lonv
      real pm(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pm_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pm_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real dndx(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmon_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmon_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real grdscl(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real umask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vmask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmask2(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real rmask_wet(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmask_wet(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real umask_wet(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vmask_wet(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rmask_wet_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real Dcrit(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real wetdry(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /mask_r_wet/rmask_wet /mask_p_wet/pmask_wet
      common /mask_u_wet/umask_wet /mask_v_wet/vmask_wet
      common /mask_r_wet_avg/rmask_wet_avg
      common /Dcrit_wet/Dcrit
      common /wetdry_wet/wetdry
      real zob(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /Z0B_VAR/zob
      real zeta(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      real ubar(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      real vbar(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
      real urhs(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vrhs(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real Duon(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real DVom(0:Lm+1+padd_X,0:Mm+1+padd_E)
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
      real Tperiod(Ntides)
      common /tides_Tperiod/ Tperiod
      real SSH_Tamp(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_SSH_Tamp/ SSH_Tamp
      real SSH_Tphase(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_SSH_Tphase/ SSH_Tphase
      real zetabry_west(0:Mm+1+padd_E),
     &    zetabry_west_dt(0:Mm+1+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real zetabry_east(0:Mm+1+padd_E),
     &    zetabry_east_dt(0:Mm+1+padd_E,2)
      common /bry_zeta_east/ zetabry_east, zetabry_east_dt
      real zetabry_south(0:Lm+1+padd_X),
     &    zetabry_south_dt(0:Lm+1+padd_X,2)
      common /bry_zeta_south/ zetabry_south, zetabry_south_dt
      real zetabry_north(0:Lm+1+padd_X),
     &    zetabry_north_dt(0:Lm+1+padd_X,2)
      common /bry_zeta_north/ zetabry_north, zetabry_north_dt
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
      integer*4 Istr,Iend,Jstr,Jend, i,j
     &       ,itide
      real init
      parameter (init=0.D0)
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
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta(i,j,1)=0.D0
          zeta(i,j,2)=init
          zeta(i,j,3)=init
          zeta(i,j,4)=init
          ubar(i,j,1)=init
          ubar(i,j,2)=init
          ubar(i,j,3)=init
          ubar(i,j,4)=init
          vbar(i,j,1)=init
          vbar(i,j,2)=init
          vbar(i,j,3)=init
          vbar(i,j,4)=init
          rmask_wet(i,j)=1
          umask_wet(i,j)=1
          vmask_wet(i,j)=1
          pmask_wet(i,j)=1
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          sustr(i,j)=init
          svstr(i,j)=init
          sustrg(i,j,1)=init
          svstrg(i,j,1)=init
          sustrg(i,j,2)=init
          svstrg(i,j,2)=init
          bustr(i,j)=init
          bvstr(i,j)=init
          bustrg(i,j,1)=init
          bvstrg(i,j,1)=init
          bustrg(i,j,2)=init
          bvstrg(i,j,2)=init
        enddo
      enddo
        do i=1,Msrc
                qbar(i)=init
                qbarg(i,1)=init
                qbarg(i,2)=init
        enddo
      if (istr.eq.1) then
        do j=JstrR,JendR
          zetabry_west(j)=init
        enddo
      endif
      if (iend.eq.Lm) then
        do j=JstrR,JendR
          zetabry_east(j)=init
        enddo
      endif
      if (jstr.eq.1) then
        do i=IstrR,IendR
          zetabry_south(i)=init
        enddo
      endif
      if (jend.eq.Mm) then
        do i=IstrR,IendR
          zetabry_north(i)=init
        enddo
      endif
      do itide=1,Ntides
        Tperiod(itide)=init
        do j=JstrR,JendR
          do i=IstrR,IendR
            SSH_Tamp(i,j,itide)=init
            SSH_Tphase(i,j,itide)=init
          enddo
        enddo
      enddo
      return
      end
