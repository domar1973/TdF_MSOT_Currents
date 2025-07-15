      subroutine init_scalars (ierr)
      implicit none
      integer*4 ierr, i,j,itrc, lvar, lenstr
      character*20 nametrc, unitt
      character*60 vname1, vname3
      integer*4 omp_get_num_threads
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
	      parameter (LLm0=360,   MMm0=359,  N=1)
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
C$OMP PARALLEL
C$OMP CRITICAL (isca_cr_rgn)
      numthreads=omp_get_num_threads()
C$OMP END CRITICAL (isca_cr_rgn)
C$OMP END PARALLEL
       write(stdout,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
      if (numthreads.gt.NPP) then
         write(stdout,'(/1x,A,I3/)')
     &    'ERROR: Requested number of threads exceeds NPP =', NPP
        ierr=ierr+1
      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
         write(stdout,
     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
     &                'wrong choice of numthreads =', numthreads,
     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
        ierr=ierr+1
      endif
      time=0.D0
      tdays=0.D0
      PREDICTOR_2D_STEP=.FALSE.
      iic=0
      kstp=1
      krhs=1
      knew=1
      ntstart=1
      nfast=1
      synchro_flag=.true.
      first_time=0
      may_day_flag=0
      do j=0,NPP
        do i=0,31
          proc(i,j)=0
          CPU_time(i,j)=0.D0
        enddo
      enddo
      trd_count=0
      nrecrst=0
      nrechis=0
      tile_count=0
      bc_count=0
      avgke=0.D0
      avgpe=0.D0
      avgkp=0.D0
      volume=0.D0
      hmin=+1.D+20
      hmax=-1.D+20
      grdmin=+1.D+20
      grdmax=-1.D+20
      Cu_min=+1.D+20
      Cu_max=-1.D+20
      lonmin=+1.D+20
      lonmax=-1.D+20
      latmin=+1.D+20
      latmax=-1.D+20
      bc_crss=0.D0
      ncidrst=-1
      ncidhis=-1
      ncidfrc=-1
      ncidbulk=-1
      ncidclm=-1
      ncidqbar=-1
      ncidbtf=-1
      call get_date (date_str)
      vname(:,:)='  '
      vname(1,indxTime)='scrum_time'
      vname(2,indxTime)='time since initialization'
      vname(3,indxTime)='second'
      vname(4,indxTime)='time, scalar, series'
      vname(5,indxTime)='time'
      vname(6,indxTime)='    '
      vname(7,indxTime)='T'
      vname(8,indxTime)='.NO.'
      vname(9,indxTime)=' '
      vname(10,indxTime)=start_date
      vname(1,indxTime2)='time'
      vname(2,indxTime2)='time since initialization'
      vname(3,indxTime2)='second'
      vname(4,indxTime2)='time, scalar, series'
      vname(5,indxTime2)='time'
      vname(6,indxTime2)='  '
      vname(7,indxTime2)='T '
      vname(8,indxTime2)='.NO.'
      vname(9,indxTime2)=' '
      vname(10,indxTime2)=start_date
      vname(1,indxZ)='zeta                                        '
      vname(2,indxZ)='free-surface                                '
      vname(3,indxZ)='meter                                       '
      vname(4,indxZ)='free-surface, scalar, series                '
      vname(5,indxZ)='sea_surface_height                          '
      vname(6,indxZ)='lat_rho lon_rho                             '
      vname(7,indxZ)='                                            '
      vname(1,indxUb)='ubar                                       '
      vname(2,indxUb)='vertically integrated u-momentum component '
      vname(3,indxUb)='meter second-1                             '
      vname(4,indxUb)='ubar-velocity, scalar, series              '
      vname(5,indxUb)='barotropic_sea_water_x_'
     &                                  // 'velocity_at_u_location'
      vname(6,indxUb)='lat_u lon_u                                '
      vname(7,indxUb)='                                           '
      vname(1,indxVb)='vbar                                       '
      vname(2,indxVb)='vertically integrated v-momentum component '
      vname(3,indxVb)='meter second-1                             '
      vname(4,indxVb)='vbar-velocity, scalar, series              '
      vname(5,indxVb)='barotropic_sea_water_y_'
     &                                  // 'velocity_at_v_location'
      vname(6,indxVb)='lat_v lon_v                                '
      vname(7,indxVb)='                                           '
      vname(11,indxVb)=' '
      vname(1,indxBostr)='bostr                                   '
      vname(2,indxBostr)='Kinematic bottom stress                 '
      vname(3,indxBostr)='N/m2                                    '
      vname(4,indxBostr)='                                        '
      vname(5,indxBostr)='                                        '
      vname(6,indxBostr)='lat_rho lon_rho                         '
      vname(7,indxBostr)='                                        '
      vname(1,indxBustr)='bustr                                   '
      vname(2,indxBustr)='Kinematic u bottom stress component     '
      vname(3,indxBustr)='N/m2                                    '
      vname(4,indxBustr)='                                        '
      vname(5,indxBustr)='                                        '
      vname(6,indxBustr)='lat_u lon_u                             '
      vname(7,indxBustr)='                                        '
      vname(1,indxBvstr)='bvstr                                   '
      vname(2,indxBvstr)='Kinematic v bottom stress component     '
      vname(3,indxBvstr)='N/m2                                    '
      vname(4,indxBvstr)='                                        '
      vname(5,indxBvstr)='                                        '
      vname(6,indxBvstr)='lat_v lon_v                             '
      vname(7,indxBvstr)='                                        '
      vname(1,indxWstr)='wstr                                     '
      vname(2,indxWstr)='Kinematic wind stress                    '
      vname(3,indxWstr)='N/m2                                     '
      vname(4,indxWstr)='                                         '
      vname(5,indxWstr)='magnitude_of_surface_downward_stress     '
      vname(6,indxWstr)='lat_rho lon_rho                          '
      vname(7,indxWstr)='                                         '
      vname(1,indxUWstr)='sustr                                   '
      vname(2,indxUWstr)='Kinematic u wind stress component       '
      vname(3,indxUWstr)='N/m2                                    '
      vname(4,indxUWstr)='                                        '
      vname(5,indxUWstr)='surface_downward_eastward_stress        '
      vname(6,indxUWstr)='lat_u lon_u                             '
      vname(7,indxUWstr)='                                        '
      vname(8,indxUWstr)='                                        '
      vname(1,indxVWstr)='svstr                                   '
      vname(2,indxVWstr)='Kinematic v wind stress component       '
      vname(3,indxVWstr)='N/m2                                    '
      vname(4,indxVWstr)='                                        '
      vname(5,indxVWstr)='surface_downward_northward_stress       '
      vname(6,indxVWstr)='lat_v lon_v                             '
      vname(7,indxVWstr)='                                        '
      vname(1,indxPAIR)='Pair                                   '
      vname(2,indxPAIR)='Surfece pressure       '
      vname(3,indxPAIR)='N/m2                                    '
      vname(4,indxPAIR)='                                        '
      vname(5,indxPAIR)='surface pressure       '
      vname(6,indxPAIR)='lat_rho lon_rho                            '
      vname(7,indxPAIR)='                                        '
      vname(1,indxQBFC)='QBFC                                   '
      vname(2,indxQBFC)='Quadratic bottom friction coef.         '
      vname(3,indxQBFC)='                                        '
      vname(4,indxQBFC)='                                        '
      vname(5,indxQBFC)='quadratic bottom friction coef.         '
      vname(6,indxQBFC)='lat_rho lon_rho                            '
      vname(7,indxQBFC)='                                        '
      vname(1,indxSSH)='SSH                                       '
      vname(2,indxSSH)='sea surface height                        '
      vname(3,indxSSH)='meter                                     '
      vname(4,indxSSH)='SSH, scalar, series                       '
      vname(5,indxSSH)='sea_surface_height_above_sea_level        '
      vname(6,indxSSH)='lat_rho lon_rho                           '
      vname(7,indxSSH)='                                          '
      vname(1,indxSUSTR)='sustr                                   '
      vname(2,indxSUSTR)='surface u-momentum stress               '
      vname(3,indxSUSTR)='Newton meter-2                          '
      vname(4,indxSUSTR)='surface u-mom. stress, scalar, series   '
      vname(5,indxSUSTR)='surface_downward_x_stress               '
      vname(6,indxSUSTR)='lat_u lon_u                             '
      vname(7,indxSUSTR)='                                        '
      vname(1,indxSVSTR)='svstr                                   '
      vname(2,indxSVSTR)='surface v-momentum stress               '
      vname(3,indxSVSTR)='Newton meter-2                          '
      vname(4,indxSVSTR)='surface v-mom. stress, scalar, series   '
      vname(5,indxSVSTR)='surface_downward_y_stress               '
      vname(6,indxSVSTR)='lat_v lon_v                             '
      vname(7,indxSVSTR)='                                        '
        vname(1,indxQBAR)='Qbar'
        vname(2,indxQBAR)='Runoff Discharge'
        vname(3,indxQBAR)='m3.s-1'
        vname(4,indxQBAR)='Qbar,scalar, series'
      return
      end
