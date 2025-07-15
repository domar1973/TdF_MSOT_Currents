      subroutine read_inp (ierr)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      integer*4 kwsize, testunit, input
      parameter (kwsize=32, testunit=40, input=15)
      character end_signal*3, keyword*32, fname*180
      parameter (end_signal='end')
      integer*4 ierr, iargc, is,ie, kwlen, lstr, lenstr
      logical dumboolean
      fname='croco.in'
      if (iargc().eq.1) call getarg(1,fname)
      wrthis(indxTime)=.false.
      ierr=0
      call setup_kwds (ierr)
      open (input,file=fname,status='old',form='formatted',err=97)
   1  keyword='                                '
      read(input,'(A)',err=1,end=99) keyword
      if (ichar(keyword(1:1)).eq.33) goto 1
      is=1
   2  if (is.eq.kwsize) then
        goto 1
      elseif (keyword(is:is).eq.' ') then
        is=is+1
        goto 2
      endif
      ie=is
   3  if (keyword(ie:ie).eq.':') then
        keyword(ie:ie)=' '
        goto 4
      elseif (keyword(ie:ie).ne.' ' .and. ie.lt.kwsize) then
        ie=ie+1
        goto 3
      endif
      goto 1
   4  kwlen=ie-is
      if (is.gt.1) keyword(1:kwlen)=keyword(is:is+kwlen-1)
      if (keyword(1:kwlen).eq.'title') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) title
        lstr=lenstr(title)
         write(stdout,'(/1x,A)') title(1:lstr)
      elseif (keyword(1:kwlen).eq.'start_date') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) start_date
        lstr=lenstr(start_date)
         write(stdout,'(/2x,A,1x,A/)')
     &   'Start date:', start_date(1:lstr)
      elseif (keyword(1:kwlen).eq.'time_stepping') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntimes,dt,ndtfast, ninfo
         write(stdout,
     & '(I10,2x,A,1x,A /F10.2,2x,A,2(/I10,2x,A,1x,A)/F10.4,2x,A)')
     &    ntimes,  'ntimes   Total number of timesteps for',
     &                                            '3D equations.',
     &    dt,      'dt       Timestep [sec] for 3D equations',
     &    ndtfast, 'ndtfast  Number of 2D timesteps within each',
     &                                                 '3D step.',
     &    ninfo,   'ninfo    Number of timesteps between',
     &                                     'runtime diagnostics.'
        dtfast=dt/float(ndtfast)
        if (NWEIGHT.lt.(2*ndtfast-1)) then
          write(stdout,'(a,i3)')
     &    ' Error - Number of 2D timesteps (2*ndtfast-1): ',
     &    2*ndtfast-1
          write(stdout,'(a,i3)')
     &    '           exceeds barotopic weight dimension: ',NWEIGHT
          goto 95
        endif
      elseif (keyword(1:kwlen).eq.'initial') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrrec
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          ininame=fname(1:lstr)
           write(stdout,'(1x,A,2x,A,4x,A,I3)')
     &     'Initial State File:', ininame(1:lstr), 'Record:',nrrec
      elseif (keyword(1:kwlen).eq.'grid') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open(testunit,file=fname(1:lstr), status='old', err=97)
        close(testunit)
        grdname=fname(1:lstr)
         write(stdout,'(10x,A,2x,A)')
     &                   'Grid File:', grdname(1:lstr)
      elseif (keyword(1:kwlen).eq.'forcing') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        frcname=fname(1:lstr)
         write(stdout,'(2x,A,2x,A)')
     &             'Forcing Data File:', frcname(1:lstr)
      elseif (keyword(1:kwlen).eq.'restart') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrst, nrpfrst
        read(input,'(A)',err=95)  fname
        lstr=lenstr(fname)
        rstname=fname(1:lstr)
         write(stdout,
     &             '(7x,A,2x,A,4x,A,I6,4x,A,I4)')
     &             'Restart File:', rstname(1:lstr),
     &             'nrst =', nrst, 'rec/file: ', nrpfrst
      elseif (keyword(1:kwlen).eq.'history') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefhis, nwrt, nrpfhis
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        hisname=fname(1:lstr)
         write(stdout,
     &             '(7x,A,2x,A,2x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &       'History File:', hisname(1:lstr),  'Create new:',
     &       ldefhis, 'nwrt =', nwrt, 'rec/file =', nrpfhis
      elseif (keyword(1:kwlen).eq.'primary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxZ),  wrthis(indxUb)
     &                                       ,  wrthis(indxVb)
     &                     , wrthis(indxQBFC)
        if ( wrthis(indxZ) .or. wrthis(indxUb) .or. wrthis(indxVb)
     &                        .or. wrthis(indxQBFC)
     &     ) wrthis(indxTime)=.true.
         write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &    , wrthis(indxZ),  'write zeta ', 'free-surface.'
     &    , wrthis(indxUb), 'write UBAR ', '2D U-momentum component.'
     &    , wrthis(indxVb), 'write VBAR ', '2D V-momentum component.'
     &    ,wrthis(indxQBFC),'write QBFC ', 'Quad. bottom fric. coef.'
      elseif (keyword(1:kwlen).eq.'rho0') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rho0
         write(stdout,'(F10.4,2x,A,1x,A)')
     &        rho0, 'rho0     Boussinesq approximation',
     &                           'mean density, kg/m3.'
      elseif (keyword(1:kwlen).eq.'bottom_drag') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rdrg, rdrg2, Zobt, Cdb_min, Cdb_max
         write(stdout,'(5(1pe10.3,2x,A/))')
     &     rdrg, 'rdrg     Linear bottom drag coefficient (m/si).',
     &    rdrg2, 'rdrg2    Quadratic bottom drag coefficient.',
     &     Zobt, 'Zobt     Bottom roughness for logarithmic law (m).',
     &  Cdb_min, 'Cdb_min  Minimum bottom drag coefficient.',
     &  Cdb_max, 'Cdb_max  Maximum bottom drag coefficient.'
      elseif (keyword(1:kwlen).eq.'gamma2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) gamma2
         write(stdout,'(f10.2,2x,A,1x,A)')
     &     gamma2, 'gamma2   Slipperiness parameter:',
     &                     'free-slip +1, or no-slip -1.'
      elseif (keyword(1:kwlen).eq.'nudg_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) tauT_in,tauT_out,tauM_in,tauM_out
          tauT_in =1.D0/(tauT_in *86400.D0)
          tauT_out=1.D0/(tauT_out*86400.D0)
          tauM_in =1.D0/(tauM_in *86400.D0)
          tauM_out=1.D0/(tauM_out*86400.D0)
           write(stdout,'(1pe10.3,2x,A)')
     &        tauT_in,'tauT_in  Nudging coefficients [sec^-1]'
           write(stdout,'(1pe10.3,2x,A)')
     &       tauT_out,'tauT_out Nudging coefficients [sec^-1]'
           write(stdout,'(1pe10.3,2x,A)')
     &        tauM_in,'tauM_in  Nudging coefficients [sec^-1]'
           write(stdout,'(1pe10.3,2x,A/)')
     &       tauM_out,'tauM_out Nudging coefficients [sec^-1]'
      elseif (keyword(1:kwlen).eq.'psource_ncfile') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        qbarname=fname(1:lstr)
         write(stdout,'(2x,A,2x,A)')
     &             'Runoff Data File:', qbarname(1:lstr)
        read(input,*,err=95) Nsrc
         write(stdout,'(/6x,i6,2x,A,1x,A)')
     &                               Nsrc, 'Number of point sources'
        do is=1,Nsrc
          read(input,*,err=95) Isrc(is),Jsrc(is),Dsrc(is),qbardir(is),
     &                                     (Lsrc(is,itrc), itrc=1,NT),
     &                                     (Tsrc0(is,itrc), itrc=1,NT)
         write(stdout,'(3(/6x,i6,2x,A),(/6x,F6.0,2x,A))')
     &    Isrc(is), 'I point source indice'
     &  , Jsrc(is), 'J point source indice'
     &  , Dsrc(is), 'Orientation of point source flow'
     &  , qbardir(is), 'Direction of point source flow'
      enddo
      else
         write(stdout,'(/3(1x,A)/)')
     &                  'WARNING: Unrecognized keyword:',
     &                   keyword(1:kwlen),' --> DISREGARDED.'
      endif
      if (keyword(1:kwlen) .eq. end_signal) goto 99
      goto 1
  95  write(stdout,'(/1x,4A/)') 'READ_INP ERROR while reading block',
     &                    ' with keyword ''', keyword(1:kwlen), '''.'
      ierr=ierr+1
      goto 99
  97  write(stdout,'(/1x,4A/)') 'READ_INP ERROR: Cannot find input ',
     &                                'file ''', fname(1:lstr), '''.'
      ierr=ierr+1
  99  close (input)
      if (ierr.eq.0) then
        call check_kwds (ierr)
        call check_srcs
        call check_switches1 (ierr)
        call check_switches2 (ierr)
      endif
      if (ierr.ne.0) then
        write(stdout,'(/1x,2A,I3,1x,A/)') 'READ_INP ERROR: ',
     & 'A total of', ierr, 'configuration errors discovered.'
        return
      endif
      return
      end
