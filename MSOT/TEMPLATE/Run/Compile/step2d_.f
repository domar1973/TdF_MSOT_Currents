      subroutine step2d (tile)
      implicit none
      integer*4 tile, trd
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
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,7,0:NPP-1)
      common/private_scratch/ A2d,A3d
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
      trd=0
C$    trd=omp_get_thread_num()
      call step2D_FB_tile ( Istr,Iend,Jstr,Jend, A2d(1,1,trd)
     &                    , A2d(1, 2,trd), A2d(1, 3,trd), A2d(1, 4,trd)
     &                    , A2d(1, 5,trd), A2d(1, 6,trd), A2d(1, 7,trd)
     &                    , A2d(1, 8,trd), A2d(1, 9,trd)
     &                    , A2d(1,10,trd), A2d(1,11,trd)
     &                    )
      return
      end
      subroutine step2D_FB_tile (Istr,Iend,Jstr,Jend, zeta_new,
     &                           Dnew,rubar,rvbar,
     &                           Drhs, UFx,UFe,
     &                           VFx,VFe
     &                          ,wrk1,wrk2
     &                          )
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
      integer*4 Istr,Iend,Jstr,Jend, i,j, kbak, kold,
     &        is,
     &        imin,imax,jmin,jmax
      real sum_c
      real    VMAX,VMAXL
      real zeta_new(Istr-2:Iend+2,Jstr-2:Jend+2),  cff,
     &         Dnew(Istr-2:Iend+2,Jstr-2:Jend+2),  cff0,
     &        rubar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff1,
     &        rvbar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff2,
     &         Drhs(Istr-2:Iend+2,Jstr-2:Jend+2),  cff3,
     &          UFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &          UFe(Istr-2:Iend+2,Jstr-2:Jend+2),  DUnew,
     &          VFx(Istr-2:Iend+2,Jstr-2:Jend+2),  DVnew,
     &          VFe(Istr-2:Iend+2,Jstr-2:Jend+2)
      real     wrk1(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         wrk2(Istr-2:Iend+2,Jstr-2:Jend+2)
      real     curvX,curvE,cffE,cffX, gamma
      parameter (gamma = -0.25D0)
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
      real cff1_WD,cff2_WD
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if (istr.eq.1) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (iend.eq.Lm) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (jstr.eq.1) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (jend.eq.Mm) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      if (iic.eq.ntstart) then
        kbak=kstp
        kold=kstp
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iic.eq.ntstart+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1= 1.781105D0
        cff2=-1.06221D0
        cff3= 0.281105D0
      endif
      imin=max(IstrU-3,0)
      imax=min(Iend+2,Lm+1)
      jmin=max(JstrV-3,0)
      jmax=min(Jend+2,Mm+1)
      do j=jmin,jmax
        do i=imin,imax
          Drhs(i,j)=cff1*zeta(i,j,kstp)+cff2*zeta(i,j,kbak)
     &                                 +cff3*zeta(i,j,kold)
     &                                             + h(i,j)
        enddo
      enddo
      do j=max(Jstr-2,0),min(Jend+2,Mm+1)
        do i=imin+1,imax
          urhs(i,j)=cff1*ubar(i,j,kstp) +cff2*ubar(i,j,kbak)
     &                                  +cff3*ubar(i,j,kold)
          DUon(i,j)=0.5D0*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
     &                                                              )
        enddo
      enddo
      do j=jmin+1,jmax
        do i=max(Istr-2,0),min(Iend+2,Lm+1)
          vrhs(i,j)=cff1*vbar(i,j,kstp) +cff2*vbar(i,j,kbak)
     &                                  +cff3*vbar(i,j,kold)
          DVom(i,j)=0.5D0*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
     &                                                              )
        enddo
      enddo
      do is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        if ( (imin+1) .le.i .and. i.le.imax .and.
     &       (jmin+1) .le.j .and. j.le.jmax) then
          if (Dsrc(is).eq.0) then
            urhs(i,j)=2.D0*Qbar(is)/( on_u(i,j)
     &                             *(Drhs(i-1,j)+Drhs(i,j)) )
            DUon(i,j)=Qbar(is)
          else
            vrhs(i,j)=2.D0*Qbar(is)/( om_v(i,j)
     &                             *(Drhs(i,j-1)+Drhs(i,j)) )
            DVom(i,j)=Qbar(is)
          endif
        endif
      enddo
      if (iic.eq.ntstart) then
        cff0=1.D0
        cff1=0.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iic.eq.ntstart+1) then
        cff0= 1.0833333333333D0
        cff1=-0.1666666666666D0
        cff2= 0.0833333333333D0
        cff3= 0.D0
      else
        cff0=0.614D0
        cff1=0.285D0
        cff2=0.088D0
        cff3=0.013D0
      endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1))
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
          cff=0.5D0+SIGN(0.5D0,Dcrit(i,j)-h(i,j))
          zeta_new(i,j)=zeta_new(i,j)+
     &                 cff*(Dcrit(i,j)-h(i,j))*(1.D0-rmask(i,j))
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          UFx(i,j)=cff0*zeta_new(i,j) +cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold)
          UFe(i,j)=UFx(i,j)
          VFe(i,j)=UFx(i,j)*UFx(i,j)
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=zeta_new(i,j)+h(i,j)
          zeta(i,j,knew)=zeta_new(i,j)
        enddo
      enddo
      call wetdry_tile (Istr,Iend,Jstr,Jend)
C$OMP BARRIER
      call zetabc_tile (Istr,Iend,Jstr,Jend)
      cff=0.5D0*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*(
     &                         (h(i-1,j)+h(i,j))*(UFe(i-1,j)
     &                        -UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &                                                              )
          rvbar(i,j)=cff*om_v(i,j)*(
     &            (h(i,j-1)+h(i,j))*(UFe(i,j-1)
     &                        -UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &                                                              )
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-(h(i-1,j)+h(i,j)+UFe(i-1,j)+UFe(i,j))*
     & (pair2d(i,j)-pair2d(i-1,j))*on_u(i,j)*0.5D0
          rvbar(i,j)=rvbar(i,j)-(h(i,j-1)+h(i,j)+UFe(i,j-1)+UFe(i,j))*
     & (pair2d(i,j)-pair2d(i,j-1))*om_v(i,j)*0.5D0
        enddo
      enddo
        do j=Jstr,Jend
          do i=max(IstrU-1,2),min(Iend+1,Lm)
            wrk1(i,j)=urhs(i-1,j)-2.D0*urhs(i,j)+urhs(i+1,j)
            wrk2(i,j)=Duon(i-1,j)-2.D0*Duon(i,j)+Duon(i+1,j)
          enddo
        enddo
        if (istr.eq.1) then
          do j=Jstr,Jend
            wrk1(1,j)=wrk1(2,j)
            wrk2(1,j)=wrk2(2,j)
          enddo
        endif
        if (iend.eq.Lm) then
          do j=Jstr,Jend
            wrk1(Lm+1,j)=wrk1(Lm,j)
            wrk2(Lm+1,j)=wrk2(Lm,j)
          enddo
        endif
        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cffX=urhs(i,j)+urhs(i+1,j)
            if (cffX.gt.0.D0) then
              curvX=wrk1(i,j)
            else
              curvX=wrk1(i+1,j)
            endif
            UFx(i,j)=0.25D0*(cffX+gamma*curvX)*( Duon(i,j)+
     &               Duon(i+1,j)-0.125D0*(wrk2(i,j)+wrk2(i+1,j)))
          enddo
        enddo
        do j=max(JstrV-1,2),min(Jend+1,Mm)
          do i=Istr,Iend
            wrk1(i,j) =vrhs(i,j-1)-2.D0*vrhs(i,j)+vrhs(i,j+1)
            wrk2(i,j)=Dvom(i,j-1)-2.D0*Dvom(i,j)+Dvom(i,j+1)
          enddo
        enddo
        if (jstr.eq.1) then
          do i=Istr,Iend
            wrk1(i,1)=wrk1(i,2)
            wrk2(i,1)=wrk2(i,2)
          enddo
        endif
        if (jend.eq.Mm) then
          do i=Istr,Iend
            wrk1(i,Mm+1)=wrk1(i,Mm)
            wrk2(i,Mm+1)=wrk2(i,Mm)
          enddo
        endif
        do j=JstrV-1,Jend
          do i=Istr,Iend
            cffE=vrhs(i,j)+vrhs(i,j+1)
            if (cffE.gt.0.D0) then
              curvE=wrk1(i,j)
            else
              curvE=wrk1(i,j+1)
            endif
            VFe(i,j)=0.25D0*(cffE+gamma*curvE)*( Dvom(i,j)+
     &               Dvom(i,j+1)-0.125D0*(wrk2(i,j)+wrk2(i,j+1)))
          enddo
        enddo
        do j=max(Jstr-1,1),min(Jend+1,Mm)
          do i=IstrU,Iend
            wrk1(i,j)=urhs(i,j-1)-2.D0*urhs(i,j)+urhs(i,j+1)
          enddo
        enddo
        if (jstr.eq.1) then
          do i=IstrU,Iend
            wrk1(i,0)=wrk1(i,1)
          enddo
        endif
        if (jend.eq.Mm) then
          do i=IstrU,Iend
            wrk1(i,Mm+1)=wrk1(i,Mm)
          enddo
        endif
        do j=Jstr,Jend+1
          do i=IstrU-1,Iend
           wrk2(i,j)=Dvom(i-1,j)-2.D0*Dvom(i,j)+Dvom(i+1,j)
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=IstrU,Iend
            cffX=urhs(i,j)+urhs(i,j-1)
            cffE=Dvom(i,j)+Dvom(i-1,j)
            if (cffE.gt.0.D0) then
              curvX=wrk1(i,j-1)
            else
              curvX=wrk1(i,j)
            endif
            UFe(i,j)=0.25D0*(cffX+gamma*curvX)*(cffE-0.125D0*(
     &                             wrk2(i,j)+wrk2(i-1,j) ))
          enddo
        enddo
        do j=JstrV,Jend
          do i=max(Istr-1,1),min(Iend+1,Lm)
            wrk1(i,j)=vrhs(i-1,j)-2.D0*vrhs(i,j)+vrhs(i+1,j)
          enddo
        enddo
        if (istr.eq.1) then
          do j=JstrV,Jend
            wrk1(0,j)=wrk1(1,j)
          enddo
        endif
        if (iend.eq.Lm) then
          do j=JstrV,Jend
            wrk1(Lm+1,j)=wrk1(Lm,j)
          enddo
        endif
        do j=JstrV-1,Jend
          do i=Istr,Iend+1
           wrk2(i,j)=Duon(i,j-1)-2.D0*Duon(i,j)+Duon(i,j+1)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend+1
            cffE=vrhs(i,j)+vrhs(i-1,j)
            cffX=Duon(i,j)+Duon(i,j-1)
            if (cffX.gt.0.D0) then
              curvE=wrk1(i-1,j)
            else
              curvE=wrk1(i,j)
            endif
            VFx(i,j)=0.25D0*(cffE+gamma*curvE)*(cffX-0.125D0*(
     &                             wrk2(i,j)+wrk2(i,j-1) ))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j)-UFx(i,j  )+UFx(i-1,j)
     &                           -UFe(i,j+1)+UFe(i  ,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j  )
     &                           -VFe(i  ,j)+VFe(i,j-1)
          enddo
        enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=Drhs(i,j)*(
     &                   fomn(i,j)
     &          +0.5D0*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &                -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
     &                   )
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar(i,j)=rubar(i,j)+0.25D0*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j)-0.25D0*(VFe(i,j)+VFe(i,j-1))
        enddo
      enddo
      if (rdrg2.gt.0.D0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25D0*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
            rdrg2=0.5D0*(QBFC(i,j) + QBFC(i-1,j))
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25D0*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
            rdrg2=0.5D0*(QBFC(i,j) + QBFC(i,j-1))
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0D0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          DUon(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
      cff=0.5D0*dtfast
      cff2=2.D0*dtfast
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew=( (DUon(i,j)+DUon(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                          *rubar(i,j)+cff2*sustr(i,j)
     &                                                    )
     &                                         *umask(i,j)
          cff1_WD=ABS(ABS(umask_wet(i,j))-1.D0)
          cff2_WD=0.5D0+SIGN(0.5D0,DUnew)*umask_wet(i,j)
          umask_wet(i,j)=0.5D0*umask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.D0-cff1_WD)
          DUnew=DUnew*umask_wet(i,j)
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          DVnew=( (DUon(i,j)+DUon(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                          *rvbar(i,j)+cff2*svstr(i,j)
     &                                                    )
     &                                         *vmask(i,j)
          cff1_WD=ABS(ABS(vmask_wet(i,j))-1.D0)
          cff2_WD=0.5D0+SIGN(0.5D0,DVnew)*vmask_wet(i,j)
          vmask_wet(i,j)=0.5D0*vmask_wet(i,j)*cff1_WD
     &                        +cff2_WD*(1.D0-cff1_WD)
          DVnew=DVnew*vmask_wet(i,j)
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
        enddo
      enddo
      call u2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
      call v2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
      do is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            ubar(i,j,knew)=2.D0*Qbar(is)/( on_u(i,j)
     &                       *(Dnew(i-1,j)+Dnew(i,j)) )
          else
            vbar(i,j,knew)=2.D0*Qbar(is)/( om_v(i,j)
     &                       *(Dnew(i,j-1)+Dnew(i,j)) )
          endif
        endif
      enddo
      call diag_tile (Istr,Iend,Jstr,Jend, UFx,UFe)
      VMAXL=100.D0
      VMAX=0.D0
      do j=Jstr,Jend
        do i=Istr,Iend
          cff1=ubar(i,j,knew)
          cff2=vbar(i,j,knew)
          cff=max(abs(cff1),abs(cff2))
          IF (cff.GE.VMAX .or. cff1.ne.cff1 .or. cff2.ne.cff2) THEN
            IF (cff.GE.VMAX .and. cff1.eq.cff1 .and. cff2.eq.cff2) THEN
              VMAX=cff
            ELSE
              VMAX=666.D0
            ENDIF
            imax=i
            jmax=j
          ENDIF
        enddo
      enddo
      IF (VMAX.GT.VMAXL) THEN
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =   STEP2D:   ABNORMAL JOB END        = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '
        if (VMAX.eq.666.D0) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')
     &                                       '  IMAX JMAX  =',imax,jmax
        write(stdout,'(A,I6/)')    '  IIC        =',iic
        may_day_flag=1
        stop
      ENDIF
      return
      end
