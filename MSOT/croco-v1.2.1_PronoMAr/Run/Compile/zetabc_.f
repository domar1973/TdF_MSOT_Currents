      subroutine zetabc_tile(Istr,Iend,Jstr,Jend)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j
      real cff,cx,dft,dfx, eps
      parameter (eps=1.D-20)
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
      if (istr.eq.1) then
        do j=JstrV-1,Jend
          cx=dtfast*pm(1,j)*sqrt(g*(h(1,j)+zeta(1,j,kstp)))
          zeta(0,j,knew)=(zeta(0,j,kstp)+cx*zeta(1,j,knew))/(1.D0+cx)
     &                      *rmask(0,j)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=JstrV-1,Jend
          cx=dtfast*pm(Iend,j)*sqrt(g*(h(Iend,j)+zeta(Iend,j,kstp)))
          zeta(Iend+1,j,knew)=(zeta(Iend+1,j,kstp)
     &                              +cx*zeta(Iend,j,knew))/(1.D0+cx)
     &                      *rmask(Iend+1,j)
        enddo
      endif
      if (jstr.eq.1) then
        do i=IstrU-1,Iend
          cx=dtfast*pn(i,1)*sqrt(g*(h(i,1)+zeta(i,1,kstp)))
          zeta(i,0,knew)=(zeta(i,0,kstp)
     &                              +cx*zeta(i,1,knew))/(1.D0+cx)
     &                      *rmask(i,0)
        enddo
      endif
      if (jend.eq.Mm) then
        do i=IstrU-1,Iend
          cx=dtfast*pn(i,Jend)*sqrt(g*(h(i,Jend)+zeta(i,Jend,kstp)))
          zeta(i,Jend+1,knew)=(zeta(i,Jend+1,kstp)
     &                              +cx*zeta(i,Jend,knew))/(1.D0+cx)
     &                      *rmask(i,Jend+1)
        enddo
      endif
      if (jstr.eq.1 .and. istr.eq.1) then
        zeta(0,0,knew)=0.5D0*(zeta(1,0 ,knew)+zeta(0 ,1,knew))
     &                       *rmask(0,0)
      endif
      if (jstr.eq.1 .and. iend.eq.Lm) then
        zeta(Iend+1,0,knew)=0.5D0*(zeta(Iend+1,1 
     &                                         ,knew)+zeta(Iend,0,knew))
     &                       *rmask(Iend+1,0)
      endif
      if (jend.eq.Mm .and. istr.eq.1) then
        zeta(0,Jend+1,knew)=0.5D0*(zeta(0,Jend,knew)+zeta(1 
     &                                                    ,Jend+1,knew))
     &                       *rmask(0,Jend+1)
      endif
      if (jend.eq.Mm .and. iend.eq.Lm) then
        zeta(Iend+1,Jend+1,knew)=0.5D0*(zeta(Iend+1,Jend,knew)+
     &                            zeta(Iend,Jend+1,knew))
     &                       *rmask(Iend+1,Jend+1)
      endif
      if (istr.eq.1) then
        DO j=JstrV-1,Jend
          IF (zeta(Istr-1,j,knew).le.
     &        (Dcrit(Istr-1,j)-h(Istr-1,j))) THEN
            zeta(Istr-1,j,knew)=Dcrit(Istr-1,j)-h(Istr-1,j)-eps
          END IF
        END DO
      END IF
      if (iend.eq.Lm) then
        DO j=JstrV-1,Jend
          IF (zeta(Iend+1,j,knew).le.
     &        (Dcrit(Iend+1,j)-h(Iend+1,j))) THEN
            zeta(Iend+1,j,knew)=Dcrit(Iend+1,j)-h(Iend+1,j)-eps
          END IF
        END DO
      END IF
      if (jstr.eq.1) then
        DO i=IstrU-1,Iend
          IF (zeta(i,Jstr-1,knew).le.
     &       (Dcrit(i,Jstr-1)-h(i,Jstr-1))) THEN
               zeta(i,Jstr-1,knew)=Dcrit(i,Jstr-1)-h(i,Jstr-1)-eps
          END IF
        END DO
      END IF
      if (jend.eq.Mm) then
        DO i=IstrU-1,Iend
          IF (zeta(i,Jend+1,knew).le.
     &       (Dcrit(i,Jend+1)-h(i,Jend+1))) THEN
               zeta(i,Jend+1,knew)=Dcrit(i,Jend+1)-h(i,Jend+1)-eps
          END IF
        END DO
      END IF
      if (jstr.eq.1 .and. istr.eq.1) then
        IF (zeta(Istr-1,Jstr-1,knew).le.
     &      (Dcrit(Istr-1,Jstr-1)-h(Istr-1,Jstr-1))) THEN
            zeta(Istr-1,Jstr-1,knew)=
     &            Dcrit(Istr-1,Jstr-1)-h(Istr-1,Jstr-1)-eps
        END IF
      END IF
      if (jstr.eq.1 .and. iend.eq.Lm) then
        IF (zeta(Iend+1,Jstr-1,knew).le.
     &      (Dcrit(Iend+1,Jstr-1)-h(Iend+1,Jstr-1))) THEN
            zeta(Iend+1,Jstr-1,knew)=
     &            Dcrit(Iend+1,Jstr-1)-h(Iend+1,Jstr-1)-eps
        END IF
      END IF
      if (jend.eq.Mm .and. istr.eq.1) then
        IF (zeta(Istr-1,Jend+1,knew).le.
     &      (Dcrit(Istr-1,Jend+1)-h(Istr-1,Jend+1))) THEN
            zeta(Istr-1,Jend+1,knew)=
     &            Dcrit(Istr-1,Jend+1)-h(Istr-1,Jend+1)-eps
        END IF
      END IF
      if (jend.eq.Mm .and. iend.eq.Lm) then
        IF (zeta(Iend+1,Jend+1,knew).le.
     &      (Dcrit(Iend+1,Jend+1)-h(Iend+1,Jend+1))) THEN
            zeta(Iend+1,Jend+1,knew)=
     &            Dcrit(Iend+1,Jend+1)-h(Iend+1,Jend+1)-eps
        END IF
      END IF
      return
      end
