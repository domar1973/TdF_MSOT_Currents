      subroutine check_switches1 (ierr)
      implicit none
      integer*4 ierr, is,ie, iexample
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
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
       write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
      iexample=iexample+1
       write(stdout,'(10x,A)') 'Argentina'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Argentina'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'OPENMP'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OPENMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'START_DATE'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_DATE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TIDES2D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDES2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'OBC_EAST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_EAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'OBC_WEST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'OBC_NORTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NORTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'OBC_SOUTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_SOUTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'CURVGRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CURVGRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'SPHERICAL'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPHERICAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'MASKING'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MASKING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'WET_DRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WET_DRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_COR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_COR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ATM_PRESS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ATM_PRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'QUADRATIC_BOTTOM_FRICTION_COEF'
      ie=is +29
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QUADRATIC_BOTTOM_FRICTION_COEF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_PSOURCE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_PSOURCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'PSOURCE'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PSOURCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'PSOURCE_NCFILE'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PSOURCE_NCFILE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'ANA_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'Z_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'SSH_TIDES'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SSH_TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'SWACS'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SWACS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'FRC_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'MPI_COMM_WORLD'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI_COMM_WORLD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'M2FILTER_POWER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_POWER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'HZR'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HZR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'M2_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'UV_VADV_SPLINES'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'NTRA_T3DMIX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NTRA_T3DMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'TS_VADV_AKIMA'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_AKIMA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'BHFLUX'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BHFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'OBC_M2CHARACT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2CHARACT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(10x,A)') 'NF_CLOBBER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_CLOBBER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
       write(stdout,'(/)')
      if (iexample.eq.0) then
         write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
         write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99   write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
