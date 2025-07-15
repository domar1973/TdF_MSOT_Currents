      subroutine setup_kwds (ierr)
      implicit none
      integer*4 ierr, is,ie
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
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      is=1
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='title'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='start_date'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='time_stepping'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='initial'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='grid'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='forcing'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='restart'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='history'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='primary_history_fields'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='rho0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='bottom_drag'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='gamma2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='nudg_cof'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='psource_ncfile'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      return
  99   write(stdout,'(/1x,A,A/14x,A)')
     &  'SETUP_KWDS ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
