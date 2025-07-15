#include "cppdefs.h"


      subroutine check_switches1 (ierr)

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if new CPP-switches  !!!!!!
!!!!!! were introduced into "cppdefs.h".  NO ACTION IS NEEDED  !!!!!!
!!!!!! if changes in "cppdefs.h" are limited to activation or  !!!!!!
!!!!!! deactivation of previously known switches.              !!!!!!
!!!!!! To refresh this file compile and execute "cppcheck.F"   !!!!!!
!!!!!! as an independent program, or use commands              !!!!!!
!!!!!! "make checkdefs" or "make depend".                      !!!!!!
!!!!!! Number of Configuration Choices:                      2 !!!!!!
!!!!!! Total number of CPP-switches:                       313 !!!!!!

      implicit none
      integer ierr, is,ie, iexample
#include "param.h"
#include "strings.h"
#ifdef MPI
# include "scalars.h"
#endif
      MPI_master_only write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
#ifdef Argentina
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'Argentina'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Argentina'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Argentina_Anidado1
      iexample=iexample+1
      MPI_master_only write(stdout,'(10x,A)') 'Argentina_Anidado1'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Argentina_Anidado1'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OPENMP
      MPI_master_only write(stdout,'(10x,A)') 'OPENMP'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OPENMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MPI
      MPI_master_only write(stdout,'(10x,A)') 'MPI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef START_DATE
      MPI_master_only write(stdout,'(10x,A)') 'START_DATE'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='START_DATE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef EXACT_RESTART
      MPI_master_only write(stdout,'(10x,A)') 'EXACT_RESTART'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EXACT_RESTART'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TIDES2D
      MPI_master_only write(stdout,'(10x,A)') 'TIDES2D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDES2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_EAST
      MPI_master_only write(stdout,'(10x,A)') 'OBC_EAST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_EAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_WEST
      MPI_master_only write(stdout,'(10x,A)') 'OBC_WEST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_NORTH
      MPI_master_only write(stdout,'(10x,A)') 'OBC_NORTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NORTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_SOUTH
      MPI_master_only write(stdout,'(10x,A)') 'OBC_SOUTH'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_SOUTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CURVGRID
      MPI_master_only write(stdout,'(10x,A)') 'CURVGRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CURVGRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPHERICAL
      MPI_master_only write(stdout,'(10x,A)') 'SPHERICAL'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPHERICAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MASKING
      MPI_master_only write(stdout,'(10x,A)') 'MASKING'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MASKING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WET_DRY
      MPI_master_only write(stdout,'(10x,A)') 'WET_DRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WET_DRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SOLVE3D
      MPI_master_only write(stdout,'(10x,A)') 'SOLVE3D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLVE3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_COR
      MPI_master_only write(stdout,'(10x,A)') 'UV_COR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_COR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_ADV
      MPI_master_only write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ATM_PRESS
      MPI_master_only write(stdout,'(10x,A)') 'ATM_PRESS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ATM_PRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'CLIMATOLOGY'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef QUADRATIC_BOTTOM_FRICTION_COEF
      MPI_master_only write(stdout,'(10x,A)') 'QUADRATIC_BOTTOM_FRICTION_COEF'
      ie=is +29
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QUADRATIC_BOTTOM_FRICTION_COEF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_PSOURCE
      MPI_master_only write(stdout,'(10x,A)') 'ANA_PSOURCE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_PSOURCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PSOURCE
      MPI_master_only write(stdout,'(10x,A)') 'PSOURCE'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PSOURCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PSOURCE_NCFILE
      MPI_master_only write(stdout,'(10x,A)') 'PSOURCE_NCFILE'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PSOURCE_NCFILE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BRY
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'Z_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SSH_TIDES
      MPI_master_only write(stdout,'(10x,A)') 'SSH_TIDES'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SSH_TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TIDERAMP
      MPI_master_only write(stdout,'(10x,A)') 'TIDERAMP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDERAMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SWACS
      MPI_master_only write(stdout,'(10x,A)') 'SWACS'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SWACS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'FRC_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'M2_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ALLOW_SINGLE_BLOCK_MODE
      MPI_master_only write(stdout,'(10x,A)') 'ALLOW_SINGLE_BLOCK_MODE'
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ALLOW_SINGLE_BLOCK_MODE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SINGLE
      MPI_master_only write(stdout,'(10x,A)') 'SINGLE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SINGLE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RVTK_DEBUG
      MPI_master_only write(stdout,'(10x,A)') 'RVTK_DEBUG'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RVTK_DEBUG'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RVTK_DEBUG_WRITE
      MPI_master_only write(stdout,'(10x,A)') 'RVTK_DEBUG_WRITE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RVTK_DEBUG_WRITE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RVTK_DEBUG_PERFRST
      MPI_master_only write(stdout,'(10x,A)') 'RVTK_DEBUG_PERFRST'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RVTK_DEBUG_PERFRST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RVTK_DEBUG_READ
      MPI_master_only write(stdout,'(10x,A)') 'RVTK_DEBUG_READ'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RVTK_DEBUG_READ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CONST_TRACERS
      MPI_master_only write(stdout,'(10x,A)') 'CONST_TRACERS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CONST_TRACERS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OA_COUPLING
      MPI_master_only write(stdout,'(10x,A)') 'OA_COUPLING'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OA_COUPLING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MPI_COMM_WORLD
      MPI_master_only write(stdout,'(10x,A)') 'MPI_COMM_WORLD'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI_COMM_WORLD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef READ_PATM
      MPI_master_only write(stdout,'(10x,A)') 'READ_PATM'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='READ_PATM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_PATM
      MPI_master_only write(stdout,'(10x,A)') 'OBC_PATM'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_PATM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OA_GRID_UV
      MPI_master_only write(stdout,'(10x,A)') 'OA_GRID_UV'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OA_GRID_UV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_FLUX
      MPI_master_only write(stdout,'(10x,A)') 'BULK_FLUX'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_FLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef QCORRECTION
      MPI_master_only write(stdout,'(10x,A)') 'QCORRECTION'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='QCORRECTION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SFLX_CORR
      MPI_master_only write(stdout,'(10x,A)') 'SFLX_CORR'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SFLX_CORR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_DIURNAL_SW
      MPI_master_only write(stdout,'(10x,A)') 'ANA_DIURNAL_SW'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_DIURNAL_SW'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SFLUX_CFB
      MPI_master_only write(stdout,'(10x,A)') 'SFLUX_CFB'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SFLUX_CFB'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OW_COUPLING
      MPI_master_only write(stdout,'(10x,A)') 'OW_COUPLING'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OW_COUPLING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WKB_WWAVE
      MPI_master_only write(stdout,'(10x,A)') 'WKB_WWAVE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WKB_WWAVE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_OFFLINE
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_OFFLINE'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_OFFLINE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MRL_WCI
      MPI_master_only write(stdout,'(10x,A)') 'MRL_WCI'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MRL_WCI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_ROLLER
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_ROLLER'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_ROLLER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_STREAMING
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_STREAMING'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_STREAMING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_RAMP
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_RAMP'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_RAMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_SMFLUX
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_SMFLUX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_SMFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef XIOS
      MPI_master_only write(stdout,'(10x,A)') 'XIOS'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='XIOS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_NONE
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_NONE'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_NONE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_POWER
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_POWER'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_POWER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_COSINE
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_COSINE'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_COSINE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2FILTER_FLAT
      MPI_master_only write(stdout,'(10x,A)') 'M2FILTER_FLAT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_FLAT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NO_TRACER
      MPI_master_only write(stdout,'(10x,A)') 'NO_TRACER'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NO_TRACER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TRACERS
      MPI_master_only write(stdout,'(10x,A)') 'TRACERS'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TRACERS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NO_TEMPERATURE
      MPI_master_only write(stdout,'(10x,A)') 'NO_TEMPERATURE'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NO_TEMPERATURE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TEMPERATURE
      MPI_master_only write(stdout,'(10x,A)') 'TEMPERATURE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TEMPERATURE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SALINITY
      MPI_master_only write(stdout,'(10x,A)') 'SALINITY'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SALINITY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ
      MPI_master_only write(stdout,'(10x,A)') 'NBQ'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3FAST
      MPI_master_only write(stdout,'(10x,A)') 'M3FAST'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3FAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_IMP
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_IMP'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_IMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_THETAIMP
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_THETAIMP'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_THETAIMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_FREESLIP
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_FREESLIP'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_FREESLIP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_HZ_PROGNOSTIC
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_HZ_PROGNOSTIC'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_HZ_PROGNOSTIC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3FAST_REINIT
      MPI_master_only write(stdout,'(10x,A)') 'M3FAST_REINIT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3FAST_REINIT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TANK
      MPI_master_only write(stdout,'(10x,A)') 'TANK'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TANK'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NOT_NBQ_AM4
      MPI_master_only write(stdout,'(10x,A)') 'NOT_NBQ_AM4'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NOT_NBQ_AM4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TRACETXT
      MPI_master_only write(stdout,'(10x,A)') 'TRACETXT'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TRACETXT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAG_CFL
      MPI_master_only write(stdout,'(10x,A)') 'DIAG_CFL'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAG_CFL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef HZR
      MPI_master_only write(stdout,'(10x,A)') 'HZR'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HZR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_PERF
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_PERF'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_PERF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_MASS
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_MASS'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_MASS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_GRID_SLOW
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_GRID_SLOW'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_GRID_SLOW'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_HZCORRECT
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_HZCORRECT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_HZCORRECT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_HADV_SPLINES
      MPI_master_only write(stdout,'(10x,A)') 'W_HADV_SPLINES'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_HADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_HADV_TVD
      MPI_master_only write(stdout,'(10x,A)') 'W_HADV_TVD'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_HADV_TVD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_HADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'W_HADV_WENO5'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_HADV_C4
      MPI_master_only write(stdout,'(10x,A)') 'W_HADV_C4'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_HADV_C4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_HADV_C2
      MPI_master_only write(stdout,'(10x,A)') 'W_HADV_C2'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_HADV_C2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_VADV_SPLINES
      MPI_master_only write(stdout,'(10x,A)') 'W_VADV_SPLINES'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_VADV_TVD
      MPI_master_only write(stdout,'(10x,A)') 'W_VADV_TVD'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_VADV_TVD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_VADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'W_VADV_WENO5'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_VADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_VADV_C2
      MPI_master_only write(stdout,'(10x,A)') 'W_VADV_C2'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_VADV_C2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_NBQ
      MPI_master_only write(stdout,'(10x,A)') 'OBC_NBQ'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NBQ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_NBQORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_NBQORLANSKI'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NBQORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_NBQSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_NBQSPECIFIED'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NBQSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_WORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_WORLANSKI'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_WSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_WSPECIFIED'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_NUDGING
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_NUDGING'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_NUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQCLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'NBQCLIMATOLOGY'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQCLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NBQ_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'NBQ_FRC_BRY'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef W_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'W_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BSTRESS_FAST
      MPI_master_only write(stdout,'(10x,A)') 'BSTRESS_FAST'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BSTRESS_FAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VAR_RHO_2D
      MPI_master_only write(stdout,'(10x,A)') 'VAR_RHO_2D'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_RHO_2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NONLIN_EOS
      MPI_master_only write(stdout,'(10x,A)') 'NONLIN_EOS'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NONLIN_EOS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef INNERSHELF
      MPI_master_only write(stdout,'(10x,A)') 'INNERSHELF'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='INNERSHELF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RESET_RHO0
      MPI_master_only write(stdout,'(10x,A)') 'RESET_RHO0'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RESET_RHO0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BASIN
      MPI_master_only write(stdout,'(10x,A)') 'BASIN'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BASIN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef EQUATOR
      MPI_master_only write(stdout,'(10x,A)') 'EQUATOR'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='EQUATOR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GRAV_ADJ
      MPI_master_only write(stdout,'(10x,A)') 'GRAV_ADJ'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GRAV_ADJ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PGF_FLAT_BOTTOM
      MPI_master_only write(stdout,'(10x,A)') 'PGF_FLAT_BOTTOM'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PGF_FLAT_BOTTOM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RIP
      MPI_master_only write(stdout,'(10x,A)') 'RIP'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RIP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PGF_BASIC_JACOBIAN
      MPI_master_only write(stdout,'(10x,A)') 'PGF_BASIC_JACOBIAN'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PGF_BASIC_JACOBIAN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WJ_GRADP
      MPI_master_only write(stdout,'(10x,A)') 'WJ_GRADP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WJ_GRADP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPLIT_EOS
      MPI_master_only write(stdout,'(10x,A)') 'SPLIT_EOS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPLIT_EOS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_UP3
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_C4
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_C4'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_C4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_C2
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_C2'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_C2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_UP5
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_UP5'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_UP5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_C6
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_C6'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_C6'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_HADV_TVD
      MPI_master_only write(stdout,'(10x,A)') 'UV_HADV_TVD'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_TVD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_MIX_S
      MPI_master_only write(stdout,'(10x,A)') 'UV_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_MIX_GEO
      MPI_master_only write(stdout,'(10x,A)') 'UV_MIX_GEO'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_MIX_GEO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VIS_SMAGO_3D
      MPI_master_only write(stdout,'(10x,A)') 'UV_VIS_SMAGO_3D'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS_SMAGO_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VIS2
      MPI_master_only write(stdout,'(10x,A)') 'UV_VIS2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_DIF2
      MPI_master_only write(stdout,'(10x,A)') 'TS_DIF2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VIS_SMAGO
      MPI_master_only write(stdout,'(10x,A)') 'UV_VIS_SMAGO'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS_SMAGO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_DIF_SMAGO
      MPI_master_only write(stdout,'(10x,A)') 'TS_DIF_SMAGO'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF_SMAGO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VIS_COEF_3D
      MPI_master_only write(stdout,'(10x,A)') 'VIS_COEF_3D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VIS_COEF_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLS_MIXING_3D
      MPI_master_only write(stdout,'(10x,A)') 'GLS_MIXING_3D'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_MIXING_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLS_MIXING
      MPI_master_only write(stdout,'(10x,A)') 'GLS_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIF_COEF_3D
      MPI_master_only write(stdout,'(10x,A)') 'DIF_COEF_3D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIF_COEF_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SOLITON
      MPI_master_only write(stdout,'(10x,A)') 'SOLITON'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLITON'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2_HADV_UP3
      MPI_master_only write(stdout,'(10x,A)') 'M2_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3CLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'M3CLIMATOLOGY'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3CLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMAT_UV_MIXH
      MPI_master_only write(stdout,'(10x,A)') 'CLIMAT_UV_MIXH'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMAT_UV_MIXH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VADV_SPLINES
      MPI_master_only write(stdout,'(10x,A)') 'UV_VADV_SPLINES'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'UV_VADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VADV_C2
      MPI_master_only write(stdout,'(10x,A)') 'UV_VADV_C2'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_C2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_VADV_TVD
      MPI_master_only write(stdout,'(10x,A)') 'UV_VADV_TVD'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_TVD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_ADAPT_IMP
      MPI_master_only write(stdout,'(10x,A)') 'VADV_ADAPT_IMP'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_ADAPT_IMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VADV_ADAPT_PRED
      MPI_master_only write(stdout,'(10x,A)') 'VADV_ADAPT_PRED'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VADV_ADAPT_PRED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_UP3
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_UP3'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_UP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_C4
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_C4'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_C4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_UP5
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_UP5'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_UP5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_C6
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_C6'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_C6'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_RSUP3
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_RSUP3'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_RSUP3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_HADV_RSUP5
      MPI_master_only write(stdout,'(10x,A)') 'TS_HADV_RSUP5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_RSUP5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_DIF4
      MPI_master_only write(stdout,'(10x,A)') 'TS_DIF4'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_DIF4'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_MIX_GEO
      MPI_master_only write(stdout,'(10x,A)') 'TS_MIX_GEO'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_GEO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_MIX_ISO
      MPI_master_only write(stdout,'(10x,A)') 'TS_MIX_ISO'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_ISO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_MIX_S
      MPI_master_only write(stdout,'(10x,A)') 'TS_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_MIX_IMP
      MPI_master_only write(stdout,'(10x,A)') 'TS_MIX_IMP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_IMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_MIX_ISO_FILT
      MPI_master_only write(stdout,'(10x,A)') 'TS_MIX_ISO_FILT'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_ISO_FILT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TCLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'TCLIMATOLOGY'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TCLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMAT_TS_MIXH
      MPI_master_only write(stdout,'(10x,A)') 'CLIMAT_TS_MIXH'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMAT_TS_MIXH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CLIMAT_TS_MIXH_FINE
      MPI_master_only write(stdout,'(10x,A)') 'CLIMAT_TS_MIXH_FINE'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLIMAT_TS_MIXH_FINE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIO_HADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'BIO_HADV_WENO5'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIO_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NTRA_T3DMIX
      MPI_master_only write(stdout,'(10x,A)') 'NTRA_T3DMIX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NTRA_T3DMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_VADV_SPLINES
      MPI_master_only write(stdout,'(10x,A)') 'TS_VADV_SPLINES'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_SPLINES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_VADV_AKIMA
      MPI_master_only write(stdout,'(10x,A)') 'TS_VADV_AKIMA'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_AKIMA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_VADV_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'TS_VADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_VADV_C2
      MPI_master_only write(stdout,'(10x,A)') 'TS_VADV_C2'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_C2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TS_VADV_FCT
      MPI_master_only write(stdout,'(10x,A)') 'TS_VADV_FCT'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_FCT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE_GRID
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE_GRID'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_GRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE_DIF2
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE_DIF2'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_DIF2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE_VIS2
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE_VIS2'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SEDIMENT
      MPI_master_only write(stdout,'(10x,A)') 'SEDIMENT'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SEDIMENT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SPONGE_SED
      MPI_master_only write(stdout,'(10x,A)') 'SPONGE_SED'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_SED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLS_KOMEGA
      MPI_master_only write(stdout,'(10x,A)') 'GLS_KOMEGA'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_KOMEGA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLS_KEPSILON
      MPI_master_only write(stdout,'(10x,A)') 'GLS_KEPSILON'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_KEPSILON'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GLS_GEN
      MPI_master_only write(stdout,'(10x,A)') 'GLS_GEN'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_GEN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CANUTO_A
      MPI_master_only write(stdout,'(10x,A)') 'CANUTO_A'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CANUTO_A'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef GibLau_78
      MPI_master_only write(stdout,'(10x,A)') 'GibLau_78'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GibLau_78'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MelYam_82
      MPI_master_only write(stdout,'(10x,A)') 'MelYam_82'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MelYam_82'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef KanCla_94
      MPI_master_only write(stdout,'(10x,A)') 'KanCla_94'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='KanCla_94'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Luyten_96
      MPI_master_only write(stdout,'(10x,A)') 'Luyten_96'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Luyten_96'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CANUTO_B
      MPI_master_only write(stdout,'(10x,A)') 'CANUTO_B'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CANUTO_B'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Cheng_02
      MPI_master_only write(stdout,'(10x,A)') 'Cheng_02'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Cheng_02'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TIDES
      MPI_master_only write(stdout,'(10x,A)') 'TIDES'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ZCLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'ZCLIMATOLOGY'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ZCLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SSH
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SSH'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SSH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef UV_TIDES
      MPI_master_only write(stdout,'(10x,A)') 'UV_TIDES'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_TIDES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M2CLIMATOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'M2CLIMATOLOGY'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2CLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_M2CLIMA
      MPI_master_only write(stdout,'(10x,A)') 'ANA_M2CLIMA'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_M2CLIMA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TIDES_MAS
      MPI_master_only write(stdout,'(10x,A)') 'TIDES_MAS'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TIDES_MAS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_MAKER
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_MAKER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_MAKER_JONSWAP
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_MAKER_JONSWAP'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_JONSWAP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_MAKER_GAUSSIAN
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_MAKER_GAUSSIAN'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_GAUSSIAN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_MAKER_SPECTRUM
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_MAKER_SPECTRUM'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_SPECTRUM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef STOKES_WAVES
      MPI_master_only write(stdout,'(10x,A)') 'STOKES_WAVES'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='STOKES_WAVES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef RIVER_RAMP
      MPI_master_only write(stdout,'(10x,A)') 'RIVER_RAMP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RIVER_RAMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ONLINE
      MPI_master_only write(stdout,'(10x,A)') 'ONLINE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ONLINE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CUBIC_INTERP
      MPI_master_only write(stdout,'(10x,A)') 'CUBIC_INTERP'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CUBIC_INTERP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_MONTH_1DIGIT
      MPI_master_only write(stdout,'(10x,A)') 'BULK_MONTH_1DIGIT'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_MONTH_1DIGIT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_ECUMEV0
      MPI_master_only write(stdout,'(10x,A)') 'BULK_ECUMEV0'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_ECUMEV0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_GUSTINESS
      MPI_master_only write(stdout,'(10x,A)') 'BULK_GUSTINESS'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_GUSTINESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_ECUMEV6
      MPI_master_only write(stdout,'(10x,A)') 'BULK_ECUMEV6'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_ECUMEV6'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BULK_WASP
      MPI_master_only write(stdout,'(10x,A)') 'BULK_WASP'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BULK_WASP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CFB_STRESS
      MPI_master_only write(stdout,'(10x,A)') 'CFB_STRESS'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CFB_STRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CFB_WIND_TRA
      MPI_master_only write(stdout,'(10x,A)') 'CFB_WIND_TRA'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CFB_WIND_TRA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_SFC_BREAK
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_SFC_BREAK'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_SFC_BREAK'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_BREAK_CT93
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_BREAK_CT93'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_BREAK_CT93'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_BREAK_TG86
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_BREAK_TG86'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_BREAK_TG86'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_BREAK_TG86A
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_BREAK_TG86A'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_BREAK_TG86A'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_BREAK_R93
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_BREAK_R93'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_BREAK_R93'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_BREAK_BJ78
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_BREAK_BJ78'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_BREAK_BJ78'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MRL_CEW
      MPI_master_only write(stdout,'(10x,A)') 'MRL_CEW'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MRL_CEW'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WKB_KZ_FILTER
      MPI_master_only write(stdout,'(10x,A)') 'WKB_KZ_FILTER'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WKB_KZ_FILTER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WKB_TIME_FILTER
      MPI_master_only write(stdout,'(10x,A)') 'WKB_TIME_FILTER'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WKB_TIME_FILTER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WKB_ADD_DIFF
      MPI_master_only write(stdout,'(10x,A)') 'WKB_ADD_DIFF'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WKB_ADD_DIFF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WKB_ADD_DIFFRACTION
      MPI_master_only write(stdout,'(10x,A)') 'WKB_ADD_DIFFRACTION'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WKB_ADD_DIFFRACTION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WKB_NUDGING
      MPI_master_only write(stdout,'(10x,A)') 'WKB_NUDGING'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WKB_NUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SHOREFACE
      MPI_master_only write(stdout,'(10x,A)') 'SHOREFACE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SHOREFACE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SANDBAR
      MPI_master_only write(stdout,'(10x,A)') 'SANDBAR'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SANDBAR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BRY_WKB
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BRY_WKB'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BRY_WKB'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef STOKES_DRIFT
      MPI_master_only write(stdout,'(10x,A)') 'STOKES_DRIFT'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='STOKES_DRIFT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_WWAVE
      MPI_master_only write(stdout,'(10x,A)') 'ANA_WWAVE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_WWAVE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef WAVE_IO
      MPI_master_only write(stdout,'(10x,A)') 'WAVE_IO'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_IO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_SKPP
      MPI_master_only write(stdout,'(10x,A)') 'LMD_SKPP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_SKPP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_SKPP2005
      MPI_master_only write(stdout,'(10x,A)') 'LMD_SKPP2005'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_SKPP2005'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_BKPP
      MPI_master_only write(stdout,'(10x,A)') 'LMD_BKPP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_BKPP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_BKPP2005
      MPI_master_only write(stdout,'(10x,A)') 'LMD_BKPP2005'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_BKPP2005'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIOLOGY
      MPI_master_only write(stdout,'(10x,A)') 'BIOLOGY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef PISCES
      MPI_master_only write(stdout,'(10x,A)') 'PISCES'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PISCES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIURNAL_INPUT_SFLX
      MPI_master_only write(stdout,'(10x,A)') 'DIURNAL_INPUT_SFLX'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIURNAL_INPUT_SFLX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DIAGNOSTICS_BIO
      MPI_master_only write(stdout,'(10x,A)') 'DIAGNOSTICS_BIO'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_BIO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_trc_diaadd
      MPI_master_only write(stdout,'(10x,A)') 'key_trc_diaadd'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_trc_diaadd'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BIO_BioEBUS
      MPI_master_only write(stdout,'(10x,A)') 'BIO_BioEBUS'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BIO_BioEBUS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef VAR_CHL_C
      MPI_master_only write(stdout,'(10x,A)') 'VAR_CHL_C'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_CHL_C'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef CARBON
      MPI_master_only write(stdout,'(10x,A)') 'CARBON'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CARBON'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef HYDROGEN_SULFIDE
      MPI_master_only write(stdout,'(10x,A)') 'HYDROGEN_SULFIDE'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HYDROGEN_SULFIDE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BTFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BTFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BTFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BHFLUX
      MPI_master_only write(stdout,'(10x,A)') 'BHFLUX'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BHFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BSFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BWFLUX
      MPI_master_only write(stdout,'(10x,A)') 'BWFLUX'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BWFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LIMIT_BSTRESS
      MPI_master_only write(stdout,'(10x,A)') 'LIMIT_BSTRESS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LIMIT_BSTRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BBL
      MPI_master_only write(stdout,'(10x,A)') 'BBL'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BBL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_BSEDIM
      MPI_master_only write(stdout,'(10x,A)') 'ANA_BSEDIM'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BSEDIM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z0_BL
      MPI_master_only write(stdout,'(10x,A)') 'Z0_BL'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z0_BL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z0_RIP
      MPI_master_only write(stdout,'(10x,A)') 'Z0_RIP'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z0_RIP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef Z0_BIO
      MPI_master_only write(stdout,'(10x,A)') 'Z0_BIO'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z0_BIO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MUSTANG
      MPI_master_only write(stdout,'(10x,A)') 'MUSTANG'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MUSTANG'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SEDIMENT
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SEDIMENT'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SEDIMENT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BED_ARMOR
      MPI_master_only write(stdout,'(10x,A)') 'BED_ARMOR'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BED_ARMOR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BED_HIDEXP
      MPI_master_only write(stdout,'(10x,A)') 'BED_HIDEXP'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BED_HIDEXP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_VANDERA
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_VANDERA'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_VANDERA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_MPM
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_MPM'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_MPM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_WULIN
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_WULIN'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_WULIN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_MARIEU
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_MARIEU'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_MARIEU'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_UP1
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_UP1'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_UP1'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_UP5
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_UP5'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_UP5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BEDLOAD_WENO5
      MPI_master_only write(stdout,'(10x,A)') 'BEDLOAD_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BEDLOAD_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SLOPE_LESSER
      MPI_master_only write(stdout,'(10x,A)') 'SLOPE_LESSER'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SLOPE_LESSER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SLOPE_NEMETH
      MPI_master_only write(stdout,'(10x,A)') 'SLOPE_NEMETH'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SLOPE_NEMETH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SLOPE_KIRWAN
      MPI_master_only write(stdout,'(10x,A)') 'SLOPE_KIRWAN'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SLOPE_KIRWAN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef DUNE
      MPI_master_only write(stdout,'(10x,A)') 'DUNE'
      ie=is + 3
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DUNE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_DUNE
      MPI_master_only write(stdout,'(10x,A)') 'ANA_DUNE'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_DUNE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef SUBSTANCE
      MPI_master_only write(stdout,'(10x,A)') 'SUBSTANCE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SUBSTANCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef USE_CALENDAR
      MPI_master_only write(stdout,'(10x,A)') 'USE_CALENDAR'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='USE_CALENDAR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_noTSdiss_insed
      MPI_master_only write(stdout,'(10x,A)') 'key_noTSdiss_insed'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_noTSdiss_insed'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_nofluxwat_IWS
      MPI_master_only write(stdout,'(10x,A)') 'key_nofluxwat_IWS'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_nofluxwat_IWS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_CROCO
      MPI_master_only write(stdout,'(10x,A)') 'key_CROCO'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_CROCO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_MORPHODYN
      MPI_master_only write(stdout,'(10x,A)') 'ANA_MORPHODYN'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_MORPHODYN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MORPHODYN
      MPI_master_only write(stdout,'(10x,A)') 'MORPHODYN'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MORPHODYN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef MORPHODYN_MUSTANG_byHYDRO
      MPI_master_only write(stdout,'(10x,A)') 'MORPHODYN_MUSTANG_byHYDRO'
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MORPHODYN_MUSTANG_byHYDRO'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2SPECIFIED'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2CHARACT
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2CHARACT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2CHARACT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M2ORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M2ORLANSKI'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_VOLCONS
      MPI_master_only write(stdout,'(10x,A)') 'OBC_VOLCONS'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_VOLCONS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M3SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M3SPECIFIED'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M3SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_M3ORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_M3ORLANSKI'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M3ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_TSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'OBC_TSPECIFIED'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_TSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef OBC_TORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'OBC_TORLANSKI'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_TORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF'
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef key_agrif
      MPI_master_only write(stdout,'(10x,A)') 'key_agrif'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='key_agrif'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_UPDATE_MIX_LOW
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_UPDATE_MIX_LOW'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_MIX_LOW'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_UPDATE_MIX
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_UPDATE_MIX'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_MIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_UPDATE_DECAL
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_UPDATE_DECAL'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_UPDATE_DECAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_CONSERV_VOL
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_CONSERV_VOL'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_CONSERV_VOL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_CONSERV_TRA
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_CONSERV_TRA'
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_CONSERV_TRA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_EAST
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_EAST'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_EAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_WEST
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_WEST'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_NORTH
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_NORTH'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_NORTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_SOUTH
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_SOUTH'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_SOUTH'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_FLUX_BC
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_FLUX_BC'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_FLUX_BC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_M2SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_M2SPECIFIED'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_M2SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_2WAY
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_2WAY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_2WAY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_M3SPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_M3SPECIFIED'
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_M3SPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_TSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_TSPECIFIED'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_TSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_M3ORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_M3ORLANSKI'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_M3ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_TORLANSKI
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_TORLANSKI'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_TORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_WSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_WSPECIFIED'
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_WSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AGRIF_OBC_NBQSPECIFIED
      MPI_master_only write(stdout,'(10x,A)') 'AGRIF_OBC_NBQSPECIFIED'
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AGRIF_OBC_NBQSPECIFIED'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef FILLVAL
      MPI_master_only write(stdout,'(10x,A)') 'FILLVAL'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FILLVAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NC4PAR
      MPI_master_only write(stdout,'(10x,A)') 'NC4PAR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NC4PAR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef NF_CLOBBER
      MPI_master_only write(stdout,'(10x,A)') 'NF_CLOBBER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_CLOBBER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef AVERAGES_K
      MPI_master_only write(stdout,'(10x,A)') 'AVERAGES_K'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AVERAGES_K'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_STFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_STFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_STFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SSFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ANA_SRFLUX
      MPI_master_only write(stdout,'(10x,A)') 'ANA_SRFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3NUDGING
      MPI_master_only write(stdout,'(10x,A)') 'M3NUDGING'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3NUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef TNUDGING
      MPI_master_only write(stdout,'(10x,A)') 'TNUDGING'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TNUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef ROBUST_DIAG
      MPI_master_only write(stdout,'(10x,A)') 'ROBUST_DIAG'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ROBUST_DIAG'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef M3_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'M3_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef T_FRC_BRY
      MPI_master_only write(stdout,'(10x,A)') 'T_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='T_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BODYFORCE
      MPI_master_only write(stdout,'(10x,A)') 'BODYFORCE'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BODYFORCE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef BVF_MIXING
      MPI_master_only write(stdout,'(10x,A)') 'BVF_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BVF_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_MIXING
      MPI_master_only write(stdout,'(10x,A)') 'LMD_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_RIMIX
      MPI_master_only write(stdout,'(10x,A)') 'LMD_RIMIX'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_RIMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
#ifdef LMD_CONVEC
      MPI_master_only write(stdout,'(10x,A)') 'LMD_CONVEC'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LMD_CONVEC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
#endif
      MPI_master_only write(stdout,'(/)')
      if (iexample.eq.0) then
        MPI_master_only write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
        MPI_master_only write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  MPI_master_only write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
