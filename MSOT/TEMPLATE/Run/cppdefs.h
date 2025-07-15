! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
/*
   This is "cppdefs.h": MODEL CONFIGURATION FILE
   ==== == ============ ===== ============= ====
*/
#define Argentina
#undef Argentina_Anidado1
# if defined Argentina
/*
!====================================================================
!    Modelo para la Simulacion de Ondas de Tormenta (MSOT)
!====================================================================
!
!----------------------
! BASIC OPTIONS
!----------------------
!
*/
 
                      /* Parallelization */
# define  OPENMP
# undef  MPI
                      /* Open Boundary Conditions */
/* Viejo MSOT
# define START_DATE
# undef EXACT_RESTART
# define TIDES2D
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH


FORZAR NIVEL Y CORRIENTES*/
# define TIDES         /* clave maestra */
# define SSH_TIDES     /* forzar elevación */
# define UV_TIDES      /* forzar corrientes */
# define POT_TIDES     /* forzar potencial gravitatorio (opcional) */
# undef TIDES2D      /* quitamos la versión analítica */
# undef ANA_BRY      /* no más frontera analítica */
# define FRC_BRY       /* frontera desde archivos netCDF */
# define Z_FRC_BRY     /* leer amplitud/fase de zeta */
# define OBC_M2CHARACT /* método característico para UV */
# undef OBC_M2ORLANSKI
# define OBC_WEST
# define OBC_EAST
# define OBC_NORTH
# define OBC_SOUTH
                      /* Grid configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
# define WET_DRY      /* que se vacie la playita */
                      /* Model dynamics */
# undef SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Surface Forcing */
# define ATM_PRESS
                      /* Lateral Forcing */
# undef CLIMATOLOGY
# define QUADRATIC_BOTTOM_FRICTION_COEF
                      /* Point Sources - Rivers */
# define ANA_PSOURCE
# define PSOURCE
# define PSOURCE_NCFILE
                      /* Open Boundary Conditions */
# ifdef TIDES2D
#   define ANA_BRY    /* OJO: modificación de Mati en clm_tifes.F */
                      /* Las corrientes se infieren internamente */ 
                      /* pero amplitudes y fases se leen de $jobID_frc.nc */
#   define Z_FRC_BRY
#   define SSH_TIDES
#   undef TIDERAMP
#   undef SWACS
# endif

#elif defined Argentina_Anidado1

/*
!====================================================================
! Modelo para la Simulacion de Ondas de Tormenta (MSOT) -- Anidado
!====================================================================
*/
                      /* Parallelization */
# define  OPENMP
# undef  MPI
                      /* Open Boundary Conditions */
# define START_DATE
# undef EXACT_RESTART
# undef TIDES2D
# define OBC_EAST
# define OBC_WEST
# define OBC_NORTH
# define OBC_SOUTH
                      /* Grid configuration */
# define CURVGRID
# define SPHERICAL
# define MASKING
# define WET_DRY      /* que se vacie la playita */
                      /* Model dynamics */
# undef SOLVE3D
# define UV_COR
# define UV_ADV
                      /* Surface Forcing */
# define ATM_PRESS
                      /* Lateral Forcing */
# undef CLIMATOLOGY
# define QUADRATIC_BOTTOM_FRICTION_COEF
                      /* Point Sources - Rivers */
# define ANA_PSOURCE
# define PSOURCE
# define PSOURCE_NCFILE
                      /* Open Boundary Conditions */
# define  FRC_BRY
# ifdef FRC_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
# endif

#endif /* END OF CONFIGURATION CHOICE */

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
