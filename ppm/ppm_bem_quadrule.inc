      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_bem_quadrule
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Sets up various rules for the quadrature (cubature)
      !                 over a unit triangle
      !
      !  Input        : rule    (I) : the rule to be returned
      !
      !  Input/output : 
      !
      !  Output       : qp      (F) : quadrature points
      !                 qw      (F) : quadrature weights
      !                 info    (I) : return status. 0 on success.
      !
      !  Remarks      : qp is a 2-dimensional array containing the (s,t)-
      !                 coordinates of quadrature points in the unit triangle
      !
      !                  * (0,1)
      !                  |\
      !                  | \
      !                  t  \
      !                  |   \
      !                  |    \
      !                  *--s--* (1,0)
      !                 (0,0)
      !
      !                 e.g. qp(1,1) contains the s-coordinate of the first
      !                 quadrature points and qp(2,1) the t-coordinate. The
      !                 array qp and qw must be big enough to hold the
      !                 quadrature points indicated by the rule.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_bem_quadrule.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.6  2004/10/01 16:08:56  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/07/27 07:56:01  ivos
      !  Corrected typo in CALL statement.
      !
      !  Revision 1.4  2004/07/27 07:06:48  ivos
      !  Renamed parameters after moving their definition to ppm_param.h.
      !
      !  Revision 1.3  2004/07/26 11:46:53  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:46:35  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.1  2004/07/16 08:33:11  oingo
      !  Initial release. Not tested
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_bem_quadrule_s(qp,qw,rule,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_bem_quadrule_d(qp,qw,rule,info)
#endif
      
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_bem_quadrule_npoints
      !-------------------------------------------------------------------------
      !  Type kind
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), POINTER       :: qp
      REAL(MK), DIMENSION(:)  , POINTER       :: qw
      INTEGER                 , INTENT(IN   ) :: rule
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: temp
      INTEGER                :: nqp
      REAL(MK)               :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_bem_quadrule',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_bem_quadrule_npoints(rule,nqp,info)
          IF (nqp .EQ. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_bem_quadrule',   &
     &             'Unknown quadrature rule',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(qp,2) .LT. nqp) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_argument,'ppm_bem_quadrule',   &
     &             'QP array not big enough',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(qw,1) .LT. nqp) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_argument,'ppm_bem_quadrule',   &
     &             'QW array not big enough',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Fill in the quadrature points according to the rule
      !-------------------------------------------------------------------------
      IF (rule .EQ. ppm_param_bem_quadrule_center) THEN
          qp(1:2,1) = (/ 1.0/3.0, 1.0/3.0 /)
          qw(1)     = 1.0
      ELSE IF (rule .EQ. ppm_param_bem_quadrule_nodes) THEN
          qp(1:2,1) = (/ 0.0, 0.0 /)
          qp(1:2,2) = (/ 1.0, 0.0 /)
          qp(1:2,3) = (/ 0.0, 1.0 /)
          qw(:)     = 1.0 / 3.0
      ELSE IF (rule .EQ. ppm_param_bem_quadrule_edges) THEN
          qp(1:2,1) = (/ 0.0, 0.5 /)
          qp(1:2,2) = (/ 0.5, 0.5 /)
          qp(1:2,3) = (/ 0.5, 0.0 /)
          qw(:)     = 1.0 / 3.0
      ELSE IF (rule .EQ. ppm_param_bem_quadrule_cne) THEN
          qp(1:2,1) = (/ 0.0, 0.0 /)
          qp(1:2,2) = (/ 1.0, 0.0 /)
          qp(1:2,3) = (/ 0.0, 1.0 /)
          qp(1:2,4) = (/ 0.0, 0.5 /)
          qp(1:2,5) = (/ 0.5, 0.5 /)
          qp(1:2,6) = (/ 0.5, 0.0 /)
          qp(1:2,7) = (/ 1.0/3.0, 1.0/3.0 /)
          qw(1:3)   = 1.0 / 20.0
          qw(4:6)   = 2.0 / 15.0
          qw(7)     = 27.0 / 60.0
      ELSE IF (rule .EQ. ppm_param_bem_quadrule_stroud) THEN  ! Stroud T2:5-1
          temp = SQRT(15.0_MK)
          qp(1:2,1) = (/ 1.0/3.0, 1.0/3.0 /)
          qp(1:2,2) = (/ (9.0 + 2.0*temp) / 21.0, (6.0 -     temp) / 21.0 /)
          qp(1:2,3) = (/ (6.0 -     temp) / 21.0, (9.0 + 2.0*temp) / 21.0 /)
          qp(1:2,4) = (/ (6.0 -     temp) / 21.0, (6.0 -     temp) / 21.0 /)
          qp(1:2,5) = (/ (9.0 - 2.0*temp) / 21.0, (6.0 +     temp) / 21.0 /)
          qp(1:2,6) = (/ (6.0 +     temp) / 21.0, (9.0 - 2.0*temp) / 21.0 /)
          qp(1:2,7) = (/ (6.0 +     temp) / 21.0, (6.0 +     temp) / 21.0 /)
          qw(1)     = 9.0 / 40.0
          qw(2:4)   = (155.0 - temp) / 1200.0
          qw(5:7)   = (155.0 + temp) / 1200.0
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_bem_quadrule',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_bem_quadrule_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_bem_quadrule_d
#endif
