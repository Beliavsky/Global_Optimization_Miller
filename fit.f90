!   PROGRAM NAME    - FIT
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-27  Time: 13:16:33


!-----------------------------------------------------------------------

!   LATEST REVISION - DECEMBER 27, 1986

!   PURPOSE     - PARAMETER ESTIMATION OF RESPIRATORY MODELS
!                 BY A GLOBAL OPTIMIZATION ALGORITHM

!   REQUIRED ROUTINES   - GLOBAL, LOCAL, URDMN, FUNCT, FUN

!-----------------------------------------------------------------------

MODULE fit_common
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /fu/ om(50), rezz(50), irel, zre(50), zim(50)

REAL (dp), SAVE :: om(50), rezz(50), zre(50), zim(50)
INTEGER, SAVE   :: irel

END MODULE fit_common

PROGRAM fit
USE fit_common
USE global_minimum
IMPLICIT NONE

REAL (dp)             :: x0(15,20), f00(20), MIN(18), MAX(18),   &
                         df, do, f0, o0, o1, tpi
CHARACTER (LEN=80)    :: label
INTEGER               :: i, in, ipr, k, m, nc, ng, nparm, npt, nsampl, nsig
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, onep5 = 1.5_dp,  &
                         eight = 8.0_dp, tausend = 1000._dp
INTEGER, ALLOCATABLE  :: seed(:)
character (len=*), parameter :: input_file = "input.txt", output_file = "output.txt"
INTERFACE
  SUBROUTINE ladder(x, f, f0, df, o0, DO, npt, ipr, np)
    USE fit_common
    IMPLICIT NONE
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN)      :: f
    REAL (dp), INTENT(IN)      :: f0
    REAL (dp), INTENT(IN)      :: df
    REAL (dp), INTENT(IN)      :: o0
    REAL (dp), INTENT(IN OUT)  :: DO
    INTEGER, INTENT(IN OUT)    :: npt
    INTEGER, INTENT(IN OUT)    :: ipr
    INTEGER, INTENT(IN OUT)    :: np
  END SUBROUTINE ladder
END INTERFACE

tpi = eight*ATAN(one)
in = 5
ipr = 20
OPEN(5, FILE=input_file, action = "read")
OPEN(ipr, FILE=output_file, action = "write")
READ(in, 901) label
901 FORMAT(A)
WRITE(ipr, 902) label
WRITE(*, 902) label
902 FORMAT(' ', A/)
nparm = 5
READ(in, *) MIN(1), MAX(1), MIN(2), MAX(2), MIN(3), MAX(3), MIN(4), MAX(4), &
            MIN(5), MAX(5)
920 FORMAT(' ', f9.4, ' ', f9.4)
WRITE(ipr, 920) MIN(1), MAX(1), MIN(2), MAX(2), MIN(3), MAX(3),  &
                MIN(4), MAX(4), MIN(5), MAX(5)
WRITE(*, 920) MIN(1), MAX(1), MIN(2), MAX(2), MIN(3), MAX(3),  &
              MIN(4), MAX(4), MIN(5), MAX(5)
WRITE(ipr, 905)
WRITE(*, 905)
905 FORMAT(// ' RUN PARAMETERS'/)
READ(in, *) irel, nsampl, ng, nsig
WRITE(ipr, 922) irel, nsampl, ng, nsig
WRITE(*, 922) irel, nsampl, ng, nsig
922 FORMAT(' ', i2, ' ', i4, ' ', i2, ' ', i1 /)
WRITE(ipr, 906)
WRITE(*, 906)
906 FORMAT(' SAMPLE'/)
m = 6
do i=1,m
   READ(in, *) om(i), zre(i), zim(i)
   WRITE(ipr, 925) om(i), zre(i), zim(i)
   WRITE(*, 925) om(i), zre(i), zim(i)
   925 FORMAT(' ', f6.3, ' ', f7.2, ' ', f8.2)
end do
WRITE(ipr, 926)
WRITE(*, 926)
926 FORMAT(//)
o0 = tausend
o1 = zero
DO = (om(2) - om(1))*tpi
DO  i=1,m
  om(i) = om(i)*tpi
  IF (om(i) < o0) o0 = om(i)
  IF (om(i) > o1) o1 = om(i)
  IF (i <= 1) CYCLE
  IF (om(i) - om(i-1) < DO) DO = om(i) - om(i-1)
END DO

! Set the random number seed

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)
WRITE(ipr, '(a / (7i11) )') ' Random number seed(s): ', seed
WRITE(ipr, * )

f0 = o0/tpi
df = DO/tpi
npt = INT((o1-o0)/DO + onep5)
DO  i=1,m
  rezz(i) = SQRT(zre(i)*zre(i) + zim(i)*zim(i))
END DO
CALL global(MIN, MAX, nparm, m, nsampl, ng, ipr, nsig, x0, nc, f00)
DO  i=1,nc
  CALL ladder(x0(1:,i), f00(i), f0, df, o0, DO, npt, ipr, nparm)
END DO
print*,"wrote to " // trim(output_file)
stop "end of program fit"
9 STOP
END PROGRAM fit
    

SUBROUTINE funct(x, value, nparm, m)

USE fit_common
IMPLICIT NONE

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: value
INTEGER, INTENT(IN)     :: nparm
INTEGER, INTENT(IN)     :: m

REAL (dp) :: f(100), zimi, zrei
INTEGER   :: i, j, kk, mm

DO  kk=1,m
  zrei = x(1) + (x(2)/(om(kk)**x(3)))
  zimi = om(kk)*x(4) - (x(5)/(om(kk)**x(3)))
  zrei = zre(kk) - zrei
  zimi = zim(kk) - zimi
  j = kk*2
  i = j-1
  IF (irel /= 0) GO TO 100
  f(i) = zrei/rezz(kk)
  f(j) = zimi/rezz(kk)
  CYCLE
  100 f(i) = zrei
  f(j) = zimi
END DO

mm = m+m
value = SUM( f(1:mm)**2 )
value = SQRT(value/m)
RETURN
END SUBROUTINE funct


SUBROUTINE ladder(x, f, f0, df, o0, DO, npt, ipr, np)

USE fit_common
IMPLICIT NONE

REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: f
REAL (dp), INTENT(IN)      :: f0
REAL (dp), INTENT(IN)      :: df
REAL (dp), INTENT(IN)      :: o0
REAL (dp), INTENT(IN OUT)  :: DO
INTEGER, INTENT(IN OUT)    :: npt
INTEGER, INTENT(IN OUT)    :: ipr
INTEGER, INTENT(IN OUT)    :: np

REAL (dp) :: zbr(50), zbi(50), fok, oj, tr1, tr2, tr3, tr4, tr5, tr6, tr7
INTEGER   :: i, kk
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, eight = 8.0_dp,  &
                        ts = 360._dp

oj = o0 - DO
fok = ts/(eight*ATAN(one))
WRITE(ipr, 901) f, x(1:np)
901 FORMAT(/////' ', g14.8, 3(/'    ', 5(g14.8, ' ')))
WRITE(ipr, 903)
903 FORMAT(//'  FREQ    REAL       IMAG       ABS       PHASE    DRE ',  &
             '   DIM   DABS    FREQ'/)
DO  kk=1,npt
  oj = oj + DO
  zbr(kk) = x(1) + (x(2)/(oj**x(3)))
  zbi(kk) = oj*x(4) - (x(5)/(oj**x(3)))
END DO
oj = f0 - df
DO  i=1,npt
  oj = oj + df
  tr1 = zbr(i)
  tr2 = zbi(i)
  tr3 = SQRT(tr1*tr1 + tr2*tr2)
  tr4 = ATAN2(tr2, tr1)*fok
  IF (zre(i) /= zero) tr5 = (zre(i)-tr1)/ABS(zre(i))
  IF (zre(i) == zero) tr5 = zre(i)-tr1
  IF (zim(i) /= zero) tr6 = (zim(i)-tr2)/ABS(zim(i))
  IF (zim(i) == zero) tr6 = zim(i)-tr2
  tr7 = SQRT(tr5*tr5 + tr6*tr6)
  WRITE(ipr, 904) oj, tr1, tr2, tr3, tr4, tr5, tr6, tr7, oj
  904 FORMAT(' ', f5.2, ' ', 3(' ', g10.4), 4(' ', f6.2), '   ', f5.2)
END DO
RETURN
END SUBROUTINE ladder
