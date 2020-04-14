

MODULE basic_types

IMPLICIT NONE

! define basic types
INTEGER, PARAMETER  :: int1  = SELECTED_INT_KIND (2)
INTEGER, PARAMETER  :: int2  = SELECTED_INT_KIND (4)
INTEGER, PARAMETER  :: int4  = SELECTED_INT_KIND (8)
INTEGER, PARAMETER  :: real4 = SELECTED_REAL_KIND(P=6)
INTEGER, PARAMETER  :: real8 = SELECTED_REAL_KIND(P=15)


END MODULE basic_types

!!=============================================================================================================================
!!=============================================================================================================================


MODULE io_files

IMPLICIT NONE

! global io file names
CHARACTER(LEN=20) :: input_file, psy_results_file, GiD_mesh_file, GiD_results_file, error_file, extra_results_file

! io units
INTEGER, PARAMETER :: default_input_unit  =  5   ! keyboard
INTEGER, PARAMETER :: default_output_unit =  6   ! screen
INTEGER, PARAMETER :: input_file_unit     = 11   ! input file (.inp)
INTEGER, PARAMETER :: psy_results_unit    = 12   ! psy results file (.rst)
INTEGER, PARAMETER :: GiD_mesh_unit       = 13   ! GiD mesh file (.post.msh)
INTEGER, PARAMETER :: GiD_results_unit    = 14   ! GiD results file (post.res)
INTEGER, PARAMETER :: error_unit          = 21   ! error file (.err)
INTEGER, PARAMETER :: extra_results_unit  = 22   ! extra results file (.xrs) (for printing additional results)


END MODULE io_files

!!=============================================================================================================================
!!=============================================================================================================================


MODULE clock

USE basic_types

IMPLICIT NONE

INTEGER(int4) :: count_rate, count_max

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE start_clock(t)
INTEGER(int4) :: t
CALL SYSTEM_CLOCK(t,count_rate,count_max)
END SUBROUTINE start_clock
!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE stop_clock(t,unit)
INTEGER(int4), INTENT(INOUT) :: t
REAL(KIND=real8) :: tt
INTEGER(int4) :: temp
CHARACTER(LEN=*),INTENT(INOUT) :: unit
CALL SYSTEM_CLOCK(temp,count_rate,count_max)
tt = REAL(temp-t,real8)/REAL(count_rate,real8)
tt = 1000.0_real8*tt
unit=" ms"
t=INT( tt )
END SUBROUTINE stop_clock
!!---------------------------------------------------------------------------------------------------------------------


END MODULE clock

!!=============================================================================================================================
!!=============================================================================================================================


MODULE tensor_class

USE basic_types

IMPLICIT NONE

INTERFACE OPERATOR(.vector.)
  MODULE PROCEDURE vector_product
END INTERFACE

INTERFACE OPERATOR(.det.)
  MODULE PROCEDURE det
END INTERFACE

INTERFACE OPERATOR(.inv.)
  MODULE PROCEDURE inv
END INTERFACE

INTERFACE OPERATOR(.skew.)
  MODULE PROCEDURE skew
END INTERFACE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION inv(a)

! dymmy arguments
REAL(KIND=real8), DIMENSION(3,3), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(3,3)             :: inv

! local variables
REAL(KIND=real8)            :: d
REAL(KIND=real8), PARAMETER :: one=1.0_real8

d = one/det(a)

inv(1,1) = d*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) 
inv(1,2) = d*(a(1,3)*a(3,2)-a(1,2)*a(3,3))
inv(1,3) = d*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
inv(2,1) = d*(a(2,3)*a(3,1)-a(2,1)*a(3,3))
inv(2,2) = d*(a(1,1)*a(3,3)-a(1,3)*a(3,1))
inv(2,3) = d*(a(1,3)*a(2,1)-a(1,1)*a(2,3))
inv(3,1) = d*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
inv(3,2) = d*(a(1,2)*a(3,1)-a(1,1)*a(3,2))
inv(3,3) = d*(a(1,1)*a(2,2)-a(1,2)*a(2,1))

END FUNCTION inv
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION skew(a)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(3,3)           :: skew

skew(1,1) = 0.0_real8
skew(2,2) = 0.0_real8
skew(3,3) = 0.0_real8

skew(3,2) = a(1)
skew(1,3) = a(2)
skew(2,1) = a(3)

skew(2,3) = -a(1)
skew(3,1) = -a(2)
skew(1,2) = -a(3)

END FUNCTION skew
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION det(a)

! dymmy arguments
REAL(KIND=real8), DIMENSION(3,3), INTENT(IN) :: a
REAL(KIND=real8)                             :: det

det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
      a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
      a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

END FUNCTION det
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION build(a,b,c)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a, b, c
REAL(KIND=real8), DIMENSION(3,3)           :: build

build(1:3,1) = a ; build(1:3,2) = b ; build(1:3,3) = c

END FUNCTION build
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION vector_product (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a,b
REAL(KIND=real8), DIMENSION(3)             :: vector_product

vector_product(1) = a(2)*b(3) - a(3)*b(2)
vector_product(2) = a(3)*b(1) - a(1)*b(3)
vector_product(3) = a(1)*b(2) - a(2)*b(1)

END FUNCTION vector_product
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION normalize(a)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(3)             :: normalize

! local variables
REAL(KIND=real8) :: norm_a

norm_a = SQRT(DOT_PRODUCT(a,a))
normalize = (1.0_real8/norm_a)*a

END FUNCTION normalize
!!---------------------------------------------------------------------------------------------------------------------


END MODULE tensor_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE rotation_class

USE basic_types
USE tensor_class

IMPLICIT NONE

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

FUNCTION rotation_tensor(a)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(3,3)           :: rotation_tensor

! local variables
REAL(KIND=real8), PARAMETER :: half=0.5_real8, one=1.0_real8, four=4.0_real8
REAL(KIND=real8)            :: norm_a, h

norm_a = DOT_PRODUCT(a,a)
h = four/(four + norm_a)
rotation_tensor = h*( skew(a) + half*(MATMUL(skew(a),skew(a))) )
rotation_tensor(1,1) = rotation_tensor(1,1) + one
rotation_tensor(2,2) = rotation_tensor(2,2) + one
rotation_tensor(3,3) = rotation_tensor(3,3) + one

END FUNCTION rotation_tensor

!!---------------------------------------------------------------------------------------------------------------------

FUNCTION xsi_tensor(a)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(3,3)           :: xsi_tensor

! local variables
REAL(KIND=real8), PARAMETER :: half=0.5_real8, one=1.0_real8, four=4.0_real8
REAL(KIND=real8)            :: norm_a, h

norm_a = DOT_PRODUCT(a,a)
h = four/(four + norm_a)
xsi_tensor = half*skew(a)
xsi_tensor(1,1) = xsi_tensor(1,1) + one
xsi_tensor(2,2) = xsi_tensor(2,2) + one
xsi_tensor(3,3) = xsi_tensor(3,3) + one
xsi_tensor = h*xsi_tensor

END FUNCTION xsi_tensor

!!---------------------------------------------------------------------------------------------------------------------

FUNCTION w_tensor(a,t)

! dummy arguments
REAL(KIND=real8), DIMENSION(3), INTENT(IN) :: a, t
REAL(KIND=real8), DIMENSION(3,3)           :: w_tensor

! local variables
REAL(KIND=real8), DIMENSION(3,3) :: Xsi, aux1
REAL(KIND=real8), DIMENSION(3)   :: aux2
REAL(KIND=real8), PARAMETER :: half=0.5_real8, four=4.0_real8
REAL(KIND=real8)            :: norm_a, h
INTEGER(KIND=int4)          :: i, j

norm_a = DOT_PRODUCT(a,a)
h = four/(four + norm_a)
Xsi = xsi_tensor(a)
aux2 = MATMUL(TRANSPOSE(Xsi),t)

DO i=1,3
  DO j=1,3
    aux1(i,j) = aux2(i)*a(j)
  END DO
END DO

w_tensor = half*h*(skew(t) - aux1)

END FUNCTION w_tensor

!!---------------------------------------------------------------------------------------------------------------------


END MODULE rotation_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE spectral_decomposition_class

USE basic_types
USE tensor_class

IMPLICIT NONE

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE eigen_values(matrix,evalues)

! dummy variables
REAL(KIND=real8), DIMENSION(3,3), INTENT(IN)  :: matrix
REAL(KIND=real8), DIMENSION(3),   INTENT(OUT) :: evalues

! local variables
REAL(KIND=real8)               :: j2, j3, p, theta, pi, sm
REAL(KIND=real8), DIMENSION(3) :: d, o, td
LOGICAL, DIMENSION(3)          :: mask

! build diagonal vector 
d(1) = matrix(1,1)
d(2) = matrix(2,2)
d(3) = matrix(3,3)

! build out-of-diagonal vector
o(3) = matrix(1,2)
o(2) = matrix(1,3)
o(1) = matrix(2,3)

! compute tensor invariants
sm = SUM(d)/3.0d0
td = d-sm
j2 = ((d(1)-d(2))**2+(d(2)-d(3))**2+(d(3)-d(1))**2)/6.0d0+SUM(o**2)
j3 = PRODUCT(td)+2.0d0*PRODUCT(o)-SUM(o*o*td)

! compute eigenvalues in decreasing order
IF (j2>1.0d-10) THEN
  p=2.0d0*SQRT(j2/3.0d0)
  pi=4.0d0*ATAN(1.0d0)
  theta=(pi/2.0d0-ASIN(3.0d0*j3/(p*j2)))/3.0d0
  mask=.TRUE.
  evalues(MAXLOC(d,mask))=sm+p*COS(theta)
  mask(MAXLOC(d,mask))=.FALSE.
  evalues(MAXLOC(d,mask))=sm+p*COS(theta-2.0d0*pi/3.0d0)
  mask(MAXLOC(d,mask))=.FALSE.
  evalues(MAXLOC(d,mask))=sm+p*COS(theta-4.0d0*pi/3.0d0)
ELSE
  evalues=sm
END IF

END SUBROUTINE eigen_values

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE eigen_vectors(matrix,evalues,evectors)

! dummy variables
REAL(KIND=real8), DIMENSION(3,3), INTENT(IN)  :: matrix
REAL(KIND=real8), DIMENSION(3,3), INTENT(OUT) :: evectors
REAL(KIND=real8), DIMENSION(3),   INTENT(IN)  :: evalues

! local variables
REAL(KIND=real8), DIMENSION(3,3) :: AA, BB
REAL(KIND=real8), DIMENSION(3)   :: d, o, alfa, ff
REAL(KIND=real8)                 :: dd, a2
INTEGER(KIND=int4)               :: i


! build diagonal vector 
d(1) = matrix(1,1)
d(2) = matrix(2,2)
d(3) = matrix(3,3)

! build out-of-diagonal vector
o(3) = matrix(1,2)
o(2) = matrix(1,3)
o(1) = matrix(2,3)

! compute alpha by Newton method (Rodrigues parameters)
CALL matrix_B(d,o,evalues,BB)
alfa=0.0_real8
DO i=0,100
  CALL matrix_A(d,o,evalues,alfa,AA)
  ff = MATMUL(0.5_real8*AA+BB,alfa) - o
  IF (SQRT(DOT_PRODUCT(ff,ff))<1.0d-10) EXIT
  dd = det(AA+BB)
  IF(ABS(dd)>1.0d-30)THEN
    alfa = alfa - (1.0d0/dd)*MATMUL(cof(AA+BB),ff)
  ELSE
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(alfa)
  END IF
END DO

! compute eigenvectors
a2=DOT_PRODUCT(alfa,alfa)
evectors(1,1) = 4.0d0-a2     +2.0d0*alfa(1)*alfa(1)
evectors(2,1) = 4.0d0*alfa(3)+2.0d0*alfa(2)*alfa(1)
evectors(3,1) =-4.0d0*alfa(2)+2.0d0*alfa(3)*alfa(1)
evectors(1,2) =-4.0d0*alfa(3)+2.0d0*alfa(2)*alfa(1)
evectors(2,2) = 4.0d0-a2     +2.0d0*alfa(2)*alfa(2)
evectors(3,2) = 4.0d0*alfa(1)+2.0d0*alfa(2)*alfa(3)
evectors(1,3) = 4.0d0*alfa(2)+2.0d0*alfa(1)*alfa(3)
evectors(2,3) =-4.0d0*alfa(1)+2.0d0*alfa(2)*alfa(3)
evectors(3,3) = 4.0d0-a2     +2.0d0*alfa(3)*alfa(3)
evectors = evectors/(4.0d0+a2)

END SUBROUTINE eigen_vectors

!!---------------------------------------------------------------------------------------------------------------------

FUNCTION cof(a)

! dummy arguments
REAL(KIND=real8), DIMENSION(3,3), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(3,3) :: cof

cof(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2) 
cof(1,2) = a(1,3)*a(3,2)-a(1,2)*a(3,3)
cof(1,3) = a(1,2)*a(2,3)-a(1,3)*a(2,2)
cof(2,1) = a(2,3)*a(3,1)-a(2,1)*a(3,3)
cof(2,2) = a(1,1)*a(3,3)-a(1,3)*a(3,1)
cof(2,3) = a(1,3)*a(2,1)-a(1,1)*a(2,3)
cof(3,1) = a(2,1)*a(3,2)-a(2,2)*a(3,1)
cof(3,2) = a(1,2)*a(3,1)-a(1,1)*a(3,2)
cof(3,3) = a(1,1)*a(2,2)-a(1,2)*a(2,1)

END FUNCTION cof

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE matrix_A(D,o,av,alfa,AA)

! dummy arguments
REAL(KIND=real8), DIMENSION(3),   INTENT(IN)  :: D, o, av, alfa
REAL(KIND=real8), DIMENSION(3,3), INTENT(OUT) :: AA

AA(1,1) = 2.0d0*o(1)*alfa(1)-o(2)*alfa(2)-o(3)*alfa(3)
AA(2,2) = 2.0d0*o(2)*alfa(2)-o(3)*alfa(3)-o(1)*alfa(1)
AA(3,3) = 2.0d0*o(3)*alfa(3)-o(1)*alfa(1)-o(2)*alfa(2)
AA(1,2) = (D(1)-av(1))*alfa(3)-o(2)*alfa(1)
AA(2,1) = (D(2)-av(2))*alfa(3)-o(1)*alfa(2)
AA(1,3) = (D(1)-av(1))*alfa(2)-o(3)*alfa(1)
AA(3,1) = (D(3)-av(3))*alfa(2)-o(1)*alfa(3)
AA(2,3) = (D(2)-av(2))*alfa(1)-o(3)*alfa(2)
AA(3,2) = (D(3)-av(3))*alfa(1)-o(2)*alfa(3)
AA = 0.250d0*AA

END SUBROUTINE matrix_A

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE matrix_B(D,o,av,BB)

! dummy arguments
REAL(KIND=real8), DIMENSION(3),   INTENT(IN)  :: D, o, av
REAL(KIND=real8), DIMENSION(3,3), INTENT(OUT) :: BB

BB(1,1) = D(2)-D(3)+av(2)-av(3)
BB(2,2) = D(3)-D(1)+av(3)-av(1)
BB(3,3) = D(1)-D(2)+av(1)-av(2)
BB(1,2) = -o(3)
BB(2,1) =  o(3)
BB(1,3) =  o(2)
BB(3,1) = -o(2)
BB(2,3) = -o(1)
BB(3,2) =  o(1)
BB=0.50d0*BB

END SUBROUTINE matrix_B

!!---------------------------------------------------------------------------------------------------------------------


END MODULE spectral_decomposition_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE jacobi_class

USE basic_types
USE tensor_class
USE io_files

IMPLICIT NONE

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------
subroutine jacobi(matrix,dd,vv)

! calcula autovalores e autovetores de um tensor simétrico 3x3

real(KIND=real8), dimension(3,3), intent(in)  :: matrix 
real(KIND=real8), dimension(3),   intent(out) :: dd  
real(KIND=real8), dimension(3,3), intent(out) :: vv

! variáveis locais
real(KIND=real8)                 :: theta,c,s,ref
real(KIND=real8), dimension(3,3) :: T,Q
real(KIND=real8), dimension(3)   :: d,o,oo
integer(KIND=int4)               :: i,j,k,iter
integer(KIND=int4), dimension(1) :: ii

! permutação cíclica de índices
integer(KIND=int4), dimension(3), parameter :: jj=(/2,3,1/), kk=(/3,1,2/)

! build diagonal vector 
d(1) = matrix(1,1)
d(2) = matrix(2,2)
d(3) = matrix(3,3)

! build out-of-diagonal vector
o(3) = matrix(1,2)
o(2) = matrix(1,3)
o(1) = matrix(2,3)

! inicialização dos autovetores
vv=0.0d0
vv(1,1)=1.0d0
vv(2,2)=1.0d0
vv(3,3)=1.0d0

! monta matriz T
T(1,1)=d(1); T(2,2)=d(2); T(3,3)=d(3)
T(1,2)=o(3); T(2,1)=o(3)
T(1,3)=o(2); T(3,1)=o(2)
T(2,3)=o(1); T(3,2)=o(1)

! valor de referência
ref=max(1.0d-12*maxval(abs(o)),1.0d-30)

oo=o
iter=0

! iterações de jacobi
do iter=1,20

  if (maxval(abs(oo))<ref) exit

  ii = maxloc(abs(oo))
  i=ii(1); k=kk(i); j=jj(i)

  theta = 0.5d0*atan2(2.0d0*oo(i),T(k,k)-T(j,j))

  c = cos(theta); s = sin(theta)

  Q      = 0.0d0
  Q(i,i) = 1.0d0
  Q(j,j) = c
  Q(k,k) = c
  Q(j,k) = s
  Q(k,j) =-s

  T=matmul(transpose(Q),matmul(T,Q))
  vv=matmul(vv,Q)

  oo(3)=T(1,2); oo(2)=T(1,3); oo(1)=T(2,3)

enddo

! output dos autovalores
dd(1)=T(1,1); dd(2)=T(2,2); dd(3)=T(3,3)

end subroutine jacobi

!!---------------------------------------------------------------------------------------------------------------------


END MODULE jacobi_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE real8_vector_class 

USE basic_types
USE io_files

IMPLICIT NONE

TYPE real8_vector
  REAL(KIND=real8), DIMENSION(:), POINTER :: r
END TYPE real8_vector


INTERFACE norm
  MODULE PROCEDURE real8_norm
END INTERFACE

INTERFACE allocate
  MODULE PROCEDURE alloc_real8_vector
END INTERFACE

INTERFACE deallocate
  MODULE PROCEDURE dealloc_real8_vector
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE assign_real8_vector
END INTERFACE 

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

FUNCTION real8_norm(v)

TYPE(real8_vector) :: v
REAL(KIND=real8)   :: real8_norm

real8_norm=DOT_PRODUCT(v%r,v%r)
real8_norm=SQRT(real8_norm)

END FUNCTION real8_norm

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE alloc_real8_vector(v,n)

! dummy variables
TYPE(real8_vector), INTENT(INOUT) :: v
INTEGER(KIND=int4), INTENT(IN)   :: n


INTEGER(KIND=int4)               :: error


ALLOCATE(v%r(n),STAT=error)
IF(error /= 0) THEN
  WRITE(error_unit,*) "allocation failure of vector"
  STOP
END IF


END SUBROUTINE alloc_real8_vector

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE dealloc_real8_vector(v)

! dummy variables
TYPE(real8_vector) :: v

INTEGER(KIND=int4) :: error

DEALLOCATE(v%r,STAT=error)

END SUBROUTINE dealloc_real8_vector 

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE assign_real8_vector(v,r)

! dummy variables
TYPE(real8_vector), INTENT(OUT) :: v
REAL(KIND=real8),   INTENT(IN)  :: r

v%r = r

END SUBROUTINE assign_real8_vector 

!!---------------------------------------------------------------------------------------------------------------------


END MODULE real8_vector_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE real8_vector_array_class 

USE basic_types
USE io_files

IMPLICIT NONE

TYPE real8_vector_array
  REAL(KIND=real8), DIMENSION(:), POINTER :: r
END TYPE real8_vector_array


INTERFACE allocate
  MODULE PROCEDURE alloc_real8_vector_array
END INTERFACE

INTERFACE deallocate
  MODULE PROCEDURE dealloc_real8_vector_array
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE assign_real8_vector_array
END INTERFACE 

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE alloc_real8_vector_array(v,lb,ub)

! dummy variables
TYPE(real8_vector_array), DIMENSION(:), POINTER     :: v
INTEGER(KIND=int4),       DIMENSION(:), INTENT(IN)  :: lb, ub

! local variables
INTEGER(KIND=int4) :: n, i, error, l, u

n = SIZE(lb)

ALLOCATE(v(1:n),STAT=error)
IF(error /= 0) THEN
  WRITE(error_unit,*) "allocation failure of vector"
  STOP
END IF

DO i=1,n
  l = lb(i) ; u = ub(i)
  ALLOCATE(v(i)%r(l:u),STAT=error)
  IF(error /= 0) THEN
    WRITE(error_unit,*) "allocation failure of vector"
    STOP
  END IF
END DO


END SUBROUTINE alloc_real8_vector_array

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE dealloc_real8_vector_array(v)

! dummy variables
TYPE(real8_vector_array), DIMENSION(:), POINTER :: v

INTEGER(KIND=int4) :: i, error

DO i=SIZE(v),1,-1
  DEALLOCATE(v(i)%r,STAT=error)
END DO
DEALLOCATE(v,STAT=error)

END SUBROUTINE dealloc_real8_vector_array

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE assign_real8_vector_array(v,r)

! dummy variables
TYPE(real8_vector_array), DIMENSION(:), INTENT(OUT) :: v
REAL(KIND=real8),                       INTENT(IN)  :: r

INTEGER(KIND=int4) :: i

DO i=1,SIZE(v)
  v(i)%r = r
END DO

END SUBROUTINE assign_real8_vector_array 

!!---------------------------------------------------------------------------------------------------------------------


END MODULE real8_vector_array_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE real8_matrix_array_class 

USE basic_types
USE io_files

IMPLICIT NONE

TYPE real8_matrix_array
  REAL(KIND=real8), DIMENSION(:,:), POINTER :: r
END TYPE real8_matrix_array

INTERFACE deallocate
  MODULE PROCEDURE dealloc_real8_matrix_array
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE assign_real8_matrix_array
END INTERFACE 

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE dealloc_real8_matrix_array(m)

! dummy variables
TYPE(real8_matrix_array), DIMENSION(:,:), POINTER :: m

INTEGER(KIND=int4) :: i, j, error

DO i=SIZE(m,DIM=1),1,-1
  DO j=SIZE(m,DIM=2),1,-1
    IF(ASSOCIATED(m(i,j)%r)) DEALLOCATE(m(i,j)%r,STAT=error)
  END DO
END DO
DEALLOCATE(m,STAT=error)

END SUBROUTINE dealloc_real8_matrix_array

!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE assign_real8_matrix_array(m,r)

! dummy variables
TYPE(real8_matrix_array), DIMENSION(:,:), INTENT(OUT) :: m
REAL(KIND=real8),                         INTENT(IN)  :: r

INTEGER(KIND=int4) :: i, j

DO i=1,SIZE(m,DIM=1)
  DO j=1,SIZE(m,DIM=2)
    IF(ASSOCIATED(m(i,j)%r)) m(i,j)%r = r
  END DO
END DO

END SUBROUTINE assign_real8_matrix_array 

!!---------------------------------------------------------------------------------------------------------------------


END MODULE real8_matrix_array_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE real4_vector_class 

USE basic_types
USE io_files

IMPLICIT NONE

TYPE real4_vector
  REAL(KIND=real4), DIMENSION(:), POINTER :: r
END TYPE real4_vector

INTERFACE norm
  MODULE PROCEDURE real4_norm
END INTERFACE

INTERFACE allocate
  MODULE PROCEDURE alloc_real4_vector
END INTERFACE

INTERFACE deallocate
  MODULE PROCEDURE dealloc_real4_vector
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE assign_real4_vector
END INTERFACE 

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------

FUNCTION real4_norm(v)

TYPE(real4_vector) :: v
REAL(KIND=real8)   ::real4_norm

real4_norm=DOT_PRODUCT(v%r,v%r)
real4_norm=SQRT(real4_norm)

END FUNCTION real4_norm
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE alloc_real4_vector(v,n)

! dummy variables
TYPE(real4_vector), INTENT(INOUT) :: v
INTEGER(KIND=int4), INTENT(IN)   :: n


INTEGER(KIND=int4)               :: error


ALLOCATE(v%r(n),STAT=error)
IF(error /= 0) THEN
  WRITE(error_unit,*) "allocation failure of vector"
  STOP
END IF


END SUBROUTINE alloc_real4_vector
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE dealloc_real4_vector(v)

! dummy variables
TYPE(real4_vector) :: v

INTEGER(KIND=int4) :: error

IF(ASSOCIATED(v%r))DEALLOCATE(v%r,STAT=error)

END SUBROUTINE dealloc_real4_vector 
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE assign_real4_vector(v,r)

! dummy variables
TYPE(real4_vector), INTENT(OUT) :: v
REAL(KIND=real4),  INTENT(IN)  :: r

v%r=r

END SUBROUTINE assign_real4_vector 
!!-------------------------------------------------------------------


END MODULE real4_vector_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE skyline_matrix_class

USE real8_vector_class
USE real8_vector_array_class

IMPLICIT NONE

TYPE skyline_matrix
  CHARACTER(LEN=20)                               :: kind
  LOGICAL :: singular, factorised, positive
  REAL(KIND=real8),         DIMENSION(:), POINTER :: d
  TYPE(real8_vector_array), DIMENSION(:), POINTER :: u
  TYPE(real8_vector_array), DIMENSION(:), POINTER :: l
END TYPE skyline_matrix


INTERFACE allocate
  MODULE PROCEDURE allocate_skyline
END INTERFACE

INTERFACE deallocate
  MODULE PROCEDURE deallocate_skyline
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE assign_skyline
END INTERFACE 

INTERFACE mult
  MODULE PROCEDURE mult_skyline
END INTERFACE

INTERFACE Crout_factors
  MODULE PROCEDURE Crout_factors_skyline
END INTERFACE

INTERFACE solve
  MODULE PROCEDURE solve_skyline
END INTERFACE 

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE allocate_skyline(m,lb,ub)

! dummy arguments
TYPE(skyline_matrix),          INTENT(INOUT) :: m
INTEGER(KIND=int4), DIMENSION(:), INTENT(IN) :: lb, ub


matrix_kind: SELECT CASE (m%kind)

  CASE ("symmetric")
    CALL alloc_sym_skyline(m,lb,ub)
  CASE ("unsymmetric")
    CALL alloc_unsym_skyline(m,lb,ub)
  CASE DEFAULT
    WRITE(error_unit,*) "matrix solver type not recognized"
    STOP

END SELECT matrix_kind


END SUBROUTINE allocate_skyline 
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE deallocate_skyline(m)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m


matrix_kind: SELECT CASE (m%kind)

  CASE ("symmetric")
    CALL dealloc_sym_skyline(m)
  CASE ("unsymmetric")
    CALL dealloc_unsym_skyline(m)
  CASE DEFAULT
    WRITE(error_unit,*) "matrix solver type not recognized"
    STOP

END SELECT matrix_kind


END SUBROUTINE deallocate_skyline 
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE assign_skyline(m,r)

! dummy arguments
TYPE(skyline_matrix), INTENT(OUT) :: m
REAL(KIND=real8),     INTENT(IN)  :: r


matrix_kind: SELECT CASE (m%kind)

  CASE ("symmetric")
    CALL assign_sym_skyline(m,r)
  CASE ("unsymmetric")
    CALL assign_unsym_skyline(m,r)
  CASE DEFAULT
    WRITE(error_unit,*) "matrix type not recognized for assignment"
    STOP

END SELECT matrix_kind


END SUBROUTINE assign_skyline 
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE mult_skyline(m,v,w)

! dummy arguments
TYPE(skyline_matrix), INTENT(IN)  :: m
TYPE(real8_vector),   INTENT(IN)  :: v
TYPE(real8_vector),   INTENT(OUT) :: w


matrix_kind: SELECT CASE (m%kind)

  CASE ("symmetric")
    CALL sym_skyline_mult_vector(m,v,w)
  CASE ("unsymmetric")
    CALL unsym_skyline_mult_vector(m,v,w)
  CASE DEFAULT
    WRITE(error_unit,*) "matrix type not recognized for multiplication"
    STOP

END SELECT matrix_kind

END SUBROUTINE mult_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE Crout_factors_skyline(m)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m


matrix_kind: SELECT CASE (m%kind)

  CASE ("symmetric")
    CALL Crout_factors_sym_skyline(m)
  CASE ("unsymmetric")
    CALL Crout_factors_unsym_skyline(m)
  CASE DEFAULT
    WRITE(error_unit,*) "matrix type not recognized for factorization"
    STOP

END SELECT matrix_kind


END SUBROUTINE Crout_factors_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE solve_skyline(m,v)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m
TYPE(real8_vector),   INTENT(INOUT) :: v


matrix_kind: SELECT CASE (m%kind)

  CASE ("symmetric")
    CALL solve_sym_skyline(m,v)
  CASE ("unsymmetric")
    CALL solve_unsym_skyline(m,v)
  CASE DEFAULT
    WRITE(error_unit,*) "matrix type not recognized for Crout solution"
    STOP

END SELECT matrix_kind

END SUBROUTINE solve_skyline
!!-------------------------------------------------------------------


SUBROUTINE alloc_sym_skyline(m,lb,ub)

! dummy arguments
TYPE(skyline_matrix),             INTENT(INOUT) :: m
INTEGER(KIND=int4), DIMENSION(:), INTENT(IN)    :: lb, ub

! local variables
INTEGER(KIND=int4) :: n, error

n = SIZE(lb)

ALLOCATE(m%d(1:n),STAT=error)
IF (error /= 0) THEN
  WRITE(error_unit,*) "allocation failure of symmetric skyline"
  STOP 
END IF

CALL allocate(m%u,lb,ub)

END SUBROUTINE alloc_sym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE alloc_unsym_skyline(m,lb,ub)

! dummy arguments
TYPE(skyline_matrix),             INTENT(INOUT) :: m
INTEGER(KIND=int4), DIMENSION(:), INTENT(IN)    :: lb, ub    ! caution: ub is present here only to compose the third argument

! local variables
INTEGER(KIND=int4) :: n,j,l,u,error

n = SIZE(lb,dim=1)


ALLOCATE(m%d(1:n),STAT=error)
IF (error /= 0) THEN
  WRITE(error_unit,*) "allocation failure of unsymmetric skyline"
  STOP
END IF

ALLOCATE(m%u(n),m%l(n),STAT=error)
IF (error /= 0) THEN
  WRITE(error_unit,*) "allocation failure of unsymmetric skyline"
  STOP
END IF

DO j=1,n
  l=lb(j); u=j-1
  ALLOCATE(m%u(j)%r(l:u),m%l(j)%r(l:u),STAT=error)
  IF (error /= 0) THEN
    WRITE(error_unit,*) "allocation failure of unsymmetric skyline"
    STOP
  END IF
END DO

END SUBROUTINE alloc_unsym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE dealloc_sym_skyline(m)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m

! local variables
INTEGER(KIND=int4) :: j, error


DO j=SIZE(m%d),1,-1
  DEALLOCATE(m%u(j)%r,STAT=error)
END DO

DEALLOCATE(m%u,STAT=error)
DEALLOCATE(m%d,STAT=error)

END SUBROUTINE dealloc_sym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE dealloc_unsym_skyline(m)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m

! local variables
INTEGER(int4) :: n,j,error


IF (ASSOCIATED(m%d).OR.ASSOCIATED(m%u).OR.ASSOCIATED(m%l)) THEN

 n = SIZE(m%d)

 DO j=n,1,-1
   DEALLOCATE(m%l(j)%r,m%u(j)%r,STAT=error)
   IF(error/=0)THEN
     WRITE(error_unit,*) "deallocation failure of unsymmetric skyline"
   END IF
 END DO

 DEALLOCATE(m%l,m%u,STAT=error)
 IF(error/=0)THEN
   WRITE(error_unit,*) "deallocation failure of unsymmetric skyline"
 END IF

 DEALLOCATE(m%d,STAT=error)
 IF(error/=0)THEN
   WRITE(error_unit,*) "deallocation failure of unsymmetric skyline"
 END IF

END IF

END SUBROUTINE dealloc_unsym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE assign_sym_skyline(m,r)

! dummy variables
TYPE(skyline_matrix), INTENT(OUT) :: m
REAL(KIND=real8),     INTENT(IN)  :: r

INTEGER(KIND=int4) :: j


m%factorised=.FALSE.

DO j=1,SIZE(m%d)
  m%u(j)%r=r
END DO

m%d=r

END SUBROUTINE assign_sym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE assign_unsym_skyline(m,r)

! dummy arguments
TYPE(skyline_matrix), INTENT(OUT) :: m
REAL(KIND=real8),     INTENT(IN)  :: r

! local variables
INTEGER(int4) :: i


m%factorised=.FALSE.

m%d=r
DO i=1,SIZE(m%d)
  m%u(i)%r=r
  m%l(i)%r=r
END DO

END SUBROUTINE assign_unsym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sym_skyline_mult_vector(m,v,w)

! dummy arguments
TYPE(skyline_matrix), INTENT(IN)  :: m
TYPE(real8_vector),   INTENT(IN)  :: v
TYPE(real8_vector),   INTENT(OUT) :: w

! local variables
INTEGER(KIND=int4) :: n,l,u,i

! array dimension
n = SIZE(m%d)

! clear result
w%r = 0.0_real8

DO i=1,n

  ! bounds
  l = LBOUND(m%u(i)%r,1) ; u = UBOUND(m%u(i)%r,1)

  ! lower triangle
  w%r(i)   = w%r(i)+DOT_PRODUCT(m%u(i)%r(l:u),v%r(l:u))

  ! diagonal
  w%r(i)   = w%r(i)+m%d(i)*v%r(i)

  ! upper triangle
  w%r(l:u) = w%r(l:u)+m%u(i)%r(l:u)*v%r(i)

END DO

END SUBROUTINE sym_skyline_mult_vector
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE unsym_skyline_mult_vector(m,v,w)

! dummy arguments
TYPE(skyline_matrix), INTENT(IN)  :: m
TYPE(real8_vector),   INTENT(IN)  :: v
TYPE(real8_vector),   INTENT(OUT) :: w

! local variables
INTEGER(int4) :: i,l,u


! diagonal
w%r = m%d*v%r

DO i=2,SIZE(m%d)

  ! column bounds
  l = LBOUND(m%u(i)%r,1) ; u = UBOUND(m%u(i)%r,1)

  ! lower triangle
  w%r(i) = w%r(i)+DOT_PRODUCT(m%l(i)%r(l:u),v%r(l:u))

  ! upper triangle
  w%r(l:u) = w%r(l:u)+m%u(i)%r(l:u)*v%r(i)

END DO

END SUBROUTINE unsym_skyline_mult_vector
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE Crout_reduction_sym_skyline (m,v)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m
TYPE(real8_vector),   INTENT(INOUT) :: v

! local variables
INTEGER(KIND=int4) :: n,i,li,ui

! system dimension
n = SIZE(m%d)

! forward reduction
DO i=2,n
  li = LBOUND(m%u(i)%r,1) ; ui = UBOUND(m%u(i)%r,1)
  v%r(i) = v%r(i) - DOT_PRODUCT(m%u(i)%r(li:ui),v%r(li:ui))
END DO

! diagonal reduction
v%r = v%r / m%d

! backward reduction
DO i=n,2,-1
  li = LBOUND(m%u(i)%r,1) ; ui = UBOUND(m%u(i)%r,1)
  v%r(li:ui) = v%r(li:ui) - m%u(i)%r(li:ui)*v%r(i)
END DO

END SUBROUTINE Crout_reduction_sym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE Crout_reduction_unsym_skyline (a,x)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: a
TYPE(real8_vector),   INTENT(INOUT) :: x

! local variables
INTEGER(KIND=int4) :: n,i,l,u


! system dimension
n = SIZE(a%d)

! forward reduction
DO i=2,n
  l = LBOUND(a%l(i)%r,1) ; u = UBOUND(a%l(i)%r,1)
  x%r(i) = x%r(i) - DOT_PRODUCT(a%l(i)%r(l:u),x%r(l:u))
END DO

! diagonal reduction
x%r = x%r / a%d

! backward reduction
DO i=n,2,-1
  l = LBOUND(a%u(i)%r,1) ; u = UBOUND(a%u(i)%r,1)
  x%r(l:u) = x%r(l:u) - a%u(i)%r(l:u)*x%r(i)
END DO

END SUBROUTINE Crout_reduction_unsym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE Crout_factors_sym_skyline(m)   

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m

! local variables
REAL(KIND=real8), PARAMETER :: tol = TINY(1.0_real4) 
INTEGER(KIND=int4) :: n,i,j,li,lj,ui,uj,lx
REAL(KIND=real8)   :: tmp


! system dimension
n = SIZE(m%d)

! initialise logical flags
m%singular = .FALSE.
m%positive = .TRUE.

! loop over columns
DO j=1,n

IF ( SIZE(m%u(j)%r)>0 ) THEN

  lj = LBOUND(m%u(j)%r,1) ; uj = j-1

  DO i=lj+1,uj
    IF ( SIZE(m%u(i)%r)>0 ) THEN
	  li = LBOUND(m%u(i)%r,1) ; ui = i-1
      lx = MAX(li,lj)
      m%u(j)%r(i) = m%u(j)%r(i) - DOT_PRODUCT(m%u(i)%r(lx:ui),m%u(j)%r(lx:ui))
    END IF
  END DO

  DO i=lj,uj
    tmp = m%u(j)%r(i)
    m%u(j)%r(i) = tmp / m%d(i)
    m%d(j) = m%d(j) - tmp*m%u(j)%r(i)
  END DO

END IF

  ! test for positiveness
  IF ( m%d(j) < -tol ) m%positive = .FALSE.

  ! test for singularity
  IF ( ABS(m%d(j)) < tol ) THEN
    m%singular = .TRUE.
    EXIT
  END IF

END DO 

m%factorised=.TRUE.

END SUBROUTINE Crout_factors_sym_skyline 
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE Crout_factors_unsym_skyline(m)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m

! local variables
REAL(KIND=real8), PARAMETER :: tol = TINY(1.0_real4) 
INTEGER(KIND=int4) :: n,i,j,li,lj,ui,uj,lx
REAL(KIND=real8)   :: tmp


! system dimension
n = SIZE(m%d)

! initialise logical flags
m%singular = .FALSE.
m%positive = .TRUE.

! loop over columns
DO j=1,n

  lj = LBOUND(m%u(j)%r,1) ; uj = UBOUND(m%u(j)%r,1)

  DO i=lj+1,uj
    li = LBOUND(m%u(i)%r,1) ; ui = UBOUND(m%u(i)%r,1)
    lx = MAX(li,lj)
    m%u(j)%r(i) = m%u(j)%r(i) - DOT_PRODUCT(m%l(i)%r(lx:ui),m%u(j)%r(lx:ui))
    m%l(j)%r(i) = m%l(j)%r(i) - DOT_PRODUCT(m%u(i)%r(lx:ui),m%l(j)%r(lx:ui))
  END DO

  DO i=lj,uj
    tmp = m%u(j)%r(i)
    m%u(j)%r(i) = m%u(j)%r(i) / m%d(i)
    m%l(j)%r(i) = m%l(j)%r(i) / m%d(i)
    m%d(j) = m%d(j) - tmp*m%l(j)%r(i)
  END DO

  ! test for positiveness
  IF ( m%d(j) < -tol ) m%positive = .FALSE.

  ! test for singularity
  IF ( ABS(m%d(j)) < tol ) THEN
    m%singular = .TRUE.
    EXIT
  END IF 

END DO 

m%factorised=.TRUE.

END SUBROUTINE Crout_factors_unsym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE solve_sym_skyline(m,v)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: m
TYPE(real8_vector),   INTENT(INOUT) :: v


! Crout factorisation
IF (.NOT.m%factorised) CALL Crout_factors_sym_skyline(m)

! Crout reduction
IF (.NOT.m%singular) CALL Crout_reduction_sym_skyline(m,v)

END SUBROUTINE solve_sym_skyline
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE solve_unsym_skyline(a,x)

! dummy arguments
TYPE(skyline_matrix), INTENT(INOUT) :: a
TYPE(real8_vector),   INTENT(INOUT) :: x

! Crout factorisation
IF (.NOT.a%factorised) CALL Crout_factors_unsym_skyline(a)

! Crout reduction
IF (.NOT.a%singular) CALL Crout_reduction_unsym_skyline(a,x)

END SUBROUTINE solve_unsym_skyline
!!---------------------------------------------------------------------------------------------------------------------


END MODULE skyline_matrix_class 

!!=============================================================================================================================
!!=============================================================================================================================


MODULE pefmat

USE tensor_class
USE skyline_matrix_class
USE spectral_decomposition_class
USE jacobi_class
USE rotation_class
USE real8_vector_class
USE real8_vector_array_class
USE real8_matrix_array_class
USE real4_vector_class


IMPLICIT NONE


INTERFACE OPERATOR(.matrix.)
  MODULE PROCEDURE matrix_eq_diagonal
END INTERFACE 

INTERFACE OPERATOR(.add.)
  MODULE PROCEDURE m_add_scalar
  MODULE PROCEDURE scalar_add_m
  MODULE PROCEDURE m_add_v
  MODULE PROCEDURE v_add_m
END INTERFACE

INTERFACE OPERATOR(.xt.)
  MODULE PROCEDURE v_xt_v
  MODULE PROCEDURE m_xt_m
END INTERFACE

INTERFACE OPERATOR(.tx.)
  MODULE PROCEDURE m_tx_v
  MODULE PROCEDURE m_tx_m
  MODULE PROCEDURE v_tx_v
END INTERFACE

INTERFACE OPERATOR(.x.)
  MODULE PROCEDURE m_x_v
  MODULE PROCEDURE m_x_m
END INTERFACE

INTERFACE OPERATOR(.ct.)
  MODULE PROCEDURE m_ct_m
END INTERFACE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION matrix_eq_diagonal(d)

! dummy arguments
REAL(KIND=real8), DIMENSION(:),   INTENT(IN) :: d

REAL(KIND=real8), DIMENSION(SIZE(d),SIZE(d)) :: matrix_eq_diagonal

INTEGER(KIND=int4) :: i

matrix_eq_diagonal = 0.0_real8

DO i=1,SIZE(d)
  matrix_eq_diagonal(i,i)=d(i)
END DO

END FUNCTION matrix_eq_diagonal
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_add_scalar (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:), INTENT(IN) :: a
REAL(KIND=real8),                 INTENT(IN) :: b

REAL(KIND=real8), DIMENSION(SIZE(a,DIM=1),SIZE(a,DIM=2)) :: m_add_scalar

! local variables
INTEGER(KIND=int4) :: i

m_add_scalar = a

DO i=1,MIN(SIZE(a,DIM=1),SIZE(a,DIM=2))
  m_add_scalar(i,i) = m_add_scalar(i,i) + b
END DO

END FUNCTION m_add_scalar
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION scalar_add_m (a,b)

! dummy arguments
REAL(KIND=real8),                 INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(:,:), INTENT(IN) :: b

REAL(KIND=real8), DIMENSION(SIZE(b,DIM=1),SIZE(b,DIM=2)) :: scalar_add_m

! local variables
INTEGER(KIND=int4) :: i

scalar_add_m = b

DO i=1,MIN(SIZE(b,DIM=1),SIZE(b,DIM=2))
  scalar_add_m(i,i) = scalar_add_m(i,i) + a
END DO

END FUNCTION scalar_add_m
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_add_v (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:), INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(:),   INTENT(IN) :: b

REAL(KIND=real8), DIMENSION(SIZE(b),SIZE(b)) :: m_add_v

! local variables
INTEGER(KIND=int4) :: i

m_add_v = a

DO i=1,SIZE(b)
  m_add_v(i,i) = m_add_v(i,i) + b(i)
END DO

END FUNCTION m_add_v
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION v_add_m (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:),   INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(:,:), INTENT(IN) :: b

REAL(KIND=real8), DIMENSION(SIZE(a),SIZE(a)) :: v_add_m

! local variables
INTEGER(KIND=int4) :: i

v_add_m = b

DO i=1,SIZE(a)
  v_add_m(i,i) = v_add_m(i,i) + a(i)
END DO

END FUNCTION v_add_m
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION v_xt_v (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:),  INTENT(IN) :: a, b

REAL(KIND=real8), DIMENSION(SIZE(a),SIZE(b)) :: v_xt_v

! local variables
INTEGER(KIND=int4) :: i, j

DO i=1,SIZE(a)
  DO j=1,SIZE(b)
    v_xt_v(i,j) = a(i)*b(j)
  END DO
END DO

END FUNCTION v_xt_v
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_xt_m (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:),  INTENT(IN) :: a, b

REAL(KIND=real8), DIMENSION(SIZE(a,DIM=1),SIZE(b,DIM=1)) :: m_xt_m

m_xt_m = MATMUL(a,TRANSPOSE(b))

END FUNCTION m_xt_m
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_tx_m (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:),  INTENT(IN) :: a, b

REAL(KIND=real8), DIMENSION(SIZE(a,DIM=2),SIZE(b,DIM=2)) :: m_tx_m

m_tx_m = MATMUL(TRANSPOSE(a),b)

END FUNCTION m_tx_m
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_tx_v (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:),  INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(:),    INTENT(IN) :: b

REAL(KIND=real8), DIMENSION(SIZE(a,DIM=2)) :: m_tx_v

m_tx_v = MATMUL(TRANSPOSE(a),b)

END FUNCTION m_tx_v
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION v_tx_v (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:),  INTENT(IN) :: a, b

REAL(KIND=real8) :: v_tx_v

v_tx_v = DOT_PRODUCT(a,b)

END FUNCTION v_tx_v
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_x_m (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:), INTENT(IN) :: a, b

REAL(KIND=real8), DIMENSION(SIZE(a,DIM=1),SIZE(b,DIM=2)) :: m_x_m

m_x_m = MATMUL(a,b)

END FUNCTION m_x_m
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_x_v (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:),  INTENT(IN) :: a
REAL(KIND=real8), DIMENSION(:),    INTENT(IN) :: b

REAL(KIND=real8), DIMENSION(SIZE(b)) :: m_x_v

m_x_v = MATMUL(a,b)

END FUNCTION m_x_v
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION m_ct_m (a,b)

! dummy arguments
REAL(KIND=real8), DIMENSION(:,:), INTENT(IN) :: a, b

REAL(KIND=real8), DIMENSION(SIZE(a,DIM=2),SIZE(a,DIM=2)) :: m_ct_m

m_ct_m = a .tx. (b .x. a)

END FUNCTION m_ct_m
!!---------------------------------------------------------------------------------------------------------------------


END MODULE pefmat 

!!=============================================================================================================================
!!=============================================================================================================================


MODULE array_of_vectors_class

USE basic_types
USE io_files

IMPLICIT NONE

TYPE array_of_vectors
  REAL(KIND=real8), DIMENSION(3) :: r
END TYPE array_of_vectors


INTERFACE norm
  MODULE PROCEDURE array_of_vectors_norm
END INTERFACE

!INTERFACE ASSIGNMENT(=)
!  MODULE PROCEDURE array_of_vectors_equality
!END INTERFACE 

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


FUNCTION array_of_vectors_norm(v)

TYPE(array_of_vectors), DIMENSION(:), POINTER :: v
REAL(KIND=real8) :: array_of_vectors_norm
INTEGER(KIND=int4) :: i


! compute norm
array_of_vectors_norm = 0.0_real8
DO i=1,SIZE(v)
  array_of_vectors_norm = array_of_vectors_norm + DOT_PRODUCT(v(i)%r,v(i)%r)
END DO
array_of_vectors_norm = SQRT(array_of_vectors_norm)


END FUNCTION array_of_vectors_norm
!!---------------------------------------------------------------------------------------------------------------------


!SUBROUTINE array_of_vectors_equality(v,x)

! dummy variables
!TYPE(array_of_vectors), DIMENSION(:), POINTER, INTENT(OUT) :: v
!TYPE(array_of_vectors), DIMENSION(:), POINTER, INTENT(IN)  :: x
!INTEGER(KIND=int4) :: i

! perform equality
!DO i=1,SIZE(x)
!  v(i)%r = x(i)%r
!END DO

!END SUBROUTINE array_of_vectors_equality 
!!-------------------------------------------------------------------


END MODULE array_of_vectors_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE contact_history_array_class

USE basic_types

IMPLICIT NONE

TYPE contact_pair_data
  INTEGER(KIND=int4) :: number
  REAL(KIND=real8)               :: normalspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: normalspring_force
  REAL(KIND=real8), DIMENSION(3) :: fricspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: fricspring_force
  REAL(KIND=real8), DIMENSION(3) :: rolresspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: contact_force
  REAL(KIND=real8), DIMENSION(3) :: friction_force  
END TYPE contact_pair_data

TYPE contact_history_array
  TYPE(contact_pair_data), DIMENSION(:), ALLOCATABLE :: pair
END TYPE contact_history_array

END MODULE contact_history_array_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE particle_data_types

USE pefmat

IMPLICIT NONE

TYPE particle_parameters
  CHARACTER(LEN=50) :: kind  ! sphere1, jet_sphere, granular_sphere, etc (one of the various particles of the program)
  CHARACTER(LEN=50) :: shape ! "sphere", "ellipsoid", "superelipsoid", etc (may be important for contact detection scheme)
END TYPE particle_parameters

TYPE particle_control_variables
  LOGICAL :: velocity, acceleration
  LOGICAL :: momenta, energy
  CHARACTER(LEN=20) :: shape, contact_type
END TYPE particle_control_variables

TYPE contact_particle_pair_data
  INTEGER(KIND=int4) :: part_number
  REAL(KIND=real8)               :: normalspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: normalspring_force
  REAL(KIND=real8), DIMENSION(3) :: fricspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: fricspring_force
  REAL(KIND=real8), DIMENSION(3) :: rolresspring_elongation
END TYPE contact_particle_pair_data 

TYPE contact_wall_pair_data
  INTEGER(KIND=int4) :: wall_number
  REAL(KIND=real8)              :: normalspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: normalspring_force
  REAL(KIND=real8), DIMENSION(3) :: fricspring_elongation
  REAL(KIND=real8), DIMENSION(3) :: fricspring_force
  REAL(KIND=real8), DIMENSION(3) :: rolresspring_elongation
END TYPE contact_wall_pair_data

TYPE particle_data
  TYPE(particle_parameters)        :: parameters
  TYPE(particle_control_variables) :: cv
  REAL(KIND=real8)                 :: mass, radius, inertia, charge, specific_heat, absorptance
  INTEGER(KIND=int4)               :: material_set_number, pressure_surface_number
  INTEGER(KIND=int4)               :: no_particles_into_verlet_list
  INTEGER(kind=INT4)               :: no_contacting_particles, no_contacting_walls
  REAL(KIND=real8)                 :: initial_temperature, temperature, temperature_dof_codes
  REAL(KIND=real8), DIMENSION(3)   :: xyz_coordinates, initial_velocity
  REAL(KIND=real8), DIMENSION(3)   :: initial_angles, initial_spin
  REAL(KIND=real8), DIMENSION(3)   :: position, velocity, acceleration
  REAL(KIND=real8), DIMENSION(3)   :: angle, spin, spin_acceleration
  REAL(KIND=real8), DIMENSION(3)   :: nearfield_force, environment_force, initial_given_force, initial_given_moment 
  REAL(KIND=real8), DIMENSION(3)   :: contact_force, friction_force, friction_moment, rolling_resistance_moment, epsilon_force
  REAL(KIND=real8), DIMENSION(3)   :: adhesion_force, adhesion_moment
  REAL(KIND=real8), DIMENSION(3)   :: linear_mom, angular_mom
  REAL(KIND=real8)                 :: conduction_heat_power
  REAL(KIND=real8), DIMENSION(6)   :: harmonic_constraints_data
  INTEGER(KIND=int4), DIMENSION(6) :: translational_dof_codes, rotational_dof_codes
  INTEGER(KIND=int4), DIMENSION(3) :: cell_address
  INTEGER(KIND=int4), DIMENSION(:), ALLOCATABLE :: verlet_list
  TYPE(contact_particle_pair_data), DIMENSION(:), ALLOCATABLE :: contacting_particles
  TYPE(contact_wall_pair_data),     DIMENSION(:), ALLOCATABLE :: contacting_walls
  REAL(KIND=real8)  :: collision_duration
  LOGICAL           :: collision
  CHARACTER(LEN=11) :: friction_status
  REAL(KIND=real8)  :: total_energy, kinetic_energy, internal_energy
END TYPE particle_data


END MODULE particle_data_types

!!=============================================================================================================================
!!=============================================================================================================================


MODULE system_data_types

USE pefmat
USE particle_data_types

IMPLICIT NONE

TYPE material_set_data
  REAL(KIND=real8) :: mass_density, charge_density
  REAL(KIND=real8) :: elasticity_modulus, poisson_coeff
  REAL(KIND=real8) :: coefficient_of_restitution
  REAL(KIND=real8) :: contact_damping_ratio, friction_damping_ratio
  REAL(KIND=real8) :: static_friction_coeff, dynamic_friction_coeff
  REAL(KIND=real8) :: rolling_resistance_coeff, rolling_resistance_damping_ratio
  REAL(KIND=real8) :: specific_heat, thermal_conductivity, degrading_temperature
  REAL(KIND=real8) :: drag_heating_efficiency, radiative_efficiency, absorptance
END TYPE material_set_data

TYPE wall_contact_list_data
  INTEGER(KIND=int4) :: part_number
  REAL(KIND=real8), DIMENSION(3) :: contact_force
  REAL(KIND=real8), DIMENSION(3) :: friction_force 
END TYPE wall_contact_list_data

TYPE wall_data
  CHARACTER(LEN=50)  :: kind, thermal_kind
  REAL(KIND=real8)   :: charge, coefficient_of_restitution, radius
  REAL(KIND=real8)   :: contact_damping_ratio, friction_damping_ratio
  REAL(KIND=real8)   :: static_friction_coeff, dynamic_friction_coeff
  REAL(KIND=real8)   :: rolling_resistance_coeff, rolling_resistance_damping_ratio
  INTEGER(kind=INT4) :: no_contacting_particles
  REAL(KIND=real8)   :: maximum_force_due_to_contacts  
  REAL(KIND=real8), DIMENSION(3) :: initial_point_position, initial_outside_normal, initial_velocity
  REAL(KIND=real8), DIMENSION(3) :: harmonic_amplitudes, harmonic_frequencies
  REAL(KIND=real8), DIMENSION(3) :: initial_acceleration, initial_driving_force
  REAL(KIND=real8), DIMENSION(3) :: point_position, velocity, acceleration, outside_normal, contact_force, friction_force
  REAL(KIND=real8)               :: initial_temperature, temperature, initial_heating_rate, heating_rate
  TYPE(wall_contact_list_data), DIMENSION(:), ALLOCATABLE :: contacting_particles
END TYPE wall_data

TYPE external_fields_data
  REAL(KIND=real8), DIMENSION(3) :: gravity_accel_vector
  REAL(KIND=real8), DIMENSION(3) :: electric_field_vector
  REAL(KIND=real8), DIMENSION(3) :: magnetic_field_vector
  REAL(KIND=real8) :: emfields_xi, emfields_xf, emfields_yi, emfields_yf, emfields_zi, emfields_zf
END TYPE external_fields_data

TYPE environment_fluid_data
  REAL(KIND=real8) :: density, viscosity
  REAL(KIND=real8) :: temperature, thermal_conductivity, specific_heat
  REAL(KIND=real8), DIMENSION(3) :: velocity
END TYPE environment_fluid_data

TYPE nearfields_data
  REAL(KIND=real8) :: kij_attr, kij_rep
  REAL(KIND=real8) :: expij_attr, expij_rep
  REAL(KIND=real8) :: cutoff_distance
END TYPE nearfields_data

TYPE adhesion_model_data
  CHARACTER(LEN=50)  :: model_type
  INTEGER(KIND=int4) :: no_parameters
  REAL(KIND=real8), DIMENSION(:), ALLOCATABLE :: model_parameters
  !REAL(KIND=real8)   :: kij, expon, epsilon_crit
END TYPE adhesion_model_data

TYPE spring_properties_set_data
  REAL(KIND=real8) :: stiffness
  REAL(KIND=real8) :: elasticity_modulus
  REAL(KIND=real8) :: critical_stress
  REAL(KIND=real8) :: dashpot_constant
END TYPE spring_properties_set_data

TYPE spring_data
  INTEGER(KIND=int4), DIMENSION(2) :: connectivities
  INTEGER(KIND=int4) :: properties_set_number
  REAL(KIND=real8)   :: initial_length
  INTEGER(KIND=int4) :: connecting_wall_number, pressure_surface_number
END TYPE spring_data

TYPE pressure_surface_data
  REAL(KIND=real8) :: initial_pressure
  REAL(KIND=real8) :: initial_influence_area
END TYPE pressure_surface_data

TYPE external_heating_device_data
  CHARACTER(LEN=50) :: kind
  REAL(KIND=real8)  :: initial_intensity, initial_intensity_increase_rate
  REAL(KIND=real8)  :: cross_sectional_area, attenuation_coeff, max_penetration
  REAL(KIND=real8)  :: intensity, intensity_increase_rate
  REAL(KIND=real8), DIMENSION(3) :: initial_striking_position, initial_velocity, axial_direction
  REAL(KIND=real8), DIMENSION(3) :: striking_position, velocity
END TYPE external_heating_device_data

TYPE step_control_variables
  REAL(KIND=real8)   :: initial_dt, dt_min, dt_max, current_dt
  REAL(KIND=real8)   :: final_time, final_factor
  CHARACTER(LEN=3)   :: adaptive_time_stepping
  INTEGER(KIND=int4) :: current_substep
  REAL(KIND=real8)   :: current_time
  REAL(KIND=real8)   :: collisions_duration_parameter
END TYPE step_control_variables

TYPE solution_control_variables
  CHARACTER(LEN=50)  :: run_name, problem_type           
  CHARACTER(LEN=80)  :: solver_type 
  CHARACTER(LEN=80)  :: contact_model, rolling_resistance_model
  CHARACTER(LEN=80)  :: contact_detection_algorithm
  CHARACTER(LEN=50)  :: particles_shape           
  REAL(KIND=real8)   :: max_no_iterations, desired_no_iterations
  REAL(KIND=real8)   :: positions_tolerance, velocities_tolerance
  REAL(KIND=real8)   :: angles_tolerance, spin_tolerance
  INTEGER(KIND=int4) :: no_steps, current_step
  TYPE(step_control_variables), DIMENSION(:), POINTER :: step  
END TYPE solution_control_variables

TYPE general_control_variables
  CHARACTER(LEN=20)  :: matrix_solver
  CHARACTER(LEN=21)  :: results_file_format
  CHARACTER(LEN=3)   :: rotational_dofs, temperature_dofs, periodic_bc
  REAL(KIND=real8)   :: dt_for_results_printing
  REAL(KIND=real8)   :: perbc_xmin, perbc_xmax, perbc_ymin, perbc_ymax, perbc_zmin, perbc_zmax
  CHARACTER(LEN=3)   :: compute_system_properties, print_rotational_dofs, print_extra_results_file
  CHARACTER(LEN=3)   :: nearfields_switch, adhesion_switch
  REAL(KIND=real8)   :: verlet_distance, dt_for_verlet_list_update
END TYPE general_control_variables

TYPE grid_cell_data
  INTEGER(KIND=int4) :: no_particles_into_cell
  INTEGER(KIND=int4), DIMENSION(:), ALLOCATABLE :: particle_list
END TYPE grid_cell_data

TYPE grid_data
  INTEGER(KIND=int4) :: ndivx, ndivy, ndivz
  REAL(KIND=real8)   :: xbeg, xend, ybeg, yend, zbeg, zend
  REAL(KIND=real8)   :: dt_for_cell_list_update
  TYPE(grid_cell_data), DIMENSION(:,:,:), POINTER :: cell
END TYPE grid_data

TYPE system_data
  CHARACTER(LEN=80)  :: name
  INTEGER(KIND=int4) :: no_particles, no_constrained_particles, no_constrained_particles_with_harmonic_constraints
  INTEGER(KIND=int4) :: no_particles_with_initial_given_forces_and_moments 
  INTEGER(KIND=int4) :: no_rigid_walls, no_springs, no_pressure_surfaces, no_external_heating_devices 
  INTEGER(KIND=int4) :: no_material_sets, no_nearfields_sets, no_spring_properties_sets 
  TYPE(particle_data),                DIMENSION(:), POINTER :: particle 
  TYPE(material_set_data),            DIMENSION(:), POINTER :: material_set 
  TYPE(wall_data),                    DIMENSION(:), POINTER :: wall 
  TYPE(nearfields_data),              DIMENSION(:), POINTER :: nearfields_set
  TYPE(spring_properties_set_data),   DIMENSION(:), POINTER :: spring_properties_set
  TYPE(spring_data),                  DIMENSION(:), POINTER :: spring
  TYPE(pressure_surface_data),        DIMENSION(:), POINTER :: pressure_surface
  TYPE(external_heating_device_data), DIMENSION(:), POINTER :: external_heating_device
  TYPE(external_fields_data)     :: external_fields
  TYPE(environment_fluid_data)   :: environment_fluid
  TYPE(adhesion_model_data)      :: adhesion_model
  REAL(KIND=real8), DIMENSION(3) :: linear_mom, angular_mom
  REAL(KIND=real8), DIMENSION(3) :: center_of_mass
  REAL(KIND=real8)               :: total_energy, kinetic_energy
  REAL(KIND=real8)               :: translational_kinetic_energy, rotational_kinetic_energy
  REAL(KIND=real8)               :: global_damping_coefficient
  TYPE(solution_control_variables) :: solution_cv
  TYPE(general_control_variables)  :: general_cv
  TYPE(grid_data) :: grid
END TYPE system_data


END MODULE system_data_types

!!=============================================================================================================================
!!=============================================================================================================================


MODULE deallocate_arrays

USE system_data_types
USE io_files

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------

SUBROUTINE deallocate_system_arrays(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: error_v, error_cp, error_cw, error_p, error_m, error_w, error_nf, error_ps, &
                      error_s, error_g, error_sps, error_scv, i, ix, iy, iz


! deallocate particles´ verlet lists
IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  DO i=1,system%no_particles
    DEALLOCATE(system%particle(i)%verlet_list,STAT=error_v)
    IF (error_v/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of verlet list of particle ", i
  END DO
END IF

! deallocate particles´ contact lists
DO i=1,system%no_particles
  IF (ALLOCATED(system%particle(i)%contacting_particles)) DEALLOCATE(system%particle(i)%contacting_particles,STAT=error_cp)
  IF (ALLOCATED(system%particle(i)%contacting_walls)) DEALLOCATE(system%particle(i)%contacting_walls,STAT=error_cw)
  IF (error_cp/=0 .OR. error_cw/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of contact lists of particle ", i
END DO

! deallocate particles array
DEALLOCATE(system%particle,STAT=error_p)
IF (error_p/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system particles array"

! deallocate materials array
DEALLOCATE(system%material_set,STAT=error_m)
IF (error_m/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system materials array"

! deallocate walls array
IF (system%no_rigid_walls>0) THEN
  DEALLOCATE(system%wall,STAT=error_w)
  IF (error_w/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system walls array"
END IF

! deallocate nearfields array
IF (system%general_cv%nearfields_switch=="on") THEN
  DEALLOCATE(system%nearfields_set,STAT=error_nf)
  IF (error_nf/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system nearfields array"
END IF

! deallocate pressure surfaces array
IF (system%no_pressure_surfaces>0) THEN
  DEALLOCATE(system%pressure_surface,STAT=error_ps)
  IF (error_ps/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system pressure_surfaces array"
END IF

! deallocate springs array
IF (system%no_springs>0) THEN
  DEALLOCATE(system%spring,STAT=error_s)
  DEALLOCATE(system%spring_properties_set,STAT=error_sps)
  IF (error_s/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system springs array"
  IF (error_sps/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system springs properties sets array"
END IF

! deallocate grid arrays
IF (system%solution_cv%contact_detection_algorithm=="binning" .OR. &
    system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  DO ix=0,system%grid%ndivx+1
    DO iy=0,system%grid%ndivy+1
      DO iz=0,system%grid%ndivz+1
        DEALLOCATE(system%grid%cell(ix,iy,iz)%particle_list)
      END DO
    END DO
  END DO
  DEALLOCATE(system%grid%cell,STAT=error_g)
  IF (error_g/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system grid arrays"
END IF

! deallocate solution cvs steps array
DEALLOCATE(system%solution_cv%step,STAT=error_scv)
IF (error_scv/=0) WRITE(error_unit,*) " CAUTION: deallocation failure of system solution cv steps array"

END SUBROUTINE deallocate_system_arrays
!!---------------------------------------------------------------------------------------------------------------------


END MODULE deallocate_arrays

!!=============================================================================================================================
!!=============================================================================================================================

! Add bellow, as a module (i.e., as a class), each new type of particle; 
! The module should build the force vector for this class of particle.
! (one module for each type of particle)

!!=============================================================================================================================
!!=============================================================================================================================


MODULE sphere1_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere1_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf                        

! local variables
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8) :: kij_attr, expij_attr, kij_rep, expij_rep, norm_rimrj, iradius, jradius, warning_dist
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, rimrj, nij

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri = R(i)%r
vi = V(i)%r

! get near-fields coefficients and exponents
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep

! clear near-field force vector
fnf = 0.0_real8

! compute near-field force vector
DO j=1,system%no_particles
  IF (j/=i) THEN
    jradius = system%particle(j)%radius
    rj = R(j)%r
    rimrj = ri - rj
    norm_rimrj = SQRT(DOT_PRODUCT(rimrj,rimrj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (norm_rimrj<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
      !STOP
    END IF  
    nij = -(one/norm_rimrj)*rimrj
    fnf = fnf + kij_attr*(norm_rimrj**expij_attr)*nij - kij_rep*(norm_rimrj**expij_rep)*nij
  END IF
END DO  

END SUBROUTINE compute_sphere1_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere1_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V  
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv                        

! local variables
INTEGER(KIND=int4) :: i
REAL(KIND=real8) :: mass, charge, ci
REAL(KIND=real8), DIMENSION(3) :: ri, vi, g, E, B


! get particle data
i = particle_number
mass = system%particle(i)%mass
charge = system%particle(i)%charge
ri = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get global damping parameters
ci = system%global_damping_coefficient

! add external fields contribution to environment force vector
fenv = mass*g + charge*E + charge*(vi .vector. B)

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_sphere1_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere1_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_sphere1_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE sphere1_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE sphere2_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere2_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V  
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf, fenv                   

! local variables
INTEGER(KIND=int4) :: i, j
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_rimrj, iradius, jradius, warning_dist
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, rimrj, nij

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri = R(i)%r
vi = V(i)%r

! get near-fields coefficients and exponents
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep

! clear near-fied force vector
fnf = 0.0_real8

! compute near-field force vector
DO j=1,system%no_particles
  jkind = system%particle(j)%parameters%kind
  IF (jkind=="sphere2" .AND. j/=i) THEN
    jradius = system%particle(j)%radius
    rj = R(j)%r
    rimrj = ri - rj
    norm_rimrj = SQRT(DOT_PRODUCT(rimrj,rimrj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (norm_rimrj<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
      !STOP
    END IF  
    nij = -(one/norm_rimrj)*rimrj
    fnf = fnf + kij_attr*(norm_rimrj**expij_attr)*nij - kij_rep*(norm_rimrj**expij_rep)*nij
  END IF
END DO  

END SUBROUTINE compute_sphere2_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere2_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i, j, np
CHARACTER(LEN=50)  :: j_kind
REAL(KIND=real8)   :: mass, charge, ci, norm_ri_m_rj, i_radius, j_radius, warning_dist, &
                      relaxation_activ_dist, bonding_activ_dist, pressure, norm_outward_dir
REAL(KIND=real8), DIMENSION(3) :: g, E, B, ri, vi, rj, ri_m_rj, nij, center, rc, tauij1, tauij2, outward_dir

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
i_radius = system%particle(i)%radius
mass = system%particle(i)%mass
charge = system%particle(i)%charge
ri = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get global damping parameters
ci = system%global_damping_coefficient

! add external fields contribution to environment force vector
fenv = mass*g + charge*E + charge*(vi .vector. B)

! set cell internal pressure value
pressure = 50.0_real8

! compute cell geometric center (MAY BE COMPUTED EXTERNALLY FROM HERE!)
center = 0.0_real8
np = 0
DO j=1,system%no_particles
  j_kind = system%particle(j)%parameters%kind
  IF (j_kind=="sphere2") THEN
    np = np + 1
    rj = R(j)%r
    center = center + rj
  END IF
END DO
center = (one/np)*center

! compute cell outward normal at particle i
outward_dir = 0.0_real8
DO j=1,system%no_particles
  j_kind = system%particle(j)%parameters%kind
  IF (j/=i .AND. j_kind=="sphere2") THEN
    j_radius = system%particle(j)%radius
    rj = R(j)%r
    ri_m_rj = ri - rj
    norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    bonding_activ_dist = (i_radius+j_radius)*1.1_real8
    IF (norm_ri_m_rj<=bonding_activ_dist) THEN
      warning_dist = MIN(1.0E-03*i_radius,1.0E-03*j_radius)
      IF (norm_ri_m_rj<=warning_dist) THEN 
        WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the NFs!"
        IF (norm_ri_m_rj<=0.00000000000001_real8) THEN
          WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
          STOP
        END IF  
      END IF  
      nij = -(one/norm_ri_m_rj)*ri_m_rj
      rc = ri - center
      tauij1 = rc .vector. nij
      tauij2 = nij .vector. tauij1
      outward_dir = outward_dir + tauij2
      !test = DOT_PRODUCT(rc,tauij2)
      !IF (test>0.0_real8) f = f + pressure*tauij2
      !IF (test<0.0_real8) f = f + pressure*(-tauij2)
      !IF (test==0.0_real8) f = 
    END IF
  END IF
END DO

! normalize outward normal
norm_outward_dir = SQRT(DOT_PRODUCT(outward_dir,outward_dir))
outward_dir = (one/norm_outward_dir)*outward_dir

! add internal pressure contribution to environment force vector
fenv = fenv + pressure*outward_dir

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_sphere2_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere2_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_sphere2_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE sphere2_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE biotissue_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_biotissue_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
INTEGER(KIND=int4) :: i, j, s, w, spr_propset
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_ri_m_rj, iradius, jradius, E, sigma_crit, &
                      L0, deltaL_crit, k, c, L, deltaL, warning_dist, aux, diw, relaxation_activ_dist
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, vj, ri0, rj0, ri0_m_rj0, ri_m_rj, nij, r1, ri1, n_out, vw

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri0 = system%particle(i)%xyz_coordinates
ri = R(i)%r
vi = V(i)%r

! get nearfields coefficients and exponents 
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep

! clear near-field force vector
fnf = 0.0_real8

! add springs contribution (bonding) to near-field force vector
springs_bonding: DO s=1,system%no_springs
  
  ! add contribution of springs that start at i (springs i-j)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0 ) THEN
    j = system%spring(s)%connectivities(2)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF

  ! add contribution of springs that end at i (springs j-i)
  IF (system%spring(s)%connectivities(2)==i) THEN
    j = system%spring(s)%connectivities(1)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF

  ! add contribution of wall springs that start at i (springs i-w)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
    w = system%spring(s)%connecting_wall_number
    r1 = system%wall(w)%point_position
    vw = system%wall(w)%velocity
    ri1 = r1 - ri
    n_out = system%wall(w)%outside_normal
    diw = DOT_PRODUCT(ri1,n_out)
    IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
    L = ABS(diw)
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      WRITE(error_unit,*) "CAUTION: negative L0 for wall spring", s, " ; L0 will be taken as the particle´s radius" 
      L0 = iradius
    END IF
    warning_dist = 1.0E-03*iradius
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: spring", s, " is extremely compressed against wall", w, "; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " length of spring", s, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = n_out !may be improved to change direction!
    aux = SQRT(DOT_PRODUCT(nij,nij))
    nij = (one/aux)*nij
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vw)
    END IF
  END IF

END DO springs_bonding

! add NFs contribution (NFs interactions with other biotissue_sphere particles) to near-field-force vector 
nfs_with_biotissue_spheres: IF (kij_attr/=0.0_real8 .OR. kij_rep/=0.0_real8) THEN

  ! add contribution from other biotissue_sphere particles
   DO j=1,system%no_particles
    jkind = system%particle(j)%parameters%kind
    IF (jkind=="biotissue_sphere") THEN
      rj = R(j)%r
      jradius = system%particle(j)%radius
      ri_m_rj = ri - rj
      norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
      warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
      IF (norm_ri_m_rj<=warning_dist) THEN 
        WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
        !STOP
      END IF  
      nij = -(one/norm_ri_m_rj)*ri_m_rj
      fnf = fnf + kij_attr*(norm_ri_m_rj**expij_attr)*nij - kij_rep*(norm_ri_m_rj**expij_rep)*nij
    END IF
  END DO

END IF nfs_with_biotissue_spheres

END SUBROUTINE compute_biotissue_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_biotissue_sphere_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i
REAL(KIND=real8) :: mass, charge, ci
REAL(KIND=real8), DIMENSION(3) :: ri, vi, g, E, B


! get particle data
i = particle_number
mass = system%particle(i)%mass
charge = system%particle(i)%charge
ri = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get global damping parameters
ci = system%global_damping_coefficient

! add external fields contribution to environment force vector
fenv = mass*g + charge*E + charge*(vi .vector. B)

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_biotissue_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_biotissue_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_biotissue_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE biotissue_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE lipid_membrane_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_lipid_membrane_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
INTEGER(KIND=int4) :: i, j, s, w, spr_propset
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_ri_m_rj, iradius, jradius, E, sigma_crit, &
                      L0, deltaL_crit, k, c, L, deltaL, warning_dist, aux, diw, relaxation_activ_dist
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, vj, ri0, rj0, ri0_m_rj0, ri_m_rj, nij, r1, ri1, n_out, vw

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri0 = system%particle(i)%xyz_coordinates
ri = R(i)%r
vi = V(i)%r

! get nearfields coefficients and exponents
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep

! clear near-field force vector
fnf = 0.0_real8

! add springs contribution (bonding) to near-field force vector
springs_bonding: DO s=1,system%no_springs
  
  ! add contribution of springs that start at i (springs i-j)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0 ) THEN
    j = system%spring(s)%connectivities(2)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF
  
  ! add contribution of springs that end at i (springs j-i)
  IF (system%spring(s)%connectivities(2)==i) THEN
    j = system%spring(s)%connectivities(1)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF
    
  ! add contribution of wall springs that start at i (springs i-w)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
    w = system%spring(s)%connecting_wall_number
    r1 = system%wall(w)%point_position
    vw = system%wall(w)%velocity
    ri1 = r1 - ri
    n_out = system%wall(w)%outside_normal
    diw = DOT_PRODUCT(ri1,n_out)
    IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
    L = ABS(diw)
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      WRITE(error_unit,*) "CAUTION: negative L0 for wall spring", s, " ; L0 will be taken as the particle´s radius" 
      L0 = iradius
    END IF
    warning_dist = 1.0E-03*iradius
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: spring", s, " is extremely compressed against wall", w, "; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " length of spring", s, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = n_out !may be improved to change direction!
    aux = SQRT(DOT_PRODUCT(nij,nij))
    nij = (one/aux)*nij
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vw)
    END IF
  END IF

END DO springs_bonding 

! add NFs contribution (NFs interactions with other lipid_membrane_sphere particles) to near-field-force vector 
nfs_with_biotissue_spheres: IF (kij_attr/=0.0_real8 .OR. kij_rep/=0.0_real8) THEN

  ! add contribution from other lipid_membrane_sphere particles
   DO j=1,system%no_particles
    jkind = system%particle(j)%parameters%kind
    IF (jkind=="lipid_membrane_sphere") THEN
      rj = R(j)%r
      jradius = system%particle(j)%radius
      ri_m_rj = ri - rj
      norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
      warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
      IF (norm_ri_m_rj<=warning_dist) THEN 
        WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
        !STOP
      END IF  
      nij = -(one/norm_ri_m_rj)*ri_m_rj
      fnf = fnf + kij_attr*(norm_ri_m_rj**expij_attr)*nij - kij_rep*(norm_ri_m_rj**expij_rep)*nij
    END IF
  END DO

END IF nfs_with_biotissue_spheres

END SUBROUTINE compute_lipid_membrane_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_lipid_membrane_sphere_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i, i_ps_number, j, np, s
CHARACTER(LEN=50)  :: j_kind
REAL(KIND=real8)   :: mass,charge, ci, pressure, area
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, g, E, B, rcg, rcgi, rij, normal_vector, tauij1, tauij2
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
mass = system%particle(i)%mass
charge = system%particle(i)%charge
i_ps_number = system%particle(i)%pressure_surface_number
ri = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get global damping parameters
ci = system%global_damping_coefficient

! add external fields contribution to environment force vector
fenv = mass*g + charge*E + charge*(vi .vector. B)

! add internal pressure contribution to environment force vector (if pressure is present)
pressure_force: IF (i_ps_number/=0) THEN

  ! get pressure data
  pressure = system%pressure_surface(i_ps_number)%initial_pressure
  area = system%pressure_surface(i_ps_number)%initial_influence_area
  
  ! compute CG of pressure surface (MAY BE COMPUTED EXTERNALLY FROM HERE...!)
  np = 0
  rcg = 0.0_real8
  DO j=1,system%no_particles
    IF (system%particle(j)%pressure_surface_number==i_ps_number) THEN
      np = np + 1
      rj = R(j)%r
      rcg = rcg + rj
    END IF
  END DO
  rcg = (one/np)*rcg

  ! clear normal vector of pressure surface at particle i
  normal_vector = 0.0_real8

  ! compute normal vector of pressure surface at particle i (from neighboring particles)
  neighbor_particles: DO s=1,system%no_springs
  
    ! compute vector connecting CG to particle i
    rcgi = ri - rcg
    
    ! compute normal vector w.r.t. neighbor particle j 
    IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0) THEN
      j = system%spring(s)%connectivities(2)
      rj = R(j)%r
      rij = rj - ri
      tauij1 = rcgi .vector. rij
      tauij2 = rij .vector. tauij1
      normal_vector = normal_vector + tauij2
    END IF
  
    ! compute normal vector w.r.t. neighbor particle j 
    IF (system%spring(s)%connectivities(2)==i) THEN
      j = system%spring(s)%connectivities(1)
      rj = R(j)%r
      rij = rj - ri
      tauij1 = rcgi .vector. rij
      tauij2 = rij .vector. tauij1
      normal_vector = normal_vector + tauij2
    END IF
    
    ! compute normal vector w.r.t. wall 
    !IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
      !w = system%spring(s)%connecting_wall_number
      !r1 = system%wall(w)%point1_coords
      !ri1 = r1 - ri
      !n_out = system%wall(w)%outside_normal
      !diw = DOT_PRODUCT(ri1,n_out)
      !IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
      !press_surf_number = system%spring(s)%pressure_surface_number
      !pressure = system%pressure_surface(press_surf_number)%initial_pressure
      !area = system%pressure_surface(press_surf_number)%initial_influence_area
      !j = system%spring(s)%connectivities(1)
      !j_radius = system%particle(j)%radius
      !rj = R(j)%r
      !rij = rj - ri
      !normal_dir = normal_dir + rij
    !END IF

  END DO neighbor_particles
  
  ! normalize normal vector
  normal_vector = (one/(SQRT(DOT_PRODUCT(normal_vector,normal_vector))))*normal_vector
  
  ! add internal pressure contribution to environment force vector
  IF (DOT_PRODUCT(rcgi,normal_vector)>0.0_real8) fenv = fenv + pressure*area*normal_vector
  IF (DOT_PRODUCT(rcgi,normal_vector)<0.0_real8) fenv = fenv + pressure*area*(-normal_vector)
  IF (DOT_PRODUCT(rcgi,normal_vector)==0.0_real8) WRITE (error_unit,*) " CAUTION: ambiguous normal vector at particle", i, " be aware of problems on the pressure forces!"

END IF pressure_force

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_lipid_membrane_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_lipid_membrane_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_lipid_membrane_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE lipid_membrane_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE internal_membrane_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_internal_membrane_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
CHARACTER(LEN=50)  :: jkind
INTEGER(KIND=int4) :: i, j, s, w, spr_propset
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_ri_m_rj, iradius, jradius, E, sigma_crit, &
                      L0, deltaL_crit, k, c, L, deltaL, warning_dist, aux, diw
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, vj, ri0, rj0, ri0_m_rj0, ri_m_rj, nij, r1, ri1, n_out, vw

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri0 = system%particle(i)%xyz_coordinates
ri = R(i)%r
vi = V(i)%r

! get nearfields coefficients and exponents (for interaction with external membrane only)
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep

! clear near-field force vector
fnf = 0.0_real8

! add springs contribution (bonding) to near-field force vector
springs_bonding: DO s=1,system%no_springs

  ! add contribution of springs that start at i (springs i-j)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0 ) THEN
    j = system%spring(s)%connectivities(2)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF
  
  ! add contribution of springs that end at i (springs j-i)
  IF (system%spring(s)%connectivities(2)==i) THEN
    j = system%spring(s)%connectivities(1)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF
    
  ! add contribution of wall springs that start at i (springs i-w)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
    w = system%spring(s)%connecting_wall_number
    r1 = system%wall(w)%point_position
    vw = system%wall(w)%velocity
    ri1 = r1 - ri
    n_out = system%wall(w)%outside_normal
    diw = DOT_PRODUCT(ri1,n_out)
    IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
    L = ABS(diw)
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      WRITE(error_unit,*) "CAUTION: negative L0 for wall spring", s, " ; L0 will be taken as the particle´s radius" 
      L0 = iradius
    END IF
    warning_dist = 1.0E-03*iradius
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: spring", s, " is extremely compressed against wall", w, "; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " length of spring", s, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = n_out !may be improved to allow change of direction!
    aux = SQRT(DOT_PRODUCT(nij,nij))
    nij = (one/aux)*nij
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vw)
    END IF
  END IF

END DO springs_bonding 

! add NFs contribution (NFs interactions with external membrane) to near-field-force vector 
nfs_with_external_membrane: IF (kij_attr/=0.0_real8 .OR. kij_rep/=0.0_real8) THEN

  ! add contribution from external particles
   DO j=1,system%no_particles
    jkind = system%particle(j)%parameters%kind
    IF (jkind=="external_membrane_sphere") THEN
      rj = R(j)%r
      jradius = system%particle(j)%radius
      ri_m_rj = ri - rj
      norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
      warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
      IF (norm_ri_m_rj<=warning_dist) THEN 
        WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
        !STOP
      END IF  
      nij = -(one/norm_ri_m_rj)*ri_m_rj
      fnf = fnf + kij_attr*(norm_ri_m_rj**expij_attr)*nij - kij_rep*(norm_ri_m_rj**expij_rep)*nij
    END IF
  END DO

END IF nfs_with_external_membrane

END SUBROUTINE compute_internal_membrane_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_internal_membrane_sphere_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i, i_ps_number, j, np, s
CHARACTER(LEN=50)  :: j_kind
REAL(KIND=real8)   :: mass, charge, ci, pressure, area
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, g, E, B, rcg, rcgi, rij, normal_vector, tauij1, tauij2

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
mass = system%particle(i)%mass
charge = system%particle(i)%charge
i_ps_number = system%particle(i)%pressure_surface_number
ri = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get global damping parameters
ci = system%global_damping_coefficient

! add external fields contribution to environment force vector
fenv = mass*g + charge*E + charge*(vi .vector. B)

! add internal pressure contribution to environment force vector (if pressure is present)
pressure_force: IF (i_ps_number/=0) THEN

  ! get pressure data
  pressure = system%pressure_surface(i_ps_number)%initial_pressure
  area = system%pressure_surface(i_ps_number)%initial_influence_area
  
  ! compute CG of pressure surface (MAY BE COMPUTED EXTERNALLY FROM HERE...!)
  np = 0
  rcg = 0.0_real8
  DO j=1,system%no_particles
    IF (system%particle(j)%pressure_surface_number==i_ps_number) THEN
      np = np + 1
      rj = R(j)%r
      rcg = rcg + rj
    END IF
  END DO
  rcg = (one/np)*rcg

  ! clear normal vector of pressure surface at particle i
  normal_vector = 0.0_real8

  ! compute normal vector of pressure surface at particle i (from neighboring particles)
  neighbor_particles: DO s=1,system%no_springs
  
    ! compute vector connecting CG to particle i
    rcgi = ri - rcg
    
    ! compute normal vector w.r.t. neighbor particle j 
    IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0) THEN
      j = system%spring(s)%connectivities(2)
      rj = R(j)%r
      rij = rj - ri
      tauij1 = rcgi .vector. rij
      tauij2 = rij .vector. tauij1
      normal_vector = normal_vector + tauij2
    END IF
  
    ! compute normal vector w.r.t. neighbor particle j 
    IF (system%spring(s)%connectivities(2)==i) THEN
      j = system%spring(s)%connectivities(1)
      rj = R(j)%r
      rij = rj - ri
      tauij1 = rcgi .vector. rij
      tauij2 = rij .vector. tauij1
      normal_vector = normal_vector + tauij2
    END IF
    
    ! compute normal vector w.r.t. wall 
    !IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
      !w = system%spring(s)%connecting_wall_number
      !r1 = system%wall(w)%point1_coords
      !ri1 = r1 - ri
      !n_out = system%wall(w)%outside_normal
      !diw = DOT_PRODUCT(ri1,n_out)
      !IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
      !press_surf_number = system%spring(s)%pressure_surface_number
      !pressure = system%pressure_surface(press_surf_number)%initial_pressure
      !area = system%pressure_surface(press_surf_number)%initial_influence_area
      !j = system%spring(s)%connectivities(1)
      !j_radius = system%particle(j)%radius
      !rj = R(j)%r
      !rij = rj - ri
      !normal_dir = normal_dir + rij
    !END IF

  END DO neighbor_particles
  
  ! normalize normal vector
  normal_vector = (one/(SQRT(DOT_PRODUCT(normal_vector,normal_vector))))*normal_vector
  
  ! add internal pressure contribution to environment force vector
  IF (DOT_PRODUCT(rcgi,normal_vector)>0.0_real8) fenv = fenv + pressure*area*normal_vector
  IF (DOT_PRODUCT(rcgi,normal_vector)<0.0_real8) fenv = fenv + pressure*area*(-normal_vector)
  IF (DOT_PRODUCT(rcgi,normal_vector)==0.0_real8) WRITE (error_unit,*) " CAUTION: ambiguous normal vector at particle", i, " be aware of problems on the pressure forces!"

END IF pressure_force

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_internal_membrane_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_internal_membrane_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_internal_membrane_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE internal_membrane_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE external_membrane_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_external_membrane_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
CHARACTER(LEN=50)  :: jkind
INTEGER(KIND=int4) :: i, j, s, w, spr_propset
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_ri_m_rj, iradius, jradius, E, sigma_crit, &
                      L0, deltaL_crit, k, c, L, deltaL, warning_dist, aux, diw
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, vj, ri0, rj0, ri0_m_rj0, ri_m_rj, nij, r1, ri1, n_out, vw

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri0 = system%particle(i)%xyz_coordinates
ri = R(i)%r
vi = V(i)%r

! get nearfields coefficients and exponents (for interaction with internal membrane only)
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep

! clear near-field force vector
fnf = 0.0_real8

! add springs contribution (bonding) to near-field force vector
springs_bonding: DO s=1,system%no_springs
  
  ! add contribution of springs that start at i (springs i-j)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0 ) THEN
    j = system%spring(s)%connectivities(2)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF
  
  ! add contribution of springs that end at i (springs j-i)
  IF (system%spring(s)%connectivities(2)==i) THEN
    j = system%spring(s)%connectivities(1)
    rj = R(j)%r
    vj = V(j)%r
    rj0 = system%particle(j)%xyz_coordinates
    jradius = system%particle(j)%radius
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      ri0_m_rj0 = ri0 - rj0
      L0 = SQRT(DOT_PRODUCT(ri0_m_rj0,ri0_m_rj0))
    END IF
    ri_m_rj = ri - rj
    L = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " distance between particles", i, " and", j, " is less than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = -(one/L)*ri_m_rj
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vj)
    END IF
  END IF
    
  ! add contribution of wall springs that start at i (springs i-w)
  IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
    w = system%spring(s)%connecting_wall_number
    r1 = system%wall(w)%point_position
    vw = system%wall(w)%velocity
    ri1 = r1 - ri
    n_out = system%wall(w)%outside_normal
    diw = DOT_PRODUCT(ri1,n_out)
    IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
    L = ABS(diw)
    spr_propset = system%spring(s)%properties_set_number
    k = system%spring_properties_set(spr_propset)%stiffness
    E = system%spring_properties_set(spr_propset)%elasticity_modulus
    sigma_crit = system%spring_properties_set(spr_propset)%critical_stress
    c = system%spring_properties_set(spr_propset)%dashpot_constant
    L0 = system%spring(s)%initial_length
    IF (L0<0) THEN
      WRITE(error_unit,*) "CAUTION: negative L0 for wall spring", s, " ; L0 will be taken as the particle´s radius" 
      L0 = iradius
    END IF
    warning_dist = 1.0E-03*iradius
    IF (L<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: spring", s, " is extremely compressed against wall", w, "; be aware of numerical errors on the spring forces!"
      IF (L<=0.00000000000001_real8) THEN
        WRITE(error_unit,*) " length of spring", s, " is smaller than 1.0E-14; program will be terminated to avoid bad numbers"
        STOP
      END IF  
    END IF  
    nij = n_out  !may be improved to allow change of direction!
    aux = SQRT(DOT_PRODUCT(nij,nij))
    nij = (one/aux)*nij
    deltaL = L - L0
    deltaL_crit = sigma_crit*L0/E
    IF (deltaL<=deltaL_crit) THEN
      fnf = fnf + k*deltaL*nij - c*(vi-vw)
    END IF
  END IF

END DO springs_bonding 


! add NFs contribution (NFs interactions with internal membrane) to near-field force vector
nfs_with_internal_membrane: IF (kij_attr/=0.0_real8 .OR. kij_rep/=0.0_real8) THEN

  ! add contribution from internal particles
  DO j=1,system%no_particles
    jkind = system%particle(j)%parameters%kind
    IF (jkind=="internal_membrane_sphere") THEN
      rj = R(j)%r
      jradius = system%particle(j)%radius
      ri_m_rj = ri - rj
      norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
      warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
      IF (norm_ri_m_rj<=warning_dist) THEN 
        WRITE(error_unit,*) " CAUTION: particles ", i, " and ", j, " are extremely overlapped; be aware of numerical errors on the NFs!"
        !STOP
      END IF  
      nij = -(one/norm_ri_m_rj)*ri_m_rj
      fnf = fnf + kij_attr*(norm_ri_m_rj**expij_attr)*nij - kij_rep*(norm_ri_m_rj**expij_rep)*nij
    END IF
  END DO

END IF nfs_with_internal_membrane

END SUBROUTINE compute_external_membrane_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_external_membrane_sphere_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i, i_ps_number, j, np, s
CHARACTER(LEN=50)  :: j_kind
REAL(KIND=real8)   :: mass, charge, ci, pressure, area
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, g, E, B, rcg, rcgi, rij, normal_vector, tauij1, tauij2

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
mass = system%particle(i)%mass
charge = system%particle(i)%charge
i_ps_number = system%particle(i)%pressure_surface_number
ri = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get global damping parameters
ci = system%global_damping_coefficient

! add external fields contribution to environment force vector
fenv = mass*g + charge*E + charge*(vi .vector. B)

! add internal pressure contribution to environment force vector (if pressure is present)
pressure_force: IF (i_ps_number/=0) THEN

  ! get pressure data
  pressure = system%pressure_surface(i_ps_number)%initial_pressure
  area = system%pressure_surface(i_ps_number)%initial_influence_area
  
  ! compute CG of pressure surface (MAY BE COMPUTED EXTERNALLY FROM HERE...!)
  np = 0
  rcg = 0.0_real8
  DO j=1,system%no_particles
    IF (system%particle(j)%pressure_surface_number==i_ps_number) THEN
      np = np + 1
      rj = R(j)%r
      rcg = rcg + rj
    END IF
  END DO
  rcg = (one/np)*rcg

  ! clear normal vector of pressure surface at particle i
  normal_vector = 0.0_real8

  ! compute normal vector of pressure surface at particle i (from neighboring particles)
  neighbor_particles: DO s=1,system%no_springs
  
    ! compute vector connecting CG to particle i
    rcgi = ri - rcg
    
    ! compute normal vector w.r.t. neighbor particle j 
    IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)/=0) THEN
      j = system%spring(s)%connectivities(2)
      rj = R(j)%r
      rij = rj - ri
      tauij1 = rcgi .vector. rij
      tauij2 = rij .vector. tauij1
      normal_vector = normal_vector + tauij2
    END IF
  
    ! compute normal vector w.r.t. neighbor particle j 
    IF (system%spring(s)%connectivities(2)==i) THEN
      j = system%spring(s)%connectivities(1)
      rj = R(j)%r
      rij = rj - ri
      tauij1 = rcgi .vector. rij
      tauij2 = rij .vector. tauij1
      normal_vector = normal_vector + tauij2
    END IF
    
    ! compute normal vector w.r.t. wall 
    !IF (system%spring(s)%connectivities(1)==i .AND. system%spring(s)%connectivities(2)==0) THEN
      !w = system%spring(s)%connecting_wall_number
      !r1 = system%wall(w)%point1_coords
      !ri1 = r1 - ri
      !n_out = system%wall(w)%outside_normal
      !diw = DOT_PRODUCT(ri1,n_out)
      !IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance for particle", i, ", spring", s, " and wall", w, " !"
      !press_surf_number = system%spring(s)%pressure_surface_number
      !pressure = system%pressure_surface(press_surf_number)%initial_pressure
      !area = system%pressure_surface(press_surf_number)%initial_influence_area
      !j = system%spring(s)%connectivities(1)
      !j_radius = system%particle(j)%radius
      !rj = R(j)%r
      !rij = rj - ri
      !normal_dir = normal_dir + rij
    !END IF

  END DO neighbor_particles
  
  ! normalize normal vector
  normal_vector = (one/(SQRT(DOT_PRODUCT(normal_vector,normal_vector))))*normal_vector
  
  ! add internal pressure contribution to environment force vector
  IF (DOT_PRODUCT(rcgi,normal_vector)>0.0_real8) fenv = fenv + pressure*area*normal_vector
  IF (DOT_PRODUCT(rcgi,normal_vector)<0.0_real8) fenv = fenv + pressure*area*(-normal_vector)
  IF (DOT_PRODUCT(rcgi,normal_vector)==0.0_real8) WRITE (error_unit,*) " CAUTION: ambiguous normal vector at particle", i, " be aware of problems on the pressure forces!"

END IF pressure_force

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_external_membrane_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_external_membrane_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_external_membrane_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE external_membrane_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE jet_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_jet_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
INTEGER(KIND=int4) :: i, j
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_rimrj, iradius, jradius, warning_dist, dcutoff
REAL(KIND=real8), DIMENSION(3) :: ri, vi, rj, rimrj, nij

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri = R(i)%r
vi = V(i)%r

! get nearfields data (coefficients, exponents and cutoff distance)
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep
dcutoff = system%nearfields_set(1)%cutoff_distance

! clear near-field force vector
fnf = 0.0_real8

! compute near-field force vector
DO j=1,system%no_particles
  jkind = system%particle(j)%parameters%kind
  IF (jkind=="jet_sphere" .AND. j/=i) THEN
    jradius = system%particle(j)%radius
    rj = R(j)%r
    rimrj = ri - rj
    norm_rimrj = SQRT(DOT_PRODUCT(rimrj,rimrj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (norm_rimrj<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
      !STOP
    END IF  
    IF (norm_rimrj<=dcutoff) THEN
      nij = -(one/norm_rimrj)*rimrj
      fnf = fnf + kij_attr*(norm_rimrj**expij_attr)*nij - kij_rep*(norm_rimrj**expij_rep)*nij
    END IF
  END IF
END DO  

END SUBROUTINE compute_jet_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_jet_sphere_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i
REAL(KIND=real8) :: mi, qi, ri, rho, eta, Ai, ci, viminusvf_norm, Re, Cd, Cstok, &
                    emfields_xi, emfields_xf, emfields_yi, emfields_yf, emfields_zi, emfields_zf
REAL(KIND=real8), DIMENSION(3) :: xi, vi, g, E, B, vf, viminusvf
REAL(KIND=real8), PARAMETER :: one=1.0_real8, pi=4.0_real8*ATAN(1.0_real8) 


! get particle data
i = particle_number
mi = system%particle(i)%mass
qi = system%particle(i)%charge
ri = system%particle(i)%radius
xi = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get electric and magnetic fields domain limits
emfields_xi = system%external_fields%emfields_xi
emfields_xf = system%external_fields%emfields_xf
emfields_yi = system%external_fields%emfields_yi
emfields_yf = system%external_fields%emfields_yf
emfields_zi = system%external_fields%emfields_zi
emfields_zf = system%external_fields%emfields_zf

! get environment fluid data
rho = system%environment_fluid%density
eta = system%environment_fluid%viscosity
vf = system%environment_fluid%velocity

! get global damping parameters
ci = system%global_damping_coefficient

! compute drag force frontal area
Ai = pi*ri*ri

! compute relative velocity between particle and environment fluid
viminusvf = vi - vf
viminusvf_norm = SQRT(DOT_PRODUCT(viminusvf,viminusvf))

! compute drag parameter in case of stokesian flow
Cstok = 6.0_real8*eta*pi*ri

! compute Reynold´s number
IF (eta/=0.0_real8) THEN
  Re = 2.0_real8*ri*rho*viminusvf_norm/eta
ELSE
  Re = -1.0_real8
END IF

! compute drag coefficient
IF (Re>0.1_real8 .AND. Re<=1.0_real8) Cd = 24.0_real8*(one/Re)
IF (Re>1.0_real8 .AND. Re<=400_real8) Cd = 24.0_real8*(one/Re)**0.646_real8
IF (Re>400.0_real8 .AND. Re<=3.0E+05) Cd = 0.5_real8
IF (Re>3.0E+05 .AND. Re<=2.0E+06) Cd = 0.000366*(Re)**0.4275_real8
IF (Re>2.0E+06) Cd = 0.18_real8

! add drag contribution to environment force vector
IF (Re<0.0_real8) THEN
  fenv = 0.0_real8
ELSE IF (Re>=0.0_real8 .AND. Re<=0.1_real8) THEN
  fenv = -Cstok*viminusvf
ELSE
  fenv = -0.5_real8*rho*Cd*Ai*viminusvf_norm*viminusvf
END IF

! add electric and magnetic fields contribution to environment force vector (only if within electromagnetic domain)
IF (xi(1)>=emfields_xi .AND. xi(1)<=emfields_xf) THEN
   IF (xi(2)>=emfields_yi .AND. xi(2)<=emfields_yf) THEN
     IF (xi(3)>=emfields_zi .AND. xi(3)<=emfields_zf) THEN
       fenv = fenv + qi*E + qi*(vi .vector. B)
     END IF
   END IF
END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Example 4.8 from LD:
!IF (xi(1)>=1.0_real8 .AND. xi(1)<=2.0_real8) THEN
!  fenv = fenv + qi*E + qi*(vi .vector. B)
!END IF
!fenv = fenv + mi*g
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! add gravity contribution to environment force vector
fenv = fenv + mi*g

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_jet_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_jet_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_jet_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE jet_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE granular_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_granular_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
INTEGER(KIND=int4) :: i, j, k, imat, jmat
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_rimrj, iradius, jradius, warning_dist, dcutoff
REAL(KIND=real8), DIMENSION(3) :: ri, rj, rimrj, nij, r1, ri1, n_out, vi, vj, vw, aux

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri = R(i)%r
vi = V(i)%r

! get nearfields data (coefficients, exponents and cutoff distance)
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep
dcutoff = system%nearfields_set(1)%cutoff_distance

! clear near-field force vector
fnf = 0.0_real8

! compute near-field force vector
DO j=1,system%no_particles
  jkind = system%particle(j)%parameters%kind
  IF (jkind=="granular_sphere" .AND. j/=i) THEN
    jradius = system%particle(j)%radius
    rj = R(j)%r
    rimrj = ri - rj
    norm_rimrj = SQRT(DOT_PRODUCT(rimrj,rimrj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (norm_rimrj<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
      !STOP
    END IF  
    IF (norm_rimrj<=dcutoff) THEN
      nij = -(one/norm_rimrj)*rimrj
      fnf = fnf + kij_attr*(norm_rimrj**expij_attr)*nij - kij_rep*(norm_rimrj**expij_rep)*nij
    END IF
  END IF
END DO

END SUBROUTINE compute_granular_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_granular_sphere_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i
REAL(KIND=real8) :: mi, qi, ri, rho, eta, Ai, ci, viminusvf_norm, Re, Cd, Cstok, &
                    emfields_xi, emfields_xf, emfields_yi, emfields_yf, emfields_zi, emfields_zf
REAL(KIND=real8), DIMENSION(3) :: xi, vi, g, E, B, vf, viminusvf
REAL(KIND=real8), PARAMETER :: one=1.0_real8, pi=4.0_real8*ATAN(1.0_real8) 


! get particle data
i = particle_number
mi = system%particle(i)%mass
qi = system%particle(i)%charge
ri = system%particle(i)%radius
xi = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get electric and magnetic fields domain limits
emfields_xi = system%external_fields%emfields_xi
emfields_xf = system%external_fields%emfields_xf
emfields_yi = system%external_fields%emfields_yi
emfields_yf = system%external_fields%emfields_yf
emfields_zi = system%external_fields%emfields_zi
emfields_zf = system%external_fields%emfields_zf

! get environment fluid data
rho = system%environment_fluid%density
eta = system%environment_fluid%viscosity
vf = system%environment_fluid%velocity

! get global damping parameters
ci = system%global_damping_coefficient

! compute drag force frontal area
Ai = pi*ri*ri

! compute relative velocity between particle and environment fluid
viminusvf = vi - vf
viminusvf_norm = SQRT(DOT_PRODUCT(viminusvf,viminusvf))

! compute drag parameter in case of stokesian flow
Cstok = 6.0_real8*eta*pi*ri

! compute Reynold´s number
IF (eta/=0.0_real8) THEN
  Re = 2.0_real8*ri*rho*viminusvf_norm/eta
ELSE
  Re = -1.0_real8
END IF

! compute drag coefficient
IF (Re>0.1_real8 .AND. Re<=1.0_real8) Cd = 24.0_real8*(one/Re)
IF (Re>1.0_real8 .AND. Re<=400_real8) Cd = 24.0_real8*(one/Re)**0.646_real8
IF (Re>400.0_real8 .AND. Re<=3.0E+05) Cd = 0.5_real8
IF (Re>3.0E+05 .AND. Re<=2.0E+06) Cd = 0.000366*(Re)**0.4275_real8
IF (Re>2.0E+06) Cd = 0.18_real8

! add drag contribution to environment force vector
IF (Re<0.0_real8) THEN
  fenv = 0.0_real8
ELSE IF (Re>=0.0_real8 .AND. Re<=0.1_real8) THEN
  fenv = -Cstok*viminusvf
ELSE
  fenv = -0.5_real8*rho*Cd*Ai*viminusvf_norm*viminusvf
END IF

! add electric and magnetic fields contribution to environment force vector (only if within electromagnetic domain)
IF (xi(1)>=emfields_xi .AND. xi(1)<=emfields_xf) THEN
   IF (xi(2)>=emfields_yi .AND. xi(2)<=emfields_yf) THEN
     IF (xi(3)>=emfields_zi .AND. xi(3)<=emfields_zf) THEN
       fenv = fenv + qi*E + qi*(vi .vector. B)
     END IF
   END IF
END IF

! add gravity contribution to environment force vector
fenv = fenv + mi*g

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

END SUBROUTINE compute_granular_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_granular_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_granular_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE granular_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE granular_sphere_in_given_fluid_flow_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_granular_sphere_igff_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
INTEGER(KIND=int4) :: i, j, k, imat, jmat
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_rimrj, iradius, jradius, warning_dist, dcutoff
REAL(KIND=real8), DIMENSION(3) :: ri, rj, rimrj, nij, r1, ri1, n_out, vi, vj, vw, aux

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri = R(i)%r
vi = V(i)%r

! get nearfields data (coefficients, exponents and cutoff distance)
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep
dcutoff = system%nearfields_set(1)%cutoff_distance

! clear near-field force vector
fnf = 0.0_real8

! compute near-field force vector
DO j=1,system%no_particles
  jkind = system%particle(j)%parameters%kind
  IF (jkind=="granular_sphere_inn_given_fluid_flow" .AND. j/=i) THEN
    jradius = system%particle(j)%radius
    rj = R(j)%r
    rimrj = ri - rj
    norm_rimrj = SQRT(DOT_PRODUCT(rimrj,rimrj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (norm_rimrj<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
      !STOP
    END IF  
    IF (norm_rimrj<=dcutoff) THEN
      nij = -(one/norm_rimrj)*rimrj
      fnf = fnf + kij_attr*(norm_rimrj**expij_attr)*nij - kij_rep*(norm_rimrj**expij_rep)*nij
    END IF
  END IF
END DO

END SUBROUTINE compute_granular_sphere_igff_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_granular_sphere_igff_environment_forces(particle_number,R,V,system,fenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv

! local variables
INTEGER(KIND=int4) :: i, ncn, iv, j
REAL(KIND=real8) :: mi, qi, ri, rho, eta, Ai, ci, viminusvflocal_norm, Re, Cd, Cstok, emfields_xi, emfields_xf, emfields_yi, emfields_yf,  &
                    emfields_zi, emfields_zf, vf_inf, rj, rpolar, costheta, sintheta, cosphi, sinphi, vr, vtheta, pressure_gradient
REAL(KIND=real8), DIMENSION(3) :: xi, vi, g, E, B, vf, vf_local, viminusvflocal, xj, vf_local_j
REAL(KIND=real8), PARAMETER :: one=1.0_real8, pi=4.0_real8*ATAN(1.0_real8) 


! get particle data
i = particle_number
mi = system%particle(i)%mass
qi = system%particle(i)%charge
ri = system%particle(i)%radius
xi = R(i)%r
vi = V(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get electric and magnetic fields domain limits
emfields_xi = system%external_fields%emfields_xi
emfields_xf = system%external_fields%emfields_xf
emfields_yi = system%external_fields%emfields_yi
emfields_yf = system%external_fields%emfields_yf
emfields_zi = system%external_fields%emfields_zi
emfields_zf = system%external_fields%emfields_zf

! get environment fluid data
rho = system%environment_fluid%density
eta = system%environment_fluid%viscosity
vf = system%environment_fluid%velocity

! get global damping parameters
ci = system%global_damping_coefficient

! compute drag force frontal area
Ai = pi*ri*ri

! compute local velocity of environment fluid (an IF would be better here: IF sys%cv%compute_fluid_vel=="yes" THEN...)
vf_local = vf
!CALL compute_local_fluid_velocity

! compute relative velocity between particle and environment fluid
viminusvflocal = vi - vf_local
viminusvflocal_norm = SQRT(DOT_PRODUCT(viminusvflocal,viminusvflocal))

! compute drag parameter in case of stokesian flow
Cstok = 6.0_real8*eta*pi*ri

! compute Reynold´s number
IF (eta/=0.0_real8) THEN
  Re = 2.0_real8*ri*rho*viminusvflocal_norm/eta
ELSE
  Re = -1.0_real8
END IF

! compute drag coefficient
IF (Re>0.1_real8 .AND. Re<=1.0_real8) Cd = 24.0_real8*(one/Re)
IF (Re>1.0_real8 .AND. Re<=400_real8) Cd = 24.0_real8*(one/Re)**0.646_real8
IF (Re>400.0_real8 .AND. Re<=3.0E+05) Cd = 0.5_real8
IF (Re>3.0E+05 .AND. Re<=2.0E+06) Cd = 0.000366*(Re)**0.4275_real8
IF (Re>2.0E+06) Cd = 0.18_real8

! add drag contribution to environment force vector
IF (Re<0.0_real8) THEN
  fenv = 0.0_real8
ELSE IF (Re>=0.0_real8 .AND. Re<=0.1_real8) THEN
  fenv = -Cstok*viminusvflocal
ELSE
  fenv = -0.5_real8*rho*Cd*Ai*viminusvflocal_norm*viminusvflocal
END IF

! add electric and magnetic fields contribution to environment force vector (only if within electromagnetic domain)
IF (xi(1)>=emfields_xi .AND. xi(1)<=emfields_xf) THEN
   IF (xi(2)>=emfields_yi .AND. xi(2)<=emfields_yf) THEN
     IF (xi(3)>=emfields_zi .AND. xi(3)<=emfields_zf) THEN
       fenv = fenv + qi*E + qi*(vi .vector. B)
     END IF
   END IF
END IF

! add gravity contribution to environment force vector
fenv = fenv + mi*g

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

! add pressure gradient contribution to environment force vector (for infiltration paper (joint paper with Tarek in IJMCE-2019))
!(pressure gradient = pressure drop per unit length; only x-direction is considered here)
!IF (xi(1)<=0.05_real8) THEN
!  pressure_gradient = 0.01E+03
!ELSE
!  pressure_gradient = 0.01E+03 !0.1E+03  !0.04E+03  !0.145E+03
!END IF
!fenv(1) = fenv(1) + pressure_gradient*(2.0_real8*ri)*Ai


CONTAINS
  SUBROUTINE compute_local_fluid_velocity

    ! NOTICE: first it is necessary to identify the neighbor particles of i that are fixed, i.e., those around which the fluid flows
    ! This is done with a loop over the close neighbors of i (the fixed particles are those that influence the fluid velocity field at i´s location)

    ! get fluid velocity at infinity (i.e. far from the fixed particles´ location) (flow is assumed only in global x-direction)
    vf_inf = system%environment_fluid%velocity(1)
    IF (system%environment_fluid%velocity(1)==0.0_real8) WRITE(error_unit,*) " CAUTION: fluid velocity is not in x-direction in input file!"
   
    ! get number of close neighbors of particle i
    ncn = system%particle(i)%no_particles_into_verlet_list

    ! initialize fluid local velocity
    vf_local = 0.0_real8
    
    ! loop over close neighbors of particle i
    DO iv=1,ncn

      ! get close neighbor´s number
      j = system%particle(i)%verlet_list(iv)

      ! compute local velocity due to fixed sphere j (BELOW IS A 2D INVISCID FLUID FLOW PAST A CYLINDER) (flow in xy-plane)
      !IF (system%particle(j)%parameters%kind=="granular_sphere") THEN
      !  rj = system%particle(j)%radius
      !  xj = R(j)%r
      !  rpolar = SQRT( (xj(1)-xi(1))**2 + (xj(2)-xi(2))**2 )
      !  costheta = (xi(1)-xj(1))/rpolar
      !  sintheta = (xi(2)-xj(2))/rpolar
      !  vr = vf_inf*(1.0_real8-(rj/rpolar)**2)*costheta
      !  vtheta = - vf_inf*(1.0_real8+(rj/rpolar)**2)*sintheta
      !  vf_local_j(1) = vr*costheta - vtheta*sintheta
      !  vf_local_j(2) = vr*sintheta + vtheta*costheta
      !  vf_local_j(3) = 0.0_real8
      !  vf_local = vf_local + vf_local_j 
      !END IF
        
      ! compute local velocity due to fixed sphere j (BELOW IS A 3D VISCOUS FLUID FLOW PAST A SPHERE)
      IF (system%particle(j)%parameters%kind=="granular_sphere") THEN
        rj = system%particle(j)%radius
        xj = R(j)%r
        rpolar = SQRT((xj(1)-xi(1))**2 + (xj(2)-xi(2))**2 + (xj(3)-xi(3))**2)
        costheta = (xi(1)-xj(1))/rpolar
        sintheta = SQRT((xj(2)-xi(2))**2 + (xj(3)-xi(3))**2)/rpolar
        cosphi = (xi(2)-xj(2))/(rpolar*sintheta)
        sinphi = (xi(3)-xj(3))/(rpolar*sintheta)
        vr = vf_inf*(1.0_real8 - 1.5_real8*(rj/rpolar) + 0.5_real8*(rj/rpolar)**3)*costheta
        vtheta = - vf_inf*(1.0_real8 - 0.75_real8*(rj/rpolar) - 0.25_real8*(rj/rpolar)**3)*sintheta
        vf_local_j(1) = vr*costheta - vtheta*sintheta
        vf_local_j(2) = vr*sintheta*cosphi + vtheta*costheta*cosphi
        vf_local_j(3) = vr*sintheta*sinphi + vtheta*costheta*sinphi
        vf_local = vf_local + vf_local_j 
      END IF
    
    END DO

    ! check for consistency
    IF (ALL(vf_local==0.0_real8)) vf_local = system%environment_fluid%velocity

  END SUBROUTINE compute_local_fluid_velocity
  

END SUBROUTINE compute_granular_sphere_igff_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_granular_sphere_igff_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_granular_sphere_igff_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


END MODULE granular_sphere_in_given_fluid_flow_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE thermomechanical_sphere_class

USE particle_data_types
USE system_data_types
USE pefmat
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_thermomechanical_sphere_nearfield_forces(particle_number,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fnf

! local variables
INTEGER(KIND=int4) :: i, j, k, imat, jmat
CHARACTER(LEN=50)  :: jkind
REAL(KIND=real8)   :: kij_attr, expij_attr, kij_rep, expij_rep, norm_rimrj, iradius, jradius, warning_dist, dcutoff
REAL(KIND=real8), DIMENSION(3) :: ri, rj, rimrj, nij, r1, ri1, n_out, vi, vj, vw, aux

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8


! get particle data
i = particle_number
iradius = system%particle(i)%radius
ri = R(i)%r
vi = V(i)%r

! get nearfields data (coefficients, exponents and cutoff distance)
kij_attr = system%nearfields_set(1)%kij_attr
expij_attr = system%nearfields_set(1)%expij_attr
kij_rep = system%nearfields_set(1)%kij_rep
expij_rep = system%nearfields_set(1)%expij_rep
dcutoff = system%nearfields_set(1)%cutoff_distance

! clear near-field force vector
fnf = 0.0_real8

! compute near-field force vector
DO j=1,system%no_particles
  jkind = system%particle(j)%parameters%kind
  IF (jkind=="thermomechanical_sphere" .AND. j/=i) THEN
    jradius = system%particle(j)%radius
    rj = R(j)%r
    rimrj = ri - rj
    norm_rimrj = SQRT(DOT_PRODUCT(rimrj,rimrj))
    warning_dist = MIN(1.0E-03*iradius,1.0E-03*jradius)
    IF (norm_rimrj<=warning_dist) THEN 
      WRITE(error_unit,*) " CAUTION: particles", i, " and", j, " are extremely overlapped; be aware of numerical errors on nearfield forces!"
      !STOP
    END IF  
    IF (norm_rimrj<=dcutoff) THEN
      nij = -(one/norm_rimrj)*rimrj
      fnf = fnf + kij_attr*(norm_rimrj**expij_attr)*nij - kij_rep*(norm_rimrj**expij_rep)*nij
    END IF
  END IF
END DO

END SUBROUTINE compute_thermomechanical_sphere_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_thermomechanical_sphere_environment_forces(particle_number,R,V,system,fenv,menv,W)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V, W
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(1:3) :: fenv, menv

! local variables
INTEGER(KIND=int4) :: i
REAL(KIND=real8) :: mi, qi, ri, rho, eta, Ai, ci, viminusvf_norm, Re, Cd, Cstok, Cm, viminusvft_norm, &
                    emfields_xi, emfields_xf, emfields_yi, emfields_yf, emfields_zi, emfields_zf
REAL(KIND=real8), DIMENSION(3) :: xi, vi, g, E, B, vf, viminusvf, wi, vft, viminusvft, fdragt
REAL(KIND=real8), PARAMETER :: one=1.0_real8, pi=4.0_real8*ATAN(1.0_real8) 


! get particle data
i = particle_number
mi = system%particle(i)%mass
qi = system%particle(i)%charge
ri = system%particle(i)%radius
xi = R(i)%r
vi = V(i)%r
wi = W(i)%r

! get external fields data
g = system%external_fields%gravity_accel_vector
E = system%external_fields%electric_field_vector
B = system%external_fields%magnetic_field_vector

! get electric and magnetic fields domain limits
emfields_xi = system%external_fields%emfields_xi
emfields_xf = system%external_fields%emfields_xf
emfields_yi = system%external_fields%emfields_yi
emfields_yf = system%external_fields%emfields_yf
emfields_zi = system%external_fields%emfields_zi
emfields_zf = system%external_fields%emfields_zf

! get environment fluid data
rho = system%environment_fluid%density
eta = system%environment_fluid%viscosity
vf = system%environment_fluid%velocity

!xxxxxxxxxxxxxxxx modification in vf for fragments´paper xxxxxxxxxxxxxxxxxx
!vf = 0.0_real8
!vf(1) = 10.0_real8*xi(2)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! get global damping parameters
ci = system%global_damping_coefficient

! compute drag force frontal area
Ai = pi*ri*ri

! compute relative velocity between particle and environment fluid
viminusvf = vi - vf
viminusvf_norm = SQRT(DOT_PRODUCT(viminusvf,viminusvf))

! compute drag parameter in case of stokesian flow
Cstok = 6.0_real8*eta*pi*ri

! compute Reynold´s number
IF (eta/=0.0_real8) THEN
  Re = 2.0_real8*ri*rho*viminusvf_norm/eta
ELSE
  Re = -1.0_real8
END IF

! compute drag coefficient
IF (Re>0.1_real8 .AND. Re<=1.0_real8) Cd = 24.0_real8*(one/Re)
IF (Re>1.0_real8 .AND. Re<=400_real8) Cd = 24.0_real8*(one/Re)**0.646_real8
IF (Re>400.0_real8 .AND. Re<=3.0E+05) Cd = 0.5_real8
IF (Re>3.0E+05 .AND. Re<=2.0E+06) Cd = 0.000366*(Re)**0.4275_real8
IF (Re>2.0E+06) Cd = 0.18_real8

! add drag contribution to environment force vector
IF (Re<0.0_real8) THEN
  fenv = 0.0_real8
ELSE IF (Re>=0.0_real8 .AND. Re<=0.1_real8) THEN
  fenv = -Cstok*viminusvf
ELSE
  fenv = -0.5_real8*rho*Cd*Ai*viminusvf_norm*viminusvf
END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx compute magnus force (for fragments´paper) xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! add magnus force contribution to environment force vector
!Cm = 1.0_real8
!wi = system%particle(i)%spin  !<=CAUTION: it would be better to pass wi as argument (as done with R and V)
!fenv = fenv + 0.5_real8*Cm*rho*Ai*ri*(wi .vector. viminusvf)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! add electric and magnetic fields contribution to environment force vector (only if within electromagnetic domain)
IF (xi(1)>=emfields_xi .AND. xi(1)<=emfields_xf) THEN
   IF (xi(2)>=emfields_yi .AND. xi(2)<=emfields_yf) THEN
     IF (xi(3)>=emfields_zi .AND. xi(3)<=emfields_zf) THEN
       fenv = fenv + qi*E + qi*(vi .vector. B)
     END IF
   END IF
END IF

! add gravity contribution to environment force vector
fenv = fenv + mi*g

! add global damping contribution to environment force vector
fenv = fenv - ci*vi

! add fluid contribution to environment moment vector (e.g. in case drag is eccentric)
menv = 0.0_real8

! add global damping contribution to environment moment vector (obs: since the same "ci" is used here as for the translational velocity, its necessary to multiply by the particle radius in order attain units consistency!) 
menv = menv - ci*ri*wi

!xxxxxxxxxxxxxxxxxxxx compute moment due to fluid for fragments´paper xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! compute fluid´s additional velocity at top part of the particle (for Tarek´s fragments example)
!vft = 0.0_real8
!vft(1) = 10.0_real8*(xi(2)+0.5_real8*ri) - vf(1)

! compute relative velocities at top and bottom parts
!viminusvft = vi - vft
!viminusvft_norm = SQRT(DOT_PRODUCT(viminusvft,viminusvft))

! compute additional part of drag force at top of the particle
!IF (Re<0.0_real8) THEN
!  fdragt = 0.0_real8
!ELSE IF (Re>=0.0_real8 .AND. Re<=0.1_real8) THEN
!  fdragt = -Cstok*viminusvft
!ELSE
!  !fdragt = -0.5_real8*rho*Cd*(0.5_real8*Ai)*viminusvft_norm*viminusvft
!  fdragt = 0.0_real8
!  fdragt(1) = 0.5_real8*rho*Cd*(0.5_real8*Ai)*vft(1)*vft(1)
!END IF

! compute moment due to environment fluid (for Tarek´s fragments example)
!menv = 0.0_real8
!menv(3) = -fdragt(1)*ri
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END SUBROUTINE compute_thermomechanical_sphere_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_thermomechanical_sphere_environment_heat_power(particle_number,R,V,T,system,qpenv)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V, T
TYPE(system_data) :: system
REAL(KIND=real8)  :: qpenv

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: ri, ti, Ki, gammai, epsiloni, tenv, kenv, cenv, rho, eta, etas, Ai, Ais, vfminusvi_norm, Re, Cd, Cstok, Pr, Nu, hi
REAL(KIND=real8), DIMENSION(3) :: xi, vi, vf, vfminusvi, fdrag
REAL(KIND=real8), PARAMETER :: one=1.0_real8, pi=4.0_real8*ATAN(1.0_real8), twothirds=0.6666666666666666_real8, St_Boltz=5.670367E-08


! get particle data
i = particle_number
imat = system%particle(i)%material_set_number
Ki = system%material_set(imat)%thermal_conductivity
gammai = system%material_set(imat)%drag_heating_efficiency
epsiloni = system%material_set(imat)%radiative_efficiency
ri = system%particle(i)%radius
xi = R(i)%r
vi = V(i)%r
ti = T(i)%r(1)

! get environment fluid data
tenv = system%environment_fluid%temperature
kenv = system%environment_fluid%thermal_conductivity
cenv = system%environment_fluid%specific_heat
rho = system%environment_fluid%density
eta = system%environment_fluid%viscosity
vf = system%environment_fluid%velocity

!xxxxxxxxxxxxxxxx modification in vf for fragments´ paper xxxxxxxxxxxxxxxxx
!vf = 0.0_real8
!vf(1) = 10.0_real8*xi(2)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute particle´s frontal area
Ai = pi*ri*ri

! compute particle´s surface area
Ais = 4.0_real8*Ai

! compute relative velocity between particle and environment fluid
vfminusvi = vf - vi
vfminusvi_norm = SQRT(DOT_PRODUCT(vfminusvi,vfminusvi))

! compute drag parameter in case of stokesian flow
Cstok = 6.0_real8*eta*pi*ri

! compute Reynold´s number
IF (eta/=0.0_real8) THEN
  Re = 2.0_real8*ri*rho*vfminusvi_norm/eta
ELSE
  Re = -1.0_real8
END IF

! compute drag coefficient
IF (Re>0.1_real8 .AND. Re<=1.0_real8) Cd = 24.0_real8*(one/Re)
IF (Re>1.0_real8 .AND. Re<=400_real8) Cd = 24.0_real8*(one/Re)**0.646_real8
IF (Re>400.0_real8 .AND. Re<=3.0E+05) Cd = 0.5_real8
IF (Re>3.0E+05 .AND. Re<=2.0E+06) Cd = 0.000366*(Re)**0.4275_real8
IF (Re>2.0E+06) Cd = 0.18_real8

! compute drag force
IF (Re<0.0_real8) THEN
  fdrag = 0.0_real8
ELSE IF (Re>=0.0_real8 .AND. Re<=0.1_real8) THEN
  fdrag = Cstok*vfminusvi
ELSE
  fdrag = 0.5_real8*rho*Cd*Ai*vfminusvi_norm*vfminusvi
END IF

! compute Prandtl number
IF (kenv/=0.0_real8) THEN
  Pr = cenv*eta/kenv
ELSE
  Pr = -1.0_real8
END IF  

! compute Nusselt number
IF ((Re>=0.0_real8 .AND. eta/=0.0_real8) .AND. Pr>=0.0_real8) THEN
  etas = eta
  Nu = 2.0_real8 + (0.4_real8*SQRT(Re) + 0.06_real8*Re**twothirds)*(Pr**0.4_real8)*((eta/etas)**0.25_real8)
ELSE 
  Nu = 0.0_real8
END IF

! compute convection coefficient
hi = Nu*kenv/(2.0_real8*ri)

! add convective contribution to environment heat power
qpenv = hi*Ais*(tenv-ti)

! add drag-induced heating contribution to environment heat power
qpenv = qpenv + gammai*DOT_PRODUCT(fdrag,vfminusvi)

! add radiative contribution to environment heat power
qpenv = qpenv + epsiloni*St_Boltz*Ais*(tenv**4 - ti**4)

END SUBROUTINE compute_thermomechanical_sphere_environment_heat_power
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_thermomechanical_sphere_external_device_heat_power(particle_number,R,V,system,qpext)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8)  :: qpext

! local variables
INTEGER(KIND=int4) :: i, j, k, imat, jmat
CHARACTER(LEN=50)  :: jkind, devkind
REAL(KIND=real8)   :: norm_rimrj, ri, jradius, warning_dist, ai, I0, devarea, areai, zi, mi, dmax, Idev, devradius, diproj_norm
REAL(KIND=real8), DIMENSION(3) :: xi, rj, rimrj, nij, r1, ri1, n_out, vi, vj, vw, aux, xdev, ndev, vdev, e1, e2, di, diproj

! local parameters
REAL(KIND=real8), PARAMETER :: one=1.0_real8, pi=4.0_real8*ATAN(1.0_real8)


! get particle´s data
i = particle_number
ri = system%particle(i)%radius
ai = system%particle(i)%absorptance
xi = R(i)%r
!vi = V(i)%r

! clear external heat power
qpext = 0.0_real8

! compute external heat power
DO j=1,system%no_external_heating_devices

  ! get external device´s data (kind, intensity [power per unit area], etc)
  devkind = system%external_heating_device(j)%kind
  I0 = system%external_heating_device(j)%intensity
  devarea = system%external_heating_device(j)%cross_sectional_area
  xdev = system%external_heating_device(j)%striking_position
  ndev = system%external_heating_device(j)%axial_direction
  vdev = system%external_heating_device(j)%velocity
  mi = system%external_heating_device(j)%attenuation_coeff
  dmax = system%external_heating_device(j)%max_penetration

  ! check for consistency
  IF (devkind=="fire_nozzle" .AND. dmax==0.0_real8) THEN
    WRITE(error_unit,*) "external heating device´s max_penetration is zero; program will be terminated to avoid division by zero"
    STOP
  END IF    
  
  ! enforce unit axial direction and compute device´s cross-sectional directions (ATTENTION: vdev cannot be coaxial to ndev)
  ndev = one/SQRT(DOT_PRODUCT(ndev,ndev))*ndev
  e2 = (ndev .vector. vdev)
  e2 = one/SQRT(DOT_PRODUCT(e2,e2))*e2
  e1 = e2 .vector. ndev

  ! compute distance between particle and device´s striking position and project it onto the device´s cross-sectional plane
  di = xdev-xi
  diproj = DOT_PRODUCT(di,e1)*e1 + DOT_PRODUCT(di,e2)*e2
  diproj_norm = SQRT(DOT_PRODUCT(diproj,diproj))
  
  ! compute device´s input power to particle i
  Idev = 0.0_real8
  devradius = SQRT(devarea/pi)
  IF (diproj_norm<=devradius) THEN
    areai = pi*ri*ri
    zi = ABS(DOT_PRODUCT(di,ndev))
    IF (devkind=="fire_nozzle") THEN
      IF (zi<=dmax) Idev = I0*areai*(one-zi/dmax)
    ELSE IF (devkind=="laser_beam") THEN
      Idev = I0*areai*exp(-mi*zi)
    ELSE
      WRITE(error_unit,*) "external heating device kind not recognized or not yet implemented"
      STOP
    END IF
  END IF  

  ! add device´s contribution to external heat power
  qpext = qpext + ai*Idev
  
END DO

END SUBROUTINE compute_thermomechanical_sphere_external_device_heat_power
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_thermomechanical_sphere_mass_and_inertia(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8) :: iradius, imass_density, ivol
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
iradius = system%particle(i)%radius
imat = system%particle(i)%material_set_number
imass_density = system%material_set(imat)%mass_density

! compute particle volume
ivol = 1.333333333333333_real8*pi*(iradius**3)

! compute particle mass
system%particle(i)%mass = imass_density*ivol

! compute particle inertia
system%particle(i)%inertia = 0.4_real8*system%particle(i)%mass*(iradius**2)

END SUBROUTINE compute_thermomechanical_sphere_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_thermomechanical_sphere_thermal_properties(particle_number,system)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: particle_number
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, imat
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get particle data
i = particle_number
imat = system%particle(i)%material_set_number

! compute particle´s thermal properties
system%particle(i)%specific_heat = system%material_set(imat)%specific_heat
system%particle(i)%absorptance = system%material_set(imat)%absorptance

END SUBROUTINE compute_thermomechanical_sphere_thermal_properties
!!---------------------------------------------------------------------------------------------------------------------


END MODULE thermomechanical_sphere_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE particles_classes

! Insert here new particle classes (each particle class is to be defined as a separate module)
USE sphere1_class
USE sphere2_class
USE biotissue_sphere_class
USE lipid_membrane_sphere_class
USE internal_membrane_sphere_class
USE external_membrane_sphere_class
USE jet_sphere_class
USE granular_sphere_class
USE granular_sphere_in_given_fluid_flow_class
USE thermomechanical_sphere_class
USE particle_data_types
USE system_data_types
USE pefmat

IMPLICIT NONE

CONTAINS

!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_particle_nearfield_forces(particle_number,particle_kind,R,V,system,fnf)

! dummy arguments
INTEGER(KIND=int4) :: particle_number
CHARACTER(LEN=50), INTENT(IN) :: particle_kind
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(3) :: fnf


! select particle kind
IF (system%general_cv%nearfields_switch=="on") THEN
  
  SELECT CASE (particle_kind)
    
    CASE ("sphere1")
      CALL compute_sphere1_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("sphere2")
      CALL compute_sphere2_nearfield_forces(particle_number,R,V,system,fnf)
    
    CASE ("biotissue_sphere")
      CALL compute_biotissue_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("lipid_membrane_sphere")
      CALL compute_lipid_membrane_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("internal_membrane_sphere")
      CALL compute_internal_membrane_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("external_membrane_sphere")
      CALL compute_external_membrane_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("jet_sphere")
      CALL compute_jet_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("granular_sphere")
      CALL compute_granular_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("granular_sphere_inn_given_fluid_flow")
      CALL compute_granular_sphere_igff_nearfield_forces(particle_number,R,V,system,fnf)

    CASE ("thermomechanical_sphere")
      CALL compute_thermomechanical_sphere_nearfield_forces(particle_number,R,V,system,fnf)

    CASE DEFAULT
      WRITE(error_unit,*) "particle kind not recognized or not yet implemented"
      STOP
  
  END SELECT
  
ELSE

  fnf = 0.0_real8

END IF

END SUBROUTINE compute_particle_nearfield_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_particle_environment_forces(particle_number,particle_kind,R,V,system,fenv,menv,W)

! dummy arguments
INTEGER(KIND=int4) :: particle_number
CHARACTER(LEN=50), INTENT(IN) :: particle_kind
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(3) :: fenv
REAL(KIND=real8), DIMENSION(3), OPTIONAL :: menv
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: W

! select particle kind
SELECT CASE (particle_kind)
  
  CASE ("sphere1")
    CALL compute_sphere1_environment_forces(particle_number,R,V,system,fenv)
  
  CASE ("sphere2")
    CALL compute_sphere2_environment_forces(particle_number,R,V,system,fenv)
  
  CASE ("biotissue_sphere")
    CALL compute_biotissue_sphere_environment_forces(particle_number,R,V,system,fenv)

  CASE ("lipid_membrane_sphere")
    CALL compute_lipid_membrane_sphere_environment_forces(particle_number,R,V,system,fenv)

  CASE ("internal_membrane_sphere")
    CALL compute_internal_membrane_sphere_environment_forces(particle_number,R,V,system,fenv)

  CASE ("external_membrane_sphere")
    CALL compute_external_membrane_sphere_environment_forces(particle_number,R,V,system,fenv)

  CASE ("jet_sphere")
    CALL compute_jet_sphere_environment_forces(particle_number,R,V,system,fenv)

  CASE ("granular_sphere")
    CALL compute_granular_sphere_environment_forces(particle_number,R,V,system,fenv)

  CASE ("granular_sphere_in_given_fluid_flow")
    CALL compute_granular_sphere_igff_environment_forces(particle_number,R,V,system,fenv)

  CASE ("thermomechanical_sphere")
    CALL compute_thermomechanical_sphere_environment_forces(particle_number,R,V,system,fenv,menv,W)

  CASE DEFAULT
    WRITE(error_unit,*) "particle kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_particle_environment_forces
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_particle_environment_heat_power(particle_number,particle_kind,R,V,T,system,qpenv)

! dummy arguments
INTEGER(KIND=int4) :: particle_number
CHARACTER(LEN=50), INTENT(IN) :: particle_kind
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V, T
TYPE(system_data) :: system
REAL(KIND=real8)  :: qpenv


! select particle kind
SELECT CASE (particle_kind)
  
  CASE ("thermomechanical_sphere")
    CALL compute_thermomechanical_sphere_environment_heat_power(particle_number,R,V,T,system,qpenv)

  CASE DEFAULT
    WRITE(error_unit,*) "particle kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_particle_environment_heat_power
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_particle_external_device_heat_power(particle_number,particle_kind,R,V,system,qpext)

! dummy arguments
INTEGER(KIND=int4) :: particle_number
CHARACTER(LEN=50), INTENT(IN) :: particle_kind
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(system_data) :: system
REAL(KIND=real8)  :: qpext


! select particle kind
SELECT CASE (particle_kind)
  
  CASE ("thermomechanical_sphere")
    CALL compute_thermomechanical_sphere_external_device_heat_power(particle_number,R,V,system,qpext)

  CASE DEFAULT
    WRITE(error_unit,*) "particle kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_particle_external_device_heat_power
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_particle_mass_and_inertia(particle_number,particle_kind,system)

! dummy arguments
INTEGER(KIND=int4) :: particle_number
CHARACTER(LEN=50), INTENT(IN) :: particle_kind
TYPE(system_data) :: system


! select particle kind
SELECT CASE (particle_kind)
  
  CASE ("sphere1")
    CALL compute_sphere1_mass_and_inertia(particle_number,system)
  
  CASE ("sphere2")
    CALL compute_sphere2_mass_and_inertia(particle_number,system)
  
  CASE ("biotissue_sphere")
    CALL compute_biotissue_sphere_mass_and_inertia(particle_number,system)

  CASE ("lipid_membrane_sphere")
    CALL compute_lipid_membrane_sphere_mass_and_inertia(particle_number,system)

  CASE ("internal_membrane_sphere")
    CALL compute_internal_membrane_sphere_mass_and_inertia(particle_number,system)

  CASE ("external_membrane_sphere")
    CALL compute_external_membrane_sphere_mass_and_inertia(particle_number,system)

  CASE ("jet_sphere")
    CALL compute_jet_sphere_mass_and_inertia(particle_number,system)

  CASE ("granular_sphere")
    CALL compute_granular_sphere_mass_and_inertia(particle_number,system)

  CASE ("granular_sphere_in_given_fluid_flow")
    CALL compute_granular_sphere_igff_mass_and_inertia(particle_number,system)
    
  CASE ("thermomechanical_sphere")
    CALL compute_thermomechanical_sphere_mass_and_inertia(particle_number,system)

  CASE DEFAULT
    WRITE(error_unit,*) "particle kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_particle_mass_and_inertia
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_particle_thermal_properties(particle_number,particle_kind,system)

! dummy arguments
INTEGER(KIND=int4) :: particle_number
CHARACTER(LEN=50), INTENT(IN) :: particle_kind
TYPE(system_data) :: system


! select particle kind
SELECT CASE (particle_kind)
    
  CASE ("thermomechanical_sphere")
    CALL compute_thermomechanical_sphere_thermal_properties(particle_number,system)

  CASE DEFAULT
    WRITE(error_unit,*) "particle kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_particle_thermal_properties
!!---------------------------------------------------------------------------------------------------------------------


END MODULE particles_classes

!!=============================================================================================================================
!!=============================================================================================================================


MODULE adhesion_class

USE particles_classes
USE array_of_vectors_class
USE particle_contact_lists_subroutines
USE wall_contact_lists_subroutines

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere_sphere_adhesion_force_and_moment(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet,khertz,kjkr)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, dnnet, khertz, kjkr
REAL(KIND=real8), DIMENSION(3) :: nij, vi, vj, wi, wj, fadh


! select adhesion model
SELECT CASE (system%adhesion_model%model_type)
  CASE ("normal_spring_dashpot")
    CALL sphere_sphere_normal_spring_dashpot_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet)
  CASE ("classical_jkr")
    CALL sphere_sphere_classical_jkr_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet,khertz,kjkr)
  CASE ("modified_jkr")
    CALL sphere_sphere_modified_jkr_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet)
  !CASE ("normal_and_tangential_spring_dashpot")
    !CALL sphere_sphere_normal_and_tangential_spring_dashpot_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh)
  CASE DEFAULT
    WRITE(error_unit,*) " Adhesion model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE compute_sphere_sphere_adhesion_force_and_moment
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere_wall_adhesion_force_and_moment(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: nout, vi, vw, wi, ww, fadh


! select adhesion model
SELECT CASE (system%adhesion_model%model_type)
  CASE ("normal_spring_dashpot")
    CALL sphere_wall_normal_spring_dashpot_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)
  CASE ("classical_jkr")
    CALL sphere_wall_classical_jkr_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)
  CASE ("modified_jkr")
    CALL sphere_wall_modified_jkr_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)
  !CASE ("normal_and_tangential_spring_dashpot")
    !CALL sphere_wall_normal_and_tangential_spring_dashpot_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh)
  CASE DEFAULT
    WRITE(error_unit,*) " Adhesion model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE compute_sphere_wall_adhesion_force_and_moment
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_normal_spring_dashpot_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij, fadh

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, ri, rj, mi, mj, Eeq, req, meq, xsini, xsinj, xsin, xsin_adh, kn_adh, dn_adh, &
                      vin, vjn, vreln, factor, overlap_exponent, epsilon_min, epsilon_max, overlap_min, overlap_max
REAL(KIND=real8), PARAMETER :: one=1.0_real8, half=0.5_real8, fourthirds=1.3333333333333333_real8


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
xsinj = system%material_set(jmat)%contact_damping_ratio

! get adhesion model parameters
factor = system%adhesion_model%model_parameters(1)
overlap_exponent = system%adhesion_model%model_parameters(2)
xsin_adh = system%adhesion_model%model_parameters(3)
epsilon_min = system%adhesion_model%model_parameters(4)
epsilon_max = system%adhesion_model%model_parameters(5)
!THINK: maybe the adhesion damping should be set equal to the normal contact damping

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute min/max overlaps for adhesion force
overlap_min = epsilon_min*(ri+rj)
overlap_max = epsilon_max*(ri+rj)

! compute contact damping ratio for the contacting pair (average)
!xsin = half*(xsini+xsinj)
!xsin_adh = xsin

! compute normal spring stiffness and dashpot constant
kn_adh = factor*overlap_exponent*fourthirds*Eeq*SQRT(req)*(overlap**(overlap_exponent-one))  !based on Hertz contact force
dn_adh = 2.0_real8*xsin_adh*SQRT(kn_adh*meq)

! compute normal relative velocity of the contacting pair
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! clear adhesion force vector
fadh = 0.0_real8

! compute normal adhesion force
IF (overlap>=overlap_min .AND. overlap<=overlap_max) THEN
  fadh = factor*fourthirds*Eeq*SQRT(req)*(overlap**overlap_exponent)*nij + dn_adh*vreln*nij   !kn_adh*overlap*nij + dn_adh*vreln*nij
  IF (DOT_PRODUCT(fadh,nij)<0.0_real8) THEN
    fadh = 0.0_real8
  END IF
END IF

! udpate particles´ total adhesion forces
!system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!system%particle(j)%adhesion_force = system%particle(j)%adhesion_force - fadh

END SUBROUTINE sphere_sphere_normal_spring_dashpot_adhesion
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_normal_spring_dashpot_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout, fadh

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, ri, mi, Eeq, req, meq, xsini, xsinw, xsin, xsin_adh, vin, vwn, vreln, kn_adh, &
                      dn_adh, factor, overlap_exponent, epsilon_min, epsilon_max, overlap_min, overlap_max
REAL(KIND=real8), PARAMETER :: one=1.0_real8, half=0.5_real8, fourthirds=1.3333333333333333_real8                     


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s and wall´s contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsinw = system%wall(k)%contact_damping_ratio

! get adhesion model parameters
factor = system%adhesion_model%model_parameters(1)
overlap_exponent = system%adhesion_model%model_parameters(2)
xsin_adh = system%adhesion_model%model_parameters(3)
epsilon_min = system%adhesion_model%model_parameters(4)
epsilon_max = system%adhesion_model%model_parameters(5)
!THINK: maybe the adhesion damping should be set equal to the normal contact damping

! compute equivalent properties for the contacting pair
Eeq = Ei/(one-ni_i*ni_i)
req = ri
meq = mi

! compute min/max overlaps for adhesion force
overlap_min = epsilon_min*ri
overlap_max = epsilon_max*ri

! compute contact damping ratio for the contacting pair (average)
!xsin = half*(xsini+xsinw)
!xsin_adh = xsin

! compute normal spring stiffness and dashpot constant
kn_adh = factor*overlap_exponent*fourthirds*Eeq*SQRT(req)*(overlap**(overlap_exponent-one))  !based on Hertz contact force
dn_adh = 2.0_real8*xsin_adh*SQRT(kn_adh*meq)

! compute normal relative velocity of the contacting pair
vin = DOT_PRODUCT(vi,nout)
vwn = DOT_PRODUCT(vw,nout)
vreln = vin - vwn

! clear adhesion force vector
fadh = 0.0_real8

! compute normal adhesion force
IF (overlap>=overlap_min .AND. overlap<=overlap_max) THEN
  fadh = factor*fourthirds*Eeq*SQRT(req)*(overlap**overlap_exponent)*nout + dn_adh*vreln*nout   !kn_adh*overlap*nout + dn_adh*vreln*nout
  IF (DOT_PRODUCT(fadh,nout)<0.0_real8) THEN
    fadh = 0.0_real8
  END IF
END IF

! udpate particles´ total adhesion force
!system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!!system%wall(k)%adhesion_force = system%wall(k)%adhesion_force - fadh

END SUBROUTINE sphere_wall_normal_spring_dashpot_adhesion
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_classical_jkr_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet,khertz,kjkr)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij, fadh

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, ri, rj, mi, mj, Eeq, req, meq, vin, vjn, vreln, Wadh, A, knnet, khertz, kjkr
REAL(KIND=real8), PARAMETER :: half=0.5_real8, two=2.0_real8, three=3.0_real8, threequarters=0.75_real8, &
                               onequarter=0.25_real8, pi=4.0_real8*ATAN(1.0_real8)


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff

! get adhesion model parameters
Wadh = system%adhesion_model%model_parameters(1)
A = system%adhesion_model%model_parameters(2)

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
req = ri*rj/(ri+rj)
!meq = mi*mj/(mi+mj)

! compute adhesion stiffness and dashpot constant
!knadh = three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!dnadh = A*knadh   !2.0_real8*xsiadh*SQRT(kadh*meq)

! compute net spring stiffness and net dashpot constant
knnet = two*Eeq*SQRT(req*overlap) - three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
khertz = two*Eeq*SQRT(req*overlap)
kjkr = - three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!ttttttttttttt Thursday 24/05/2018 ttttttttttttttttttttttttt
!vin = DOT_PRODUCT(vi,nij)
!vjn = DOT_PRODUCT(vj,nij)
!vreln = vin - vjn
!IF (vreln>0.0_real8) THEN
!  IF (knnet<0.0_real8) THEN
!    knnet = -three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter)) !two*Eeq*SQRT(req*overlap)
!  END IF
!  !knnet = two*Eeq*SQRT(req*overlap) + three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!ELSE IF (vreln<0.0_real8) THEN
!  IF (knnet<0.0_real8) THEN
!    knnet = -three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter)) !two*Eeq*SQRT(req*overlap)
!  END IF
!  !knnet = two*Eeq*SQRT(req*overlap) + three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!END IF
!ttttttttttttt Thursday 24/05/2018 ttttttttttttttttttttttttt
dnnet = A*knnet

! compute normal relative velocity of the contacting pair
!vin = DOT_PRODUCT(vi,nij)
!vjn = DOT_PRODUCT(vj,nij)
!vreln = vin - vjn

! compute normal adhesion force (includes damping contribution)
!fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nij - dadh*vreln*nij  !+ dadh*vreln*nij
fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nij !- dnet*vreln*nij
!IF (dnet<0.0_real8) THEN   ! (DOT_PRODUCT(fadh,nij)<0.0_real8) THEN
!  fadh = 0.0_real8  ! fadh = elastic part only?
!END IF

! udpate particles´ total adhesion forces
!system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!system%particle(j)%adhesion_force = system%particle(j)%adhesion_force - fadh

END SUBROUTINE sphere_sphere_classical_jkr_adhesion
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_classical_jkr_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout, fadh

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, ri, mi, Eeq, req, meq, vin, vwn, vreln, Wadh, A, knnet
REAL(KIND=real8), PARAMETER :: half=0.5_real8, two=2.0_real8, three=3.0_real8, threequarters=0.75_real8, &
                               onequarter=0.25_real8, pi=4.0_real8*ATAN(1.0_real8)


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff

! get adhesion model parameters
Wadh = system%adhesion_model%model_parameters(1)
A = system%adhesion_model%model_parameters(2)

! compute equivalent properties for the contacting pair
Eeq = Ei/(1.0_real8-ni_i*ni_i)
req = ri
!meq = mi

! compute adhesion stiffness and dashpot constant
!knadh = three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!dnadh = A*knadh   !2.0_real8*xsiadh*SQRT(kadh*meq)

! compute net spring stiffness and net dashpot constant
knnet = two*Eeq*SQRT(req*overlap) - three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
dnnet = A*knnet

! compute normal relative velocity of the contacting pair
!vin = DOT_PRODUCT(vi,nout)
!vwn = DOT_PRODUCT(vw,nout)
!vreln = vin - vwn

! compute normal adhesion force (includes damping contribution)
fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nout !- dnet*vreln*nout  ! + dadh*vreln*nout
!fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nout - dadh*vreln*nout
!IF (DOT_PRODUCT(fadh,nout)<0.0_real8) THEN
!  fadh = 0.0_real8
!END IF

! udpate particles´ total adhesion force
!system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!!system%wall(k)%adhesion_force = system%wall(k)%adhesion_force - fadh

END SUBROUTINE sphere_wall_classical_jkr_adhesion
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_modified_jkr_adhesion(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij, fadh

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, ri, rj, mi, mj, Eeq, req, meq, vin, vjn, vreln, Wadh, A, &
                      xsini, xsinj, xsin, kncon, knadh, knnet
REAL(KIND=real8), PARAMETER :: half=0.5_real8, two=2.0_real8, three=3.0_real8, threequarters=0.75_real8, &
                               onequarter=0.25_real8, pi=4.0_real8*ATAN(1.0_real8)


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
xsinj = system%material_set(jmat)%contact_damping_ratio

! get adhesion model parameters
Wadh = system%adhesion_model%model_parameters(1)
!A = system%adhesion_model%model_parameters(2)

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute contact damping ratio for the contacting pair (average)
xsin = half*(xsini+xsinj)

! compute adhesion stiffness and dashpot constant
!kadh = three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!dadh = A*kadh   !2.0_real8*xsiadh*SQRT(kadh*meq)

! compute net spring stiffness and net dashpot constant
kncon = two*Eeq*SQRT(req*overlap)
knadh = three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!knnet =   !kncon - three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
knnet = kncon - knadh
IF (knnet<0.0_real8) THEN
  knnet = 0.0_real8
END IF
dnnet = 2.0_real8*xsin*SQRT(knnet*meq)   !2.0_real8*xsin*SQRT(kncon*meq) - A*knadh   !2.0_real8*xsin*SQRT(knnet*meq)   !A*knnet

! compute normal relative velocity of the contacting pair
!vin = DOT_PRODUCT(vi,nij)
!vjn = DOT_PRODUCT(vj,nij)
!vreln = vin - vjn

! compute normal adhesion force (includes damping contribution)
fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nij !- dadh*vreln*nij  !+ dadh*vreln*nij
!fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nij - dnet*vreln*nij
!IF (DOT_PRODUCT(fadh,nij)<0.0_real8) THEN
!  fadh = 0.0_real8
!END IF

! udpate particles´ total adhesion forces
!system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!system%particle(j)%adhesion_force = system%particle(j)%adhesion_force - fadh

END SUBROUTINE sphere_sphere_modified_jkr_adhesion
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_modified_jkr_adhesion(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, dnnet
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout, fadh

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, ri, mi, xsini, xsinw, xsin, Eeq, req, meq, vin, vwn, vreln, kncon, knadh, knnet, Wadh
REAL(KIND=real8), PARAMETER :: half=0.5_real8, two=2.0_real8, three=3.0_real8, threequarters=0.75_real8, &
                               onequarter=0.25_real8, pi=4.0_real8*ATAN(1.0_real8)


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s and wall´s material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsinw = system%wall(k)%contact_damping_ratio

! get adhesion model parameters
Wadh = system%adhesion_model%model_parameters(1)

! compute equivalent properties for the contacting pair
Eeq = Ei/(1.0_real8-ni_i*ni_i)
req = ri
meq = mi

! compute contact damping ratio for the contacting pair (average)
xsin = half*(xsini+xsinw)

! compute adhesion stiffness and dashpot constant
!kadh = three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!dadh = A*kadh   !2.0_real8*xsiadh*SQRT(kadh*meq)

! compute net spring stiffness and net dashpot constant
kncon = two*Eeq*SQRT(req*overlap)
knadh = three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!knet = two*Eeq*SQRT(req*overlap) - three*(req**threequarters)*SQRT(half*pi*Wadh*Eeq)*(overlap**(-onequarter))
!dnet = A*knet
knnet = kncon - knadh
IF (knnet<0.0_real8) THEN
  knnet = 0.0_real8
END IF
dnnet = 2.0_real8*xsin*SQRT(knnet*meq)

! compute normal relative velocity of the contacting pair
!vin = DOT_PRODUCT(vi,nout)
!vwn = DOT_PRODUCT(vw,nout)
!vreln = vin - vwn

! compute normal adhesion force
fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nout !- dadh*vreln*nout  ! + dadh*vreln*nout
!fadh = two*(req**threequarters)*SQRT(two*pi*Wadh*Eeq)*(overlap**threequarters)*nout - dnet*vreln*nout
!IF (DOT_PRODUCT(fadh,nout)<0.0_real8) THEN
!  fadh = 0.0_real8
!END IF

! udpate particles´ total adhesion force
!system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!!system%wall(k)%adhesion_force = system%wall(k)%adhesion_force - fadh

END SUBROUTINE sphere_wall_modified_jkr_adhesion
!!---------------------------------------------------------------------------------------------------------------------

! add here subroutines for other adhesion models

END MODULE adhesion_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE flat_rigid_wall

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_flat_rigid_wall_outside_normal(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system

! local variables
REAL(KIND=real8), DIMENSION(3) :: in_nout


! Note: currently, the normal for flat walls does not change over time; changes will be possible once wall spin is coded

! compute outside normal
system%wall(wnum)%outside_normal = normalize(system%wall(wnum)%initial_outside_normal)

END SUBROUTINE compute_flat_rigid_wall_outside_normal
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_flat_rigid_wall_position_and_velocity(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set wall velocity (flat rigid walls are considered with uniform motion)
system%wall(wnum)%velocity = system%wall(wnum)%initial_velocity

! compute wall position (flat rigid walls are considered with uniform motion)
system%wall(wnum)%point_position = system%wall(wnum)%point_position + system%wall(wnum)%velocity*dt 

END SUBROUTINE compute_flat_rigid_wall_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_flat_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set wall heating rate (rate is considered constant currently; this may be changed in the future)
system%wall(wnum)%heating_rate = system%wall(wnum)%initial_heating_rate

! compute wall temperature (flat rigid walls are considered with uniform heating; this may be changed in the future)
system%wall(wnum)%temperature = system%wall(wnum)%temperature + system%wall(wnum)%heating_rate*dt 

END SUBROUTINE compute_flat_rigid_wall_temperature_and_heating_rate
!!---------------------------------------------------------------------------------------------------------------------


END MODULE flat_rigid_wall

!!=============================================================================================================================
!!=============================================================================================================================


MODULE harmonic_flat_rigid_wall

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_harmonic_flat_rigid_wall_outside_normal(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system

! local variables
REAL(KIND=real8), DIMENSION(3) :: in_nout


! Note: currently, the normal for flat walls does not change over time; changes will be possible once wall spin is coded

! compute outside normal
system%wall(wnum)%outside_normal = normalize(system%wall(wnum)%initial_outside_normal)

END SUBROUTINE compute_harmonic_flat_rigid_wall_outside_normal
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_harmonic_flat_rigid_wall_position_and_velocity(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system

! local variables
REAL(KIND=real8), DIMENSION(3) :: a, w


! get wall harmonic parameters
a = system%wall(wnum)%harmonic_amplitudes
w = system%wall(wnum)%harmonic_frequencies

! compute wall position
system%wall(wnum)%point_position = system%wall(wnum)%initial_point_position + a*sin(w*t) 

! compute wall velocity
system%wall(wnum)%velocity = a*w*cos(w*t)

END SUBROUTINE compute_harmonic_flat_rigid_wall_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_harmonic_flat_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set wall heating rate (rate is considered constant currently; this may be changed in the future)
system%wall(wnum)%heating_rate = system%wall(wnum)%initial_heating_rate

! compute wall temperature (harmonic_flat_rigid_walls are considered with uniform heating; this may be changed in the future)
system%wall(wnum)%temperature = system%wall(wnum)%temperature + system%wall(wnum)%heating_rate*dt 

END SUBROUTINE compute_harmonic_flat_rigid_wall_temperature_and_heating_rate
!!---------------------------------------------------------------------------------------------------------------------


END MODULE harmonic_flat_rigid_wall

!!=============================================================================================================================
!!=============================================================================================================================


MODULE forced_flat_rigid_wall

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_forced_flat_rigid_wall_outside_normal(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system

! local variables
REAL(KIND=real8), DIMENSION(3) :: nout_i


! Note: currently, the normal for forced flat walls does not change over time; changes will be possible once wall spin is coded

! compute outside normal
system%wall(wnum)%outside_normal = normalize(system%wall(wnum)%initial_outside_normal)

END SUBROUTINE compute_forced_flat_rigid_wall_outside_normal
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_forced_flat_rigid_wall_position_and_velocity(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system

! local variables
REAL(KIND=real8) :: xsiw, kw, dw, ai_norm, fdrivi_norm, mass, norm_fcont, norm_fdriv, deltax_norm, v_norm, vi_norm
REAL(KIND=real8), DIMENSION(3) :: ai, vi, nout, a, v, fdrivi, fdriv, fcont, x_old, v_old, deltax


! get wall´s data
ai = system%wall(wnum)%initial_acceleration
vi = system%wall(wnum)%initial_velocity
fdrivi = system%wall(wnum)%initial_driving_force
nout = system%wall(wnum)%outside_normal
xsiw = system%wall(wnum)%contact_damping_ratio

! get wall´s position and velocity at the beginning of the substep
x_old = system%wall(wnum)%point_position
v_old = system%wall(wnum)%velocity

! get wall´s current total contact force
fcont = system%wall(wnum)%contact_force

! compute initial acceleration´s norm
ai_norm = SQRT(DOT_PRODUCT(ai,ai))

! compute forces norms
fdrivi_norm = SQRT(DOT_PRODUCT(fdrivi,fdrivi))
norm_fcont = SQRT(DOT_PRODUCT(fcont,fcont))

! compute wall accelerated mass
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
IF (ai_norm==0.0_real8) THEN 
  mass = 0.4_real8 !40.0_real8
ELSE
  mass = fdrivi_norm/ai_norm
END IF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute current driving force
fdriv = fdrivi + fcont

! compute current acceleration
IF (norm_fcont==0.0_real8) THEN
  a = ai
ELSE
  a = fdriv/mass
END IF

! compute current velocity
!v = v_old + a*dt

! add damping contribution when in contact     !!if in decompression stage
!IF (norm_fcont>0.0_real8 .AND. DOT_PRODUCT(v,nout)>0.0_real8) THEN
!IF (DOT_PRODUCT(v,fdrivi)<0.0_real8) THEN
IF (norm_fcont>0.0_real8) THEN
  deltax = x_old - system%wall(wnum)%initial_point_position
  deltax_norm = SQRT(DOT_PRODUCT(deltax,deltax))
  IF (deltax_norm==0.0_real8) kw = 0.0_real8
  IF (deltax_norm/=0.0_real8) kw = SQRT(DOT_PRODUCT(fdrivi,fdrivi))/deltax_norm !<=notice: I´m using fdrivi instead of fdriv on kw!
  dw = 2.0_real8*xsiw*SQRT(mass*kw)
  !IF (norm_fcont>=fdrivi_norm) THEN
    fdriv = fdriv - dw*v_old  !- dw*v
  !ELSE IF (norm_fcont<fdrivi_norm .AND. norm_fcont>0.0_real8) THEN
  !  fdriv = fdriv + dw*v
  !ELSE IF (norm_fcont==0.0_real8) THEN  !<=this never happens because the outer IF is just for fcont>0 => move it to elsewhere
  !  fdriv = fdrivi
  !END IF
  !a = fdriv/mass
  !v = v_old + a*dt
END IF

a = fdriv/mass
v = v_old + a*dt

! check for consistency
! reset wall acceleration after damping contribution in case the initial acceleration is zero (i.e., piston at constant speed)
v_norm = SQRT(DOT_PRODUCT(v,v))
vi_norm = SQRT(DOT_PRODUCT(vi,vi))
IF ( ai_norm==0.0_real8 .AND. (DOT_PRODUCT(v,vi)>0.0_real8 .AND. v_norm>vi_norm) ) THEN 
  v = vi
  a = ai
END IF

! correct wall velocity in case of contact release while rebouncing
IF (DOT_PRODUCT(v,fdrivi)<0.0_real8 .AND. norm_fcont==0.0_real8) THEN
  v = system%wall(wnum)%initial_velocity !0.0_real8
END IF

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!IF (norm_fcont==0.0_real8 .AND. ai_norm==0.0_real8) THEN
!  v = system%wall(wnum)%initial_velocity
!  a = 0.0_real8
!END IF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! update wall velocity
system%wall(wnum)%velocity = v  !v_old + a*dt

! update wall position
system%wall(wnum)%point_position = x_old + v_old*dt + 0.5_real8*a*dt*dt

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! check for consistency (infiltration problem)
IF (system%wall(wnum)%point_position(1)>=0.05_real8) THEN
  system%wall(wnum)%point_position(1) = 0.05_real8
  system%wall(wnum)%velocity = 0.0_real8
END IF  
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END SUBROUTINE compute_forced_flat_rigid_wall_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_forced_flat_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set wall heating rate (rate is considered constant currently; this may be changed in the future)
system%wall(wnum)%heating_rate = system%wall(wnum)%initial_heating_rate

! compute wall temperature (these types of walls are considered with uniform heating; this may be changed in the future)
system%wall(wnum)%temperature = system%wall(wnum)%temperature + system%wall(wnum)%heating_rate*dt 

END SUBROUTINE compute_forced_flat_rigid_wall_temperature_and_heating_rate
!!---------------------------------------------------------------------------------------------------------------------


END MODULE forced_flat_rigid_wall

!!=============================================================================================================================
!!=============================================================================================================================


MODULE user_defined_rigid_wall

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_user_defined_rigid_wall_outside_normal(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system

! local variables
REAL(KIND=real8), DIMENSION(3) :: face1_startpoint, face1_endpoint, face2_startpoint, face2_endpoint, &
                                  face3_startpoint, face3_endpoint, face4_startpoint, face4_endpoint

! Note: currently, the normal for these walls does not change over time; changes will be possible once wall spin is coded

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! A "user-defined wall" is a multifaced wall, each of which with a finite length; the data for face 1 is given in the input file, as if it were
! a common rigid_wall, whereas the data for the remaining faces must be entered here
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!! enter points defining each face (points must belong to the faces)
!face1_point = system%wall(wnum)%point_position
!face2_point(1) = -0.0275_real8
!face2_point(2) = 0.57_real8
!face2_point(3) = 0.0_real8
!face3_point(1) = 0.0275_real8
!face3_point(2) = 0.57_real8
!face3_point(3) = 0.0_real8
!face4_point(1) = 0.0275_real8
!face4_point(2) = 0.57_real8
!face4_point(3) = 0.0_real8

!! enter coordinate limits for each face
!face1_xi = -0.3775_real8
!face1_xf = -0.0275_real8
!face1_yi = 0.57_real8
!face1_yf = 1.01_real8
!face2_yi = 0.35_real8
!face2_yf = 0.57_real8
!face3_yi = 0.35_real8
!face3_yf = 0.57_real8
!face4_xi = 0.0275_real8
!face4_xf = 0.3775_real8
!face4_yi = 0.57_real8
!face4_yf = 1.01_real8

!face1_startpoint(1) = -0.3775_real8
!face1_startpoint(2) = 1.01_real8
!face1_startpoint(3) = 0.0_real8
!face1_endpoint(1) = -0.0275_real8
!face1_endpoint(2) = -0.0275_real8

! enter outside normals for each face
!!system%wall(wnum)%outside_normal = normalize(system%wall(wnum)%initial_outside_normal)
!face1_outside_normal = normalize(system%wall(wnum)%initial_outside_normal)
!face2_outside_normal(1) = -1.0_real8
!face2_outside_normal(2) = 0.0_real8
!face2_outside_normal(3) = 0.0_real8
!face3_outside_normal(1) = 1.0_real8
!face3_outside_normal(2) = 0.0_real8
!face3_outside_normal(3) = 0.0_real8
!face4_outside_normal(1) = 1.0_real8
!face4_outside_normal(2) = -1.0_real8
!face4_outside_normal(3) = 0.0_real8



END SUBROUTINE compute_user_defined_rigid_wall_outside_normal
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_user_defined_rigid_wall_position_and_velocity(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set wall velocity (flat rigid walls are considered with uniform motion)
system%wall(wnum)%velocity = system%wall(wnum)%initial_velocity

! compute wall position (flat rigid walls are considered with uniform motion)
system%wall(wnum)%point_position = system%wall(wnum)%point_position + system%wall(wnum)%velocity*dt 

END SUBROUTINE compute_user_defined_rigid_wall_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_user_defined_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: wnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set wall heating rate (rate is considered constant currently; this may be changed in the future)
system%wall(wnum)%heating_rate = system%wall(wnum)%initial_heating_rate

! compute wall temperature (these types of walls are considered with uniform heating; this may be changed in the future)
system%wall(wnum)%temperature = system%wall(wnum)%temperature + system%wall(wnum)%heating_rate*dt 

END SUBROUTINE compute_user_defined_rigid_wall_temperature_and_heating_rate
!!---------------------------------------------------------------------------------------------------------------------


END MODULE user_defined_rigid_wall

!!=============================================================================================================================
!!=============================================================================================================================


MODULE rigid_walls_classes

USE flat_rigid_wall
USE harmonic_flat_rigid_wall
USE forced_flat_rigid_wall
USE user_defined_rigid_wall
!USE cylindrical_rigid_wall
!etc

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_rigid_wall_position_and_velocity(wnum,wkind,t,dt,system)

! dummy arguments
INTEGER(KIND=int4)            :: wnum
CHARACTER(LEN=50), INTENT(IN) :: wkind
REAL(KIND=real8)              :: t, dt
TYPE(system_data)             :: system


! select wall kind
SELECT CASE (wkind)

  CASE ("flat_rigid_wall")
    CALL compute_flat_rigid_wall_position_and_velocity(wnum,t,dt,system)
      
  CASE ("harmonic_flat_rigid_wall")
    CALL compute_harmonic_flat_rigid_wall_position_and_velocity(wnum,t,dt,system)

  CASE ("forced_flat_rigid_wall")
    CALL compute_forced_flat_rigid_wall_position_and_velocity(wnum,t,dt,system)

  CASE ("user_defined_rigid_wall")
    CALL compute_user_defined_rigid_wall_position_and_velocity(wnum,t,dt,system)

  CASE DEFAULT
    WRITE(error_unit,*) " wall kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_rigid_wall_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_rigid_wall_outside_normal(wnum,wkind,t,dt,system)

! dummy arguments
INTEGER(KIND=int4)            :: wnum
CHARACTER(LEN=50), INTENT(IN) :: wkind
REAL(KIND=real8)              :: t, dt
TYPE(system_data)             :: system


! select wall kind
SELECT CASE (wkind)

  CASE ("flat_rigid_wall")
    CALL compute_flat_rigid_wall_outside_normal(wnum,t,dt,system)
      
  CASE ("harmonic_flat_rigid_wall")
    CALL compute_harmonic_flat_rigid_wall_outside_normal(wnum,t,dt,system)

  CASE ("forced_flat_rigid_wall")
    CALL compute_forced_flat_rigid_wall_outside_normal(wnum,t,dt,system)

  CASE ("user_defined_rigid_wall")
    CALL compute_user_defined_rigid_wall_outside_normal(wnum,t,dt,system)

  CASE DEFAULT
    WRITE(error_unit,*) " wall kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_rigid_wall_outside_normal
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_rigid_wall_temperature_and_heating_rate(wnum,wkind,t,dt,system)

! dummy arguments
INTEGER(KIND=int4)            :: wnum
CHARACTER(LEN=50), INTENT(IN) :: wkind
REAL(KIND=real8)              :: t, dt
TYPE(system_data)             :: system


! select wall kind
SELECT CASE (wkind)

  CASE ("flat_rigid_wall")
    CALL compute_flat_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)
      
  CASE ("harmonic_flat_rigid_wall")
    CALL compute_harmonic_flat_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

  CASE ("forced_flat_rigid_wall")
    CALL compute_forced_flat_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

  CASE ("user_defined_rigid_wall")
    CALL compute_user_defined_rigid_wall_temperature_and_heating_rate(wnum,t,dt,system)

  CASE DEFAULT
    WRITE(error_unit,*) " wall kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_rigid_wall_temperature_and_heating_rate
!!---------------------------------------------------------------------------------------------------------------------


END MODULE rigid_walls_classes

!!=============================================================================================================================
!!=============================================================================================================================


MODULE external_heating_devices_class

USE system_data_types

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_external_heating_device_position_and_velocity(dnum,dkind,t,dt,system)

! dummy arguments
INTEGER(KIND=int4)            :: dnum
CHARACTER(LEN=50), INTENT(IN) :: dkind
REAL(KIND=real8)              :: t, dt
TYPE(system_data)             :: system


! select device kind
SELECT CASE (dkind)

  CASE ("fire_nozzle")
    CALL compute_heating_device_position_and_velocity(dnum,t,dt,system)
      
  CASE ("laser_beam")
    CALL compute_heating_device_position_and_velocity(dnum,t,dt,system)
    
  !CASE ("harmonic_fire_nozzle")
  !  CALL compute_harmonic_heating_device_position_and_velocity(dnum,t,dt,system)

  !CASE ("harmonic_laser_beam")
  !  CALL compute_harmonic_heating_device_position_and_velocity(dnum,t,dt,system)
    
  !CASE ("user_defined_heating_device")
  !  CALL compute_user_defined_heating_device_position_and_velocity(dnum,t,dt,system)

  CASE DEFAULT
    WRITE(error_unit,*) " external heating device kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_external_heating_device_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_external_heating_device_intensity(dnum,dkind,t,dt,system)

! dummy arguments
INTEGER(KIND=int4)            :: dnum
CHARACTER(LEN=50), INTENT(IN) :: dkind
REAL(KIND=real8)              :: t, dt
TYPE(system_data)             :: system


! select device kind
SELECT CASE (dkind)

  CASE ("fire_nozzle")
    CALL compute_heating_device_intensity(dnum,t,dt,system)
      
  CASE ("laser_beam")
    CALL compute_heating_device_intensity(dnum,t,dt,system)

  !CASE ("harmonic_fire_nozzle")
  !  CALL compute_harmonic_heating_device_intensity(dnum,t,dt,system)

    !CASE ("harmonic_laser_beam")
  !  CALL compute_harmonic_heating_device_intensity(dnum,t,dt,system)

  !CASE ("user_defined_heating_device")
  !  CALL compute_user_defined_heating_device_intensity(dnum,t,dt,system)

  CASE DEFAULT
    WRITE(error_unit,*) " external heating device kind not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE compute_external_heating_device_intensity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_heating_device_position_and_velocity(dnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: dnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set device velocity (these kinds of heating devices are considered with uniform motion)
system%external_heating_device(dnum)%velocity = system%external_heating_device(dnum)%initial_velocity

! compute device position (these kinds of heating devices are considered with uniform motion)
system%external_heating_device(dnum)%striking_position = system%external_heating_device(dnum)%striking_position &
                                                       + system%external_heating_device(dnum)%velocity*dt 
END SUBROUTINE compute_heating_device_position_and_velocity
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_heating_device_intensity(dnum,t,dt,system)

! dummy arguments
INTEGER(KIND=int4) :: dnum
REAL(KIND=real8)   :: t, dt
TYPE(system_data)  :: system


! set device intensity increase rate (rate is considered constant currently; this may be changed in the future)
system%external_heating_device(dnum)%intensity_increase_rate = system%external_heating_device(dnum)%initial_intensity_increase_rate

! compute device intensity (these types of devices are considered with uniformly increasing intensity; this may be changed in the future)
system%external_heating_device(dnum)%intensity = system%external_heating_device(dnum)%intensity &
                                               + system%external_heating_device(dnum)%intensity_increase_rate*dt 

END SUBROUTINE compute_heating_device_intensity
!!---------------------------------------------------------------------------------------------------------------------


END MODULE external_heating_devices_class

!!=============================================================================================================================
!!=============================================================================================================================
    
    
MODULE periodic_boundary_conditions

USE particles_classes
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE enforce_periodic_boundary_conditions(R,system)

! dummy arguments
TYPE(system_data) :: system
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R

! local variables
INTEGER(KIND=int4) :: i
REAL(KIND=real8) :: xmin, xmax, ymin, ymax, zmin, zmax, xi, yi, zi


! get domain limits for periodic BCs
xmin=system%general_cv%perbc_xmin
xmax=system%general_cv%perbc_xmax
ymin=system%general_cv%perbc_ymin
ymax=system%general_cv%perbc_ymax
zmin=system%general_cv%perbc_zmin
zmax=system%general_cv%perbc_zmax

! enforce periodic BCs on particles
DO i=1,system%no_particles
  xi = system%particle(i)%position(1)
  yi = system%particle(i)%position(2)
  zi = system%particle(i)%position(3)
  IF (xi>xmax) system%particle(i)%position(1) = xi - (xmax - xmin)
  IF (xi<xmin) system%particle(i)%position(1) = xi + (xmax - xmin)
  IF (yi>ymax) system%particle(i)%position(2) = yi - (ymax - ymin)
  IF (yi<ymin) system%particle(i)%position(2) = yi + (ymax - ymin)
  IF (zi>zmax) system%particle(i)%position(3) = zi - (zmax - zmin)
  IF (zi<zmin) system%particle(i)%position(3) = zi + (zmax - zmin)
  R(i)%r = system%particle(i)%position
END DO

END SUBROUTINE enforce_periodic_boundary_conditions
!!---------------------------------------------------------------------------------------------------------------------


END MODULE periodic_boundary_conditions

!!=============================================================================================================================
!!=============================================================================================================================


MODULE grid_subroutines

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sort_particles_into_grid_cells(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, ndivx, ndivy, ndivz, ix, iy, iz, cellx, celly, cellz, cell_no_particles
REAL(KIND=real8)   :: lx, ly, lz, xbegin, ybegin, zbegin, xi, yi, zi
INTEGER(KIND=int4), DIMENSION(:), ALLOCATABLE :: aux_list


! get grid data
lx = system%grid%xend - system%grid%xbeg
ly = system%grid%yend - system%grid%ybeg
lz = system%grid%zend - system%grid%zbeg
ndivx = system%grid%ndivx
ndivy = system%grid%ndivy
ndivz = system%grid%ndivz

! check for data inconsistency
IF (lx<0.0_real8 .OR. (ly<0.0_real8 .OR. lz<0.0_real8)) THEN
  WRITE(error_unit,*) " CAUTION: grid coodinates yield negative grid length. Please correct"
  WRITE(default_output_unit,*)" CAUTION: grid coodinates yield negative grid length. Please correct"
END IF

! allocate grid cells
ALLOCATE(system%grid%cell(0:(ndivx+1),0:(ndivy+1),0:(ndivz+1)))

! initialize cell lists
DO ix=0,ndivx+1
  DO iy=0,ndivy+1
    DO iz=0,ndivz+1
      system%grid%cell(ix,iy,iz)%no_particles_into_cell = 0
      ALLOCATE(system%grid%cell(ix,iy,iz)%particle_list(1:1))
      system%grid%cell(ix,iy,iz)%particle_list(1) = 0
    END DO
  END DO
END DO

! get grid initial boundary coordinates
xbegin = system%grid%xbeg
ybegin = system%grid%ybeg
zbegin = system%grid%zbeg

! initialize number of particles into each cell
cell_no_particles = 0

! sort particles
DO i=1,system%no_particles
  xi=system%particle(i)%position(1)
  yi=system%particle(i)%position(2)
  zi=system%particle(i)%position(3)
  cellx = INT(AINT(((xi-xbegin)/lx)*ndivx))+1
  celly = INT(AINT(((yi-ybegin)/ly)*ndivy))+1
  cellz = INT(AINT(((zi-zbegin)/lz)*ndivz))+1
  IF (cellx<=0) cellx=1
  IF (cellx>=(ndivx+1)) cellx=ndivx
  IF (celly<=0) celly=1
  IF (celly>=(ndivy+1)) celly=ndivy
  IF (cellz<=0) cellz=1
  IF (cellz>=(ndivz+1)) cellz=ndivz
  system%particle(i)%cell_address(1) = cellx
  system%particle(i)%cell_address(2) = celly
  system%particle(i)%cell_address(3) = cellz
  system%grid%cell(cellx,celly,cellz)%no_particles_into_cell = system%grid%cell(cellx,celly,cellz)%no_particles_into_cell + 1
  cell_no_particles = system%grid%cell(cellx,celly,cellz)%no_particles_into_cell
  ALLOCATE(aux_list(1:cell_no_particles))
  aux_list(1:(cell_no_particles-1)) = system%grid%cell(cellx,celly,cellz)%particle_list
  aux_list(cell_no_particles) = i
  CALL MOVE_ALLOC(aux_list,system%grid%cell(cellx,celly,cellz)%particle_list)
END DO

END SUBROUTINE sort_particles_into_grid_cells
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE update_grid_cells_lists(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, ndivx, ndivy, ndivz, ix, iy, iz, cellx, celly, cellz, cell_no_particles
REAL(KIND=real8)   :: lx, ly, lz, xbegin, ybegin, zbegin, xi, yi, zi
INTEGER(KIND=int4), DIMENSION(:), ALLOCATABLE :: aux_list


! get grid data
lx = system%grid%xend - system%grid%xbeg
ly = system%grid%yend - system%grid%ybeg
lz = system%grid%zend - system%grid%zbeg
ndivx = system%grid%ndivx
ndivy = system%grid%ndivy
ndivz = system%grid%ndivz

! deallocate old lists
DO ix=0,ndivx+1
  DO iy=0,ndivy+1
    DO iz=0,ndivz+1
      DEALLOCATE(system%grid%cell(ix,iy,iz)%particle_list)
    END DO
  END DO
END DO

! initialize new lists
DO ix=0,ndivx+1
  DO iy=0,ndivy+1
    DO iz=0,ndivz+1
      system%grid%cell(ix,iy,iz)%no_particles_into_cell = 0
      ALLOCATE(system%grid%cell(ix,iy,iz)%particle_list(1:1))
      system%grid%cell(ix,iy,iz)%particle_list(1) = 0
    END DO
  END DO
END DO

! get grid initial boundary coordinates
xbegin = system%grid%xbeg
ybegin = system%grid%ybeg
zbegin = system%grid%zbeg

! initialize number of particles into each cell
cell_no_particles = 0

! sort particles
DO i=1,system%no_particles
  xi=system%particle(i)%position(1)
  yi=system%particle(i)%position(2)
  zi=system%particle(i)%position(3)
  cellx = INT(AINT(((xi-xbegin)/lx)*ndivx))+1
  celly = INT(AINT(((yi-ybegin)/ly)*ndivy))+1
  cellz = INT(AINT(((zi-zbegin)/lz)*ndivz))+1
  IF (cellx<=0) cellx=1
  IF (cellx>=(ndivx+1)) cellx=ndivx
  IF (celly<=0) celly=1
  IF (celly>=(ndivy+1)) celly=ndivy
  IF (cellz<=0) cellz=1
  IF (cellz>=(ndivz+1)) cellz=ndivz
  system%particle(i)%cell_address(1) = cellx
  system%particle(i)%cell_address(2) = celly
  system%particle(i)%cell_address(3) = cellz
  system%grid%cell(cellx,celly,cellz)%no_particles_into_cell = system%grid%cell(cellx,celly,cellz)%no_particles_into_cell + 1
  cell_no_particles = system%grid%cell(cellx,celly,cellz)%no_particles_into_cell
  ALLOCATE(aux_list(1:cell_no_particles))
  aux_list(1:(cell_no_particles-1)) = system%grid%cell(cellx,celly,cellz)%particle_list
  aux_list(cell_no_particles) = i
  CALL MOVE_ALLOC(aux_list,system%grid%cell(cellx,celly,cellz)%particle_list)
END DO

END SUBROUTINE update_grid_cells_lists
!!---------------------------------------------------------------------------------------------------------------------


END MODULE grid_subroutines

!!=============================================================================================================================
!!=============================================================================================================================


MODULE verlet_lists_subroutines

USE particles_classes
USE grid_subroutines

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE build_particles_verlet_lists(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: ndivx, ndivy, ndivz
INTEGER(KIND=int4) :: i, j, ix, iy, iz, cell_no_particles, other_cell_no_particles, pi, pj, iix, iiy, iiz, i_list_no_particles
REAL(KIND=real8)   :: distance, neighbor_search_distance
REAL(KIND=real8), DIMENSION(3) :: xi_minus_xj
INTEGER(KIND=int4), DIMENSION(:), ALLOCATABLE :: aux_list


! get grid data
ndivx = system%grid%ndivx
ndivy = system%grid%ndivy
ndivz = system%grid%ndivz

! initialize particles´ lists
DO i=1,system%no_particles
  system%particle(i)%no_particles_into_verlet_list = 0
  ALLOCATE(system%particle(i)%verlet_list(1:1))
  system%particle(i)%verlet_list(1) = 0
END DO

! initialize number of particles of list of particle i
i_list_no_particles = 0

! build particles´ verlet lists (loop over grid cells)
DO ix=1,ndivx
  DO iy=1,ndivy
    DO iz=1,ndivz
      cell_no_particles = system%grid%cell(ix,iy,iz)%no_particles_into_cell
      DO pi=1,cell_no_particles
        i=system%grid%cell(ix,iy,iz)%particle_list(pi)
        DO iix=ix-1,ix+1
          DO iiy=iy-1,iy+1
            DO iiz=iz-1,iz+1
              other_cell_no_particles = system%grid%cell(iix,iiy,iiz)%no_particles_into_cell
              DO pj=1,other_cell_no_particles
                j=system%grid%cell(iix,iiy,iiz)%particle_list(pj)
                xi_minus_xj = system%particle(i)%position - system%particle(j)%position
                distance = SQRT(DOT_PRODUCT(xi_minus_xj,xi_minus_xj))
                neighbor_search_distance = system%particle(i)%radius + system%general_cv%verlet_distance + system%particle(j)%radius
                IF (distance<=neighbor_search_distance .AND. i/=j) THEN
                  system%particle(i)%no_particles_into_verlet_list = system%particle(i)%no_particles_into_verlet_list + 1
                  i_list_no_particles = system%particle(i)%no_particles_into_verlet_list
                  ALLOCATE(aux_list(1:i_list_no_particles)) 
                  aux_list(1:(i_list_no_particles-1)) = system%particle(i)%verlet_list
                  aux_list(i_list_no_particles) = j
                  CALL MOVE_ALLOC(aux_list,system%particle(i)%verlet_list)
                END IF
              END DO
            END DO  
          END DO    
        END DO                
      END DO
    END DO
  END DO
END DO

END SUBROUTINE build_particles_verlet_lists
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE update_particles_verlet_lists(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: ndivx, ndivy, ndivz
INTEGER(KIND=int4) :: i, j, ix, iy, iz, cell_no_particles, other_cell_no_particles, pi, pj, iix, iiy, iiz, i_list_no_particles
REAL(KIND=real8)   :: distance, neighbor_search_distance
REAL(KIND=real8), DIMENSION(3) :: xi_minus_xj
INTEGER(KIND=int4), DIMENSION(:), ALLOCATABLE :: aux_list


! get grid data
ndivx = system%grid%ndivx
ndivy = system%grid%ndivy
ndivz = system%grid%ndivz

! deallocate old lists
DO i=1,system%no_particles
  DEALLOCATE(system%particle(i)%verlet_list)
END DO

! initialize new lists
DO i=1,system%no_particles
  system%particle(i)%no_particles_into_verlet_list = 0
  ALLOCATE(system%particle(i)%verlet_list(1:1))
  system%particle(i)%verlet_list(1) = 0
END DO

! initialize number of particles of new list of particle i
i_list_no_particles = 0

! build new verlet lists (loop over grid cells)
DO ix=1,ndivx
  DO iy=1,ndivy
    DO iz=1,ndivz
      cell_no_particles = system%grid%cell(ix,iy,iz)%no_particles_into_cell
      DO pi=1,cell_no_particles
        i=system%grid%cell(ix,iy,iz)%particle_list(pi)
        DO iix=ix-1,ix+1
          DO iiy=iy-1,iy+1
            DO iiz=iz-1,iz+1
              other_cell_no_particles = system%grid%cell(iix,iiy,iiz)%no_particles_into_cell
              DO pj=1,other_cell_no_particles
                j=system%grid%cell(iix,iiy,iiz)%particle_list(pj)
                xi_minus_xj = system%particle(i)%position - system%particle(j)%position
                distance = SQRT(DOT_PRODUCT(xi_minus_xj,xi_minus_xj))
                neighbor_search_distance = system%particle(i)%radius + system%general_cv%verlet_distance + system%particle(j)%radius
                IF (distance<=neighbor_search_distance .AND. i/=j) THEN
                  system%particle(i)%no_particles_into_verlet_list = system%particle(i)%no_particles_into_verlet_list + 1
                  i_list_no_particles = system%particle(i)%no_particles_into_verlet_list
                  ALLOCATE(aux_list(1:i_list_no_particles)) 
                  aux_list(1:(i_list_no_particles-1)) = system%particle(i)%verlet_list
                  aux_list(i_list_no_particles) = j
                  CALL MOVE_ALLOC(aux_list,system%particle(i)%verlet_list)
                END IF
              END DO
            END DO  
          END DO    
        END DO                
      END DO
    END DO
  END DO
END DO

END SUBROUTINE update_particles_verlet_lists
!!---------------------------------------------------------------------------------------------------------------------


END MODULE verlet_lists_subroutines

!!=============================================================================================================================
!!=============================================================================================================================


MODULE particle_contact_lists_subroutines

USE particles_classes
USE contact_history_array_class
USE io_files

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE initialize_particle_contact_lists(i,system)

! dummy arguments
INTEGER(KIND=int4) :: i
TYPE(system_data)  :: system


! allocate contact lists of particle i
system%particle(i)%no_contacting_particles = 0
system%particle(i)%no_contacting_walls = 0
ALLOCATE(system%particle(i)%contacting_particles(1:1))
ALLOCATE(system%particle(i)%contacting_walls(1:1))
system%particle(i)%contacting_particles(1)%part_number = 0
system%particle(i)%contacting_walls(1)%wall_number = 0

! initialize normal springs elongations of particle i
system%particle(i)%contacting_particles(1)%normalspring_elongation = 0.0_real8
system%particle(i)%contacting_walls(1)%normalspring_elongation = 0.0_real8

! initialize normal springs forces of particle i
system%particle(i)%contacting_particles(1)%normalspring_force = 0.0_real8
system%particle(i)%contacting_walls(1)%normalspring_force = 0.0_real8

! initialize friction springs elongations of particle i
system%particle(i)%contacting_particles(1)%fricspring_elongation = 0.0_real8
system%particle(i)%contacting_walls(1)%fricspring_elongation = 0.0_real8

! initialize friction springs forces of particle i
system%particle(i)%contacting_particles(1)%fricspring_force = 0.0_real8
system%particle(i)%contacting_walls(1)%fricspring_force = 0.0_real8

! initialize rolling resistance springs elongations of particle i
system%particle(i)%contacting_particles(1)%rolresspring_elongation = 0.0_real8
system%particle(i)%contacting_walls(1)%rolresspring_elongation = 0.0_real8

END SUBROUTINE initialize_particle_contact_lists
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE add_particle_to_particle_contact_list(i,j,system)

! dummy arguments
INTEGER(KIND=int4) :: i, j
TYPE(system_data)  :: system

! local variables
INTEGER(KIND=int4) :: ncpi, ncpj
TYPE(contact_particle_pair_data), DIMENSION(:), ALLOCATABLE :: iaux_list, jaux_list


! add particle j to contact list of particle i
system%particle(i)%no_contacting_particles = system%particle(i)%no_contacting_particles + 1  ! CAUTION: this must not accumulate within the iterations!
ncpi = system%particle(i)%no_contacting_particles
ALLOCATE(iaux_list(1:ncpi))
iaux_list(1:(ncpi-1)) = system%particle(i)%contacting_particles
iaux_list(ncpi)%part_number = j
iaux_list(ncpi)%normalspring_elongation = 0.0_real8
iaux_list(ncpi)%normalspring_force = 0.0_real8
iaux_list(ncpi)%fricspring_elongation = 0.0_real8
iaux_list(ncpi)%fricspring_force = 0.0_real8
iaux_list(ncpi)%rolresspring_elongation = 0.0_real8
CALL MOVE_ALLOC(iaux_list,system%particle(i)%contacting_particles)

! add particle i to contact list of particle j
system%particle(j)%no_contacting_particles = system%particle(j)%no_contacting_particles + 1  ! CAUTION: this must not accumulate within the iterations!
ncpj = system%particle(j)%no_contacting_particles
ALLOCATE(jaux_list(1:ncpj))
jaux_list(1:(ncpj-1)) = system%particle(j)%contacting_particles
jaux_list(ncpj)%part_number = i
jaux_list(ncpj)%normalspring_elongation = 0.0_real8
jaux_list(ncpj)%normalspring_force = 0.0_real8
jaux_list(ncpj)%fricspring_elongation = 0.0_real8
jaux_list(ncpj)%fricspring_force = 0.0_real8
jaux_list(ncpj)%rolresspring_elongation = 0.0_real8
CALL MOVE_ALLOC(jaux_list,system%particle(j)%contacting_particles)

END SUBROUTINE add_particle_to_particle_contact_list
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE add_wall_to_particle_contact_list(i,w,system)

! dummy arguments
INTEGER(KIND=int4) :: i, w
TYPE(system_data)  :: system

! local variables
INTEGER(KIND=int4) :: no_contwalls
TYPE(contact_wall_pair_data), DIMENSION(:), ALLOCATABLE :: aux_list


! add wall w to contact list of particle i
system%particle(i)%no_contacting_walls = system%particle(i)%no_contacting_walls + 1  ! CAUTION: this must not accumulate within the iterations!
no_contwalls = system%particle(i)%no_contacting_walls
ALLOCATE(aux_list(1:no_contwalls))
aux_list(1:(no_contwalls-1)) = system%particle(i)%contacting_walls
aux_list(no_contwalls)%wall_number = w
aux_list(no_contwalls)%normalspring_elongation = 0.0_real8
aux_list(no_contwalls)%normalspring_force = 0.0_real8
aux_list(no_contwalls)%fricspring_elongation = 0.0_real8
aux_list(no_contwalls)%fricspring_force = 0.0_real8
aux_list(no_contwalls)%rolresspring_elongation = 0.0_real8
CALL MOVE_ALLOC(aux_list,system%particle(i)%contacting_walls)

END SUBROUTINE add_wall_to_particle_contact_list
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE remove_particle_from_particle_contact_list(i,j,system)

! dummy arguments
INTEGER(KIND=int4) :: i, j
TYPE(system_data)  :: system

! local variables
INTEGER(KIND=int4) :: cp, jaddr, iaddr, ncpi, ncpj
TYPE(contact_particle_pair_data), DIMENSION(:), ALLOCATABLE :: iaux_list, jaux_list


! get particle j´s address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! remove particle j from contact list of particle i
system%particle(i)%no_contacting_particles = system%particle(i)%no_contacting_particles - 1
ncpi = system%particle(i)%no_contacting_particles
ALLOCATE(iaux_list(1:ncpi))
iaux_list(1:(jaddr-1)) = system%particle(i)%contacting_particles(1:(jaddr-1))
iaux_list(jaddr:ncpi) = system%particle(i)%contacting_particles((jaddr+1):(ncpi+1))
CALL MOVE_ALLOC(iaux_list,system%particle(i)%contacting_particles)

! get particle i´s address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! remove particle i from contact list of particle j
system%particle(j)%no_contacting_particles = system%particle(j)%no_contacting_particles - 1
ncpj = system%particle(j)%no_contacting_particles
ALLOCATE(jaux_list(1:ncpj))
jaux_list(1:(iaddr-1)) = system%particle(j)%contacting_particles(1:(iaddr-1))
jaux_list(iaddr:ncpj) = system%particle(j)%contacting_particles((iaddr+1):(ncpj+1))
CALL MOVE_ALLOC(jaux_list,system%particle(j)%contacting_particles)

END SUBROUTINE remove_particle_from_particle_contact_list
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE remove_wall_from_particle_contact_list(i,w,system)

! dummy arguments
INTEGER(KIND=int4) :: i, w
TYPE(system_data)  :: system

! local variables
INTEGER(KIND=int4) :: cw, waddr, ncw
TYPE(contact_wall_pair_data), DIMENSION(:), ALLOCATABLE :: aux_list


! get wall address at contact list of particle i
DO cw=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(cw)%wall_number==w) waddr = cw
END DO

! remove wall w from contact list of particle i
system%particle(i)%no_contacting_walls = system%particle(i)%no_contacting_walls - 1
ncw = system%particle(i)%no_contacting_walls
ALLOCATE(aux_list(1:ncw))
aux_list(1:(waddr-1)) = system%particle(i)%contacting_walls(1:(waddr-1))
aux_list(waddr:ncw) = system%particle(i)%contacting_walls((waddr+1):(ncw+1))
CALL MOVE_ALLOC(aux_list,system%particle(i)%contacting_walls)

END SUBROUTINE remove_wall_from_particle_contact_list
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE reset_particle_contact_lists(i,PPCHt,PWCHt,system)

! dummy arguments
INTEGER(KIND=int4) :: i
TYPE(system_data)  :: system
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHt, PWCHt  ! INTENT(IN)

! local variables
INTEGER(KIND=int4) :: j


! deallocate old contact lists
IF (ALLOCATED(system%particle(i)%contacting_particles)) DEALLOCATE(system%particle(i)%contacting_particles)
IF (ALLOCATED(system%particle(i)%contacting_walls)) DEALLOCATE(system%particle(i)%contacting_walls)

! allocate new contact lists
system%particle(i)%no_contacting_particles = SIZE(PPCHt(i)%pair)
system%particle(i)%no_contacting_walls = SIZE(PWCHt(i)%pair)
ALLOCATE(system%particle(i)%contacting_particles(1:SIZE(PPCHt(i)%pair)))
ALLOCATE(system%particle(i)%contacting_walls(1:SIZE(PWCHt(i)%pair)))

! reset particle i´s contact lists (set as the lists at time t)
DO j=1,system%particle(i)%no_contacting_particles
  system%particle(i)%contacting_particles(j)%part_number = PPCHt(i)%pair(j)%number
  system%particle(i)%contacting_particles(j)%normalspring_elongation = PPCHt(i)%pair(j)%normalspring_elongation
  system%particle(i)%contacting_particles(j)%normalspring_force = PPCHt(i)%pair(j)%normalspring_force
  system%particle(i)%contacting_particles(j)%fricspring_elongation = PPCHt(i)%pair(j)%fricspring_elongation
  system%particle(i)%contacting_particles(j)%fricspring_force = PPCHt(i)%pair(j)%fricspring_force
  system%particle(i)%contacting_particles(j)%rolresspring_elongation = PPCHt(i)%pair(j)%rolresspring_elongation
END DO
DO j=1,system%particle(i)%no_contacting_walls
  system%particle(i)%contacting_walls(j)%wall_number = PWCHt(i)%pair(j)%number
  system%particle(i)%contacting_walls(j)%normalspring_elongation = PWCHt(i)%pair(j)%normalspring_elongation
  system%particle(i)%contacting_walls(j)%normalspring_force = PWCHt(i)%pair(j)%normalspring_force
  system%particle(i)%contacting_walls(j)%fricspring_elongation = PWCHt(i)%pair(j)%fricspring_elongation
  system%particle(i)%contacting_walls(j)%fricspring_force = PWCHt(i)%pair(j)%fricspring_force
  system%particle(i)%contacting_walls(j)%rolresspring_elongation = PWCHt(i)%pair(j)%rolresspring_elongation
END DO

END SUBROUTINE reset_particle_contact_lists
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE update_particle_contact_lists(i,PPCHtpdt,PWCHtpdt,system)

! dummy arguments
INTEGER(KIND=int4) :: i
TYPE(system_data)  :: system
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHtpdt, PWCHtpdt  ! INTENT(IN)

! local variables
INTEGER(KIND=int4) :: j

! update particle i´s contact lists
DO j=1,system%particle(i)%no_contacting_particles
  system%particle(i)%contacting_particles(j)%part_number = PPCHtpdt(i)%pair(j)%number
  system%particle(i)%contacting_particles(j)%normalspring_elongation = PPCHtpdt(i)%pair(j)%normalspring_elongation
  system%particle(i)%contacting_particles(j)%normalspring_force = PPCHtpdt(i)%pair(j)%normalspring_force
  system%particle(i)%contacting_particles(j)%fricspring_elongation = PPCHtpdt(i)%pair(j)%fricspring_elongation
  system%particle(i)%contacting_particles(j)%fricspring_force = PPCHtpdt(i)%pair(j)%fricspring_force
  system%particle(i)%contacting_particles(j)%rolresspring_elongation = PPCHtpdt(i)%pair(j)%rolresspring_elongation
END DO
DO j=1,system%particle(i)%no_contacting_walls
  system%particle(i)%contacting_walls(j)%wall_number = PWCHtpdt(i)%pair(j)%number
  system%particle(i)%contacting_walls(j)%normalspring_elongation = PWCHtpdt(i)%pair(j)%normalspring_elongation
  system%particle(i)%contacting_walls(j)%normalspring_force = PWCHtpdt(i)%pair(j)%normalspring_force
  system%particle(i)%contacting_walls(j)%fricspring_elongation = PWCHtpdt(i)%pair(j)%fricspring_elongation
  system%particle(i)%contacting_walls(j)%fricspring_force = PWCHtpdt(i)%pair(j)%fricspring_force
  system%particle(i)%contacting_walls(j)%rolresspring_elongation = PWCHtpdt(i)%pair(j)%rolresspring_elongation
END DO

END SUBROUTINE update_particle_contact_lists
!!---------------------------------------------------------------------------------------------------------------------


END MODULE particle_contact_lists_subroutines
    
!!=============================================================================================================================
!!=============================================================================================================================
    
    
MODULE wall_contact_lists_subroutines

USE particles_classes
USE contact_history_array_class
USE io_files

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE initialize_wall_contact_lists(k,system)

! dummy arguments
INTEGER(KIND=int4) :: k
TYPE(system_data)  :: system


! allocate contact lists of wall k
system%wall(k)%no_contacting_particles = 0
ALLOCATE(system%wall(k)%contacting_particles(1:1))
system%wall(k)%contacting_particles(1)%part_number = 0

! initialize contact and friction forces of wall k
system%wall(k)%contacting_particles(1)%contact_force = 0.0_real8
system%wall(k)%contacting_particles(1)%friction_force = 0.0_real8

END SUBROUTINE initialize_wall_contact_lists
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE add_particle_to_wall_contact_list(k,i,system)

! dummy arguments
INTEGER(KIND=int4) :: k, i
TYPE(system_data)  :: system

! local variables
INTEGER(KIND=int4) :: ncpk
TYPE(wall_contact_list_data), DIMENSION(:), ALLOCATABLE :: kaux_list


! add particle i to contact list of wall k
system%wall(k)%no_contacting_particles = system%wall(k)%no_contacting_particles + 1  ! CAUTION: this must not accumulate within the iterations!
ncpk = system%wall(k)%no_contacting_particles
ALLOCATE(kaux_list(1:ncpk))
kaux_list(1:(ncpk-1)) = system%wall(k)%contacting_particles
kaux_list(ncpk)%part_number = i
kaux_list(ncpk)%contact_force = 0.0_real8
kaux_list(ncpk)%friction_force = 0.0_real8
CALL MOVE_ALLOC(kaux_list,system%wall(k)%contacting_particles)

END SUBROUTINE add_particle_to_wall_contact_list
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE remove_particle_from_wall_contact_list(k,i,system)

! dummy arguments
INTEGER(KIND=int4) :: k, i
TYPE(system_data)  :: system

! local variables
INTEGER(KIND=int4) :: cp, iaddr, ncpk
TYPE(wall_contact_list_data), DIMENSION(:), ALLOCATABLE :: kaux_list


! get particle i´s address at contact list of wall k
DO cp=1,system%wall(k)%no_contacting_particles
  IF (system%wall(k)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! remove particle i from contact list of wall k
system%wall(k)%no_contacting_particles = system%wall(k)%no_contacting_particles - 1
ncpk = system%wall(k)%no_contacting_particles
ALLOCATE(kaux_list(1:ncpk))
kaux_list(1:(iaddr-1)) = system%wall(k)%contacting_particles(1:(iaddr-1))
kaux_list(iaddr:ncpk) = system%wall(k)%contacting_particles((iaddr+1):(ncpk+1))
CALL MOVE_ALLOC(kaux_list,system%wall(k)%contacting_particles)

END SUBROUTINE remove_particle_from_wall_contact_list
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE reset_wall_contact_lists(k,WPCHt,system)

! dummy arguments
INTEGER(KIND=int4) :: k
TYPE(system_data)  :: system
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: WPCHt  ! INTENT(IN)

! local variables
INTEGER(KIND=int4) :: i


! deallocate old contact list
IF (ALLOCATED(system%wall(k)%contacting_particles)) DEALLOCATE(system%wall(k)%contacting_particles)

! allocate new contact list
system%wall(k)%no_contacting_particles = SIZE(WPCHt(k)%pair)
ALLOCATE(system%wall(k)%contacting_particles(1:SIZE(WPCHt(k)%pair)))

! reset wall k´s contact list (set as the list at time t)
DO i=1,system%wall(k)%no_contacting_particles
  system%wall(k)%contacting_particles(i)%part_number = WPCHt(k)%pair(i)%number
  system%wall(k)%contacting_particles(i)%contact_force = WPCHt(k)%pair(i)%contact_force
  system%wall(k)%contacting_particles(i)%friction_force = WPCHt(k)%pair(i)%friction_force

END DO

END SUBROUTINE reset_wall_contact_lists
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE update_wall_contact_lists(k,WPCHtpdt,system)

! dummy arguments
INTEGER(KIND=int4) :: k
TYPE(system_data)  :: system
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: WPCHtpdt  ! INTENT(IN)

! local variables
INTEGER(KIND=int4) :: i

! update wall k´s contact list
DO i=1,system%wall(k)%no_contacting_particles
  system%wall(k)%contacting_particles(i)%part_number = WPCHtpdt(k)%pair(i)%number
  system%Wall(k)%contacting_particles(i)%contact_force = WPCHtpdt(k)%pair(i)%contact_force
  system%Wall(k)%contacting_particles(i)%friction_force = WPCHtpdt(k)%pair(i)%friction_force
END DO

END SUBROUTINE update_wall_contact_lists
!!---------------------------------------------------------------------------------------------------------------------


END MODULE wall_contact_lists_subroutines
    
!!=============================================================================================================================
!!=============================================================================================================================
    
    
MODULE rolling_resistance_class

USE particles_classes
USE array_of_vectors_class
USE particle_contact_lists_subroutines
USE wall_contact_lists_subroutines

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere_sphere_rolling_resistance_moment(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, norm_fcon
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij, mirol, vrelctang


! select rolling resistance model
SELECT CASE (system%solution_cv%rolling_resistance_model)
  CASE ("constant_torque")
    CALL sphere_sphere_constant_torque_rolling_resistance(i,j,wi,wj,norm_fcon,system,mirol)
  CASE ("damper_torque")
    CALL sphere_sphere_damper_torque_rolling_resistance(i,j,vi,vj,wi,wj,nij,norm_fcon,system,mirol)
  CASE ("elastic_plastic_spring_damper_torque")
    CALL sphere_sphere_epsd_torque_rolling_resistance(i,j,wi,wj,overlap,norm_fcon,system,mirol)
  CASE DEFAULT
    WRITE(error_unit,*) " Rolling resistance model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE compute_sphere_sphere_rolling_resistance_moment
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sphere_wall_rolling_resistance_moment(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, norm_fcon
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout, mirol, vrelctang


! select rolling resistance model
SELECT CASE (system%solution_cv%rolling_resistance_model)
  CASE ("constant_torque")
    CALL sphere_wall_constant_torque_rolling_resistance(i,k,wi,ww,norm_fcon,system,mirol)
  CASE ("damper_torque")
    CALL sphere_wall_damper_torque_rolling_resistance(i,k,vi,vw,wi,ww,nout,norm_fcon,system,mirol)
  CASE ("elastic_plastic_spring_damper_torque")
    CALL sphere_wall_epsd_torque_rolling_resistance(i,k,wi,ww,overlap,norm_fcon,system,mirol)
  CASE DEFAULT
    WRITE(error_unit,*) " Rolling resistance model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE compute_sphere_wall_rolling_resistance_moment
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_constant_torque_rolling_resistance(i,j,wi,wj,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: wi, wj
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat, jmat
REAL(KIND=real8)   :: Ri, Rj, Req, miri, mirj, mir, norm_wrel
REAL(KIND=real8), DIMENSION(3) :: wrel


! get particles´ data
Ri = system%particle(i)%radius
Rj = system%particle(j)%radius

! get particles material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirj = system%material_set(jmat)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirj)

! compute equivalent radius for the pair i-j
Req = Ri*Rj/(Ri+Rj)

! compute relative spin
wrel = wi - wj
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*(1.0_real8/norm_wrel)*wrel
END IF

END SUBROUTINE sphere_sphere_constant_torque_rolling_resistance
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_constant_torque_rolling_resistance(i,k,wi,ww,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: wi, ww
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat
REAL(KIND=real8)   :: Ri, Req, miri, mirw, mir, norm_wrel
REAL(KIND=real8), DIMENSION(3) :: wrel


! get particle´s data
Ri = system%particle(i)%radius

! get particle´s material-set number
imat = system%particle(i)%material_set_number

! get rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirw = system%wall(k)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirw)

! compute equivalent radius for the pair i-k
Req = Ri

! compute relative spin
wrel = wi !- ww
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*(1.0_real8/norm_wrel)*wrel
END IF

END SUBROUTINE sphere_wall_constant_torque_rolling_resistance
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_damper_torque_rolling_resistance(i,j,vi,vj,wi,wj,nij,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat, jmat
REAL(KIND=real8)   :: Ri, Rj, Req, miri, mirj, mir, norm_wrel, norm_vrelctang
REAL(KIND=real8), DIMENSION(3) :: wrel, rci, rcj, vrelc, vrelctang


! get particles´ data
Ri = system%particle(i)%radius
Rj = system%particle(j)%radius

! get particles material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirj = system%material_set(jmat)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirj)

! compute equivalent radius for the pair i-j
Req = Ri*Rj/(Ri+Rj)

! compute relative spin
wrel = wi - wj
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! build contact point radial vectors
rci = Ri*nij
rcj = -Rj*nij
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (Ri-overlap)*nij
!rcj = -(Rj-overlap)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute tangential relative velocity at contact point (due to rotation only)
!vrelc = vi + (wi .vector. rci) - vj - (wj .vector. rcj)
vrelc = (wi .vector. rci) - (wj .vector. rcj)
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*norm_vrelctang*(1.0_real8/norm_wrel)*wrel  !((wi .vector. rci) - (wj .vector. rcj))
END IF

END SUBROUTINE sphere_sphere_damper_torque_rolling_resistance
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_damper_torque_rolling_resistance(i,k,vi,vw,wi,ww,nout,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat
REAL(KIND=real8)   :: Ri, Req, miri, mirw, mir, norm_wrel, norm_vrelctang
REAL(KIND=real8), DIMENSION(3) :: wrel, rci, vrelc, vrelctang


! get particle´s data
Ri = system%particle(i)%radius

! get particle´s material-set number
imat = system%particle(i)%material_set_number

! get rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirw = system%wall(k)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirw)

! compute equivalent radius for the pair i-k
Req = Ri

! compute relative spin
wrel = wi !- ww
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! build contact point radial vector
rci = Ri*nout
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (Ri-overlap)*nout
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute tangential relative velocity at contact point (due to rotation only)
!vrelc = vi + (wi .vector. rci) - vw
vrelc = (wi .vector. rci)
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*norm_vrelctang*(1.0_real8/norm_wrel)*wrel  !*(-vrelc .vector. nout)/SQRT(DOT_PRODUCT( (vrelc .vector. nout) , (vrelc .vector. nout) ))    !mir*Req*norm_fcon*vrelctang   !(wi .vector. rci) !-ww*Rw
END IF

END SUBROUTINE sphere_wall_damper_torque_rolling_resistance
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_epsd_torque_rolling_resistance(i,j,wi,wj,overlap,norm_fcon,system,mrol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, norm_fcon
REAL(KIND=real8), DIMENSION(3) :: wi, wj
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mrol

! local variables
INTEGER(KIND=int4) :: imat, jmat, s, cp, jaddr, iaddr
REAL(KIND=real8)   :: Ri, mi, ji, Rj, mj, jj, Req, Ei, Ej, ni_i, ni_j, miri, mirj, mir, xsiri, xsirj, xsir, dt, &
                      Gi, Gj, Geq, norm_wrel, kr, jr, dr, aux_ji, aux_jj, norm_mrol_trial    
REAL(KIND=real8), DIMENSION(3) :: wrel, sij_trial, sij, dxij_trial, xijr, xij_trial, xij, mrol_trial
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, two=2.0_real8


! get particles´ data
Ri = system%particle(i)%radius
mi = system%particle(i)%mass
ji = system%particle(i)%inertia
Rj = system%particle(j)%radius
mj = system%particle(j)%mass
jj = system%particle(j)%inertia

! get particles´ material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
miri = system%material_set(imat)%rolling_resistance_coeff
xsiri = system%material_set(imat)%rolling_resistance_damping_ratio
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
mirj = system%material_set(jmat)%rolling_resistance_coeff
xsirj = system%material_set(jmat)%rolling_resistance_damping_ratio

! get particle j´s address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i´s address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! get time step size
s = system%solution_cv%current_step
dt = system%solution_cv%step(s)%current_dt

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute equivalent properties for the pair i-j
Req = Ri*Rj/(Ri+Rj)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)

! compute rolling_resistance parameters for the contacting pair (average)
mir = 0.5_real8*(miri+mirj)
xsir = 0.5_real8*(xsiri+xsirj)

! compute rolling resistance spring stiffness
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!kr = Req*norm_fcon
!kr = 10000_real8*kr
!kr = mir*Req*norm_fcon/0.05_real8
kr = 8.0_real8*Geq*SQRT(Req*overlap)*Req*Req
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute rolling resistance dashpot constant
aux_ji = ji + mi*Ri*Ri
aux_jj = jj + mj*Rj*Rj
jr = one/((one/aux_ji) + (one/aux_jj))
dr = 2.0_real8*xsir*SQRT(jr*kr)

! compute particles´ relative spin
wrel = wi - wj   ! CAUTION: think more about this!
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute trial incremental elongation of rolling resistance spring
dxij_trial = wrel*dt

! compute trial total elongation of rolling resistance spring
xijr = system%particle(i)%contacting_particles(jaddr)%rolresspring_elongation
!xt_accum_proj = xt_accum - DOT_PRODUCT(xt_accum,nij)*nij
!xt_trial = xt_accum_proj + dxt_trial
xij_trial = xijr + dxij_trial

! compute trial rolling resistance moment
mrol_trial = - kr*xij_trial - dr*wrel

! compute trial rolling resistance moment direction
norm_mrol_trial = SQRT(DOT_PRODUCT(mrol_trial,mrol_trial))
IF (norm_mrol_trial/=0.0_real8) THEN
  sij_trial = (one/norm_mrol_trial)*mrol_trial
ELSE
  sij_trial = 0.0_real8
END IF

! check for stick or slip (verify if trial rolling resistance moment violates the maximum possible value)
IF (norm_mrol_trial<=(mir*Req*norm_fcon)) THEN
  mrol = mrol_trial
  xij = xij_trial
  sij = sij_trial
  !"stick region"
ELSE
  sij = sij_trial
  mrol = mir*Req*norm_fcon*sij
  xij = -(one/kr)*(mrol - dr*wrel)
  !"slip region"
END IF

! update pair´s rolling resistance spring elongation
system%particle(i)%contacting_particles(jaddr)%rolresspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%rolresspring_elongation = xij

END SUBROUTINE sphere_sphere_epsd_torque_rolling_resistance
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_epsd_torque_rolling_resistance(i,k,wi,ww,overlap,norm_fcon,system,mrol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, norm_fcon
REAL(KIND=real8), DIMENSION(3) :: wi, ww
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mrol

! local variables
INTEGER(KIND=int4) :: imat, s, j, waddr
REAL(KIND=real8)   :: Ri, mi, ji, Req, Ei, ni_i, miri, xsiri, mirw, xsirw, Gi, Geq, mir, xsir, &
                      dt, kr, norm_wrel, jr, dr, norm_mrol_trial
REAL(KIND=real8), DIMENSION(3) :: wrel, sij_trial, sij, dxiw_trial, xiwr, xiw_trial, xiw, mrol_trial
REAL(KIND=real8), PARAMETER    :: one=1.0_real8


! get particle´s data
Ri = system%particle(i)%radius
mi = system%particle(i)%mass
ji = system%particle(i)%inertia

! get particle´s material-set number
imat = system%particle(i)%material_set_number

! get particle´s and wall´s material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
miri = system%material_set(imat)%rolling_resistance_coeff
xsiri = system%material_set(imat)%rolling_resistance_damping_ratio
mirw = system%wall(k)%rolling_resistance_coeff
xsirw = system%wall(k)%rolling_resistance_damping_ratio

! get wall´s address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! get time step size
s = system%solution_cv%current_step
dt = system%solution_cv%step(s)%current_dt

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(1.0_real8+ni_i))

! compute equivalent properties for the pair i-k
Req = Ri
Geq = Gi/(2.0_real8-ni_i)

! compute rolling_resistance parameters for the contacting pair (average)
mir = 0.5_real8*(miri+mirw)
xsir = 0.5_real8*(xsiri+xsirw)

! compute rolling resistance spring stiffness
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!kr = Req*norm_fcon
!kr = 10000_real8*kr
!kr = mir*Req*norm_fcon/0.05_real8
kr = 8.0_real8*Geq*SQRT(Req*overlap)*Req*Req
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute rolling resistance dashpot constant
jr = ji + mi*Ri*Ri
dr = 2.0_real8*xsir*SQRT(jr*kr)

! compute particle-wall relative spin
wrel = wi !- ww   ! CAUTION: think more about this!
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute trial incremental elongation of rolling resistance spring
dxiw_trial = wrel*dt

! compute trial total elongation of rolling resistance spring
xiwr = system%particle(i)%contacting_walls(waddr)%rolresspring_elongation
!xt_accum_proj = xt_accum - DOT_PRODUCT(xt_accum,nij)*nij
!xt_trial = xt_accum_proj + dxt_trial
xiw_trial = xiwr + dxiw_trial

! compute trial rolling resistance moment
mrol_trial = - kr*xiw_trial - dr*wrel

! compute trial rolling resistance moment direction
norm_mrol_trial = SQRT(DOT_PRODUCT(mrol_trial,mrol_trial))
IF (norm_mrol_trial/=0.0_real8) THEN
  sij_trial = (one/norm_mrol_trial)*mrol_trial
ELSE
  sij_trial = 0.0_real8
END IF

! check for stick or slip (verify if trial rolling resistance moment violates the maximum possible value)
IF (norm_mrol_trial<=(mir*Req*norm_fcon)) THEN
  mrol = mrol_trial
  xiw = xiw_trial
  sij = sij_trial
  !"stick region"
ELSE
  sij = sij_trial
  mrol = mir*Req*norm_fcon*sij
  xiw = -(one/kr)*(mrol - dr*wrel)
  !"slip region"
END IF

! update pair´s rolling resistance spring elongation
system%particle(i)%contacting_walls(waddr)%rolresspring_elongation = xiw

END SUBROUTINE sphere_wall_epsd_torque_rolling_resistance
!!---------------------------------------------------------------------------------------------------------------------

! add here subroutines for other rolling resistance models

!!---------------------------------------------------------------------------------------------------------------------


! below are the rolling resistance subroutines for the case when the elasticity modulus varies with temperature


!!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE compute_sph_sph_roll_resist_mom_w_thcond(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,Ei,Ej,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, norm_fcon, Ei, Ej
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij, mirol, vrelctang


! select rolling resistance model
SELECT CASE (system%solution_cv%rolling_resistance_model)
  CASE ("constant_torque")
    CALL sph_sph_const_torq_roll_resist_w_thcond(i,j,wi,wj,norm_fcon,system,mirol)
  CASE ("damper_torque")
    CALL sph_sph_damp_torq_roll_resist_w_thcond(i,j,vi,vj,wi,wj,nij,norm_fcon,system,mirol)
  CASE ("elastic_plastic_spring_damper_torque")
    CALL sph_sph_epsd_torq_roll_resist_w_thcond(i,j,wi,wj,overlap,norm_fcon,Ei,Ej,system,mirol)
  CASE DEFAULT
    WRITE(error_unit,*) " Rolling resistance model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE compute_sph_sph_roll_resist_mom_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_sph_wall_roll_resist_mom_w_thcond(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,Ei,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, norm_fcon, Ei
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout, mirol, vrelctang


! select rolling resistance model
SELECT CASE (system%solution_cv%rolling_resistance_model)
  CASE ("constant_torque")
    CALL sph_wall_const_torq_roll_resist_w_thcond(i,k,wi,ww,norm_fcon,system,mirol)
  CASE ("damper_torque")
    CALL sph_wall_damp_torq_roll_resist_w_thcond(i,k,vi,vw,wi,ww,nout,norm_fcon,system,mirol)
  CASE ("elastic_plastic_spring_damper_torque")
    CALL sph_wall_epsd_torq_roll_resist_w_thcond(i,k,wi,ww,overlap,norm_fcon,Ei,system,mirol)
  CASE DEFAULT
    WRITE(error_unit,*) " Rolling resistance model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE compute_sph_wall_roll_resist_mom_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_sph_const_torq_roll_resist_w_thcond(i,j,wi,wj,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: wi, wj
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat, jmat
REAL(KIND=real8)   :: Ri, Rj, Req, miri, mirj, mir, norm_wrel
REAL(KIND=real8), DIMENSION(3) :: wrel


! get particles´ data
Ri = system%particle(i)%radius
Rj = system%particle(j)%radius

! get particles material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirj = system%material_set(jmat)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirj)

! compute equivalent radius for the pair i-j
Req = Ri*Rj/(Ri+Rj)

! compute relative spin
wrel = wi - wj
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*(1.0_real8/norm_wrel)*wrel
END IF

END SUBROUTINE sph_sph_const_torq_roll_resist_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_wall_const_torq_roll_resist_w_thcond(i,k,wi,ww,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: wi, ww
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat
REAL(KIND=real8)   :: Ri, Req, miri, mirw, mir, norm_wrel
REAL(KIND=real8), DIMENSION(3) :: wrel


! get particle´s data
Ri = system%particle(i)%radius

! get particle´s material-set number
imat = system%particle(i)%material_set_number

! get rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirw = system%wall(k)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirw)

! compute equivalent radius for the pair i-k
Req = Ri

! compute relative spin
wrel = wi !- ww
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*(1.0_real8/norm_wrel)*wrel
END IF

END SUBROUTINE sph_wall_const_torq_roll_resist_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_sph_damp_torq_roll_resist_w_thcond(i,j,vi,vj,wi,wj,nij,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat, jmat
REAL(KIND=real8)   :: Ri, Rj, Req, miri, mirj, mir, norm_wrel, norm_vrelctang
REAL(KIND=real8), DIMENSION(3) :: wrel, rci, rcj, vrelc, vrelctang


! get particles´ data
Ri = system%particle(i)%radius
Rj = system%particle(j)%radius

! get particles material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirj = system%material_set(jmat)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirj)

! compute equivalent radius for the pair i-j
Req = Ri*Rj/(Ri+Rj)

! compute relative spin
wrel = wi - wj
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! build contact point radial vectors
rci = Ri*nij
rcj = -Rj*nij
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (Ri-overlap)*nij
!rcj = -(Rj-overlap)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute tangential relative velocity at contact point (due to rotation only)
!vrelc = vi + (wi .vector. rci) - vj - (wj .vector. rcj)
vrelc = (wi .vector. rci) - (wj .vector. rcj)
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*norm_vrelctang*(1.0_real8/norm_wrel)*wrel  !((wi .vector. rci) - (wj .vector. rcj))
END IF

END SUBROUTINE sph_sph_damp_torq_roll_resist_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_wall_damp_torq_roll_resist_w_thcond(i,k,vi,vw,wi,ww,nout,norm_fcon,system,mirol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: norm_fcon
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mirol

! local variables
INTEGER(KIND=int4) :: imat
REAL(KIND=real8)   :: Ri, Req, miri, mirw, mir, norm_wrel, norm_vrelctang
REAL(KIND=real8), DIMENSION(3) :: wrel, rci, vrelc, vrelctang


! get particle´s data
Ri = system%particle(i)%radius

! get particle´s material-set number
imat = system%particle(i)%material_set_number

! get rolling resistance coefficients
miri = system%material_set(imat)%rolling_resistance_coeff
mirw = system%wall(k)%rolling_resistance_coeff

! compute rolling_resistance coefficient for the contacting pair (average)
mir = 0.5_real8*(miri+mirw)

! compute equivalent radius for the pair i-k
Req = Ri

! compute relative spin
wrel = wi !- ww
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! build contact point radial vector
rci = Ri*nout
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (Ri-overlap)*nout
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute tangential relative velocity at contact point (due to rotation only)
!vrelc = vi + (wi .vector. rci) - vw
vrelc = (wi .vector. rci)
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute rolling resistance moment
IF (norm_wrel==0.0_real8) THEN
  mirol = 0.0_real8
ELSE
  mirol = - mir*Req*norm_fcon*norm_vrelctang*(1.0_real8/norm_wrel)*wrel  !*(-vrelc .vector. nout)/SQRT(DOT_PRODUCT( (vrelc .vector. nout) , (vrelc .vector. nout) ))    !mir*Req*norm_fcon*vrelctang   !(wi .vector. rci) !-ww*Rw
END IF

END SUBROUTINE sph_wall_damp_torq_roll_resist_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_sph_epsd_torq_roll_resist_w_thcond(i,j,wi,wj,overlap,norm_fcon,Ei,Ej,system,mrol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, norm_fcon, Ei, Ej
REAL(KIND=real8), DIMENSION(3) :: wi, wj
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mrol

! local variables
INTEGER(KIND=int4) :: imat, jmat, s, cp, jaddr, iaddr
REAL(KIND=real8)   :: Ri, mi, ji, Rj, mj, jj, Req, ni_i, ni_j, miri, mirj, mir, xsiri, xsirj, xsir, dt, &
                      Gi, Gj, Geq, norm_wrel, kr, jr, dr, aux_ji, aux_jj, norm_mrol_trial    
REAL(KIND=real8), DIMENSION(3) :: wrel, sij_trial, sij, dxij_trial, xijr, xij_trial, xij, mrol_trial
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, two=2.0_real8


! get particles´ data
Ri = system%particle(i)%radius
mi = system%particle(i)%mass
ji = system%particle(i)%inertia
Rj = system%particle(j)%radius
mj = system%particle(j)%mass
jj = system%particle(j)%inertia

! get particles´ material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ material properties
!Ei = system%material_set(imat)%elasticity_modulus  <= this comes as argument (since it´s temperature-dependent)
ni_i = system%material_set(imat)%poisson_coeff
miri = system%material_set(imat)%rolling_resistance_coeff
xsiri = system%material_set(imat)%rolling_resistance_damping_ratio
!Ej = system%material_set(jmat)%elasticity_modulus  <= this comes as argument (since it´s temperature-dependent)
ni_j = system%material_set(jmat)%poisson_coeff
mirj = system%material_set(jmat)%rolling_resistance_coeff
xsirj = system%material_set(jmat)%rolling_resistance_damping_ratio

! get particle j´s address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i´s address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! get time step size
s = system%solution_cv%current_step
dt = system%solution_cv%step(s)%current_dt

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute equivalent properties for the pair i-j
Req = Ri*Rj/(Ri+Rj)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)

! compute rolling_resistance parameters for the contacting pair (average)
mir = 0.5_real8*(miri+mirj)
xsir = 0.5_real8*(xsiri+xsirj)

! compute rolling resistance spring stiffness
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!kr = Req*norm_fcon
!kr = 10000_real8*kr
!kr = mir*Req*norm_fcon/0.05_real8
kr = 8.0_real8*Geq*SQRT(Req*overlap)*Req*Req
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute rolling resistance dashpot constant
aux_ji = ji + mi*Ri*Ri
aux_jj = jj + mj*Rj*Rj
jr = one/((one/aux_ji) + (one/aux_jj))
dr = 2.0_real8*xsir*SQRT(jr*kr)

! compute particles´ relative spin
wrel = wi - wj   ! CAUTION: think more about this!
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute trial incremental elongation of rolling resistance spring
dxij_trial = wrel*dt

! compute trial total elongation of rolling resistance spring
xijr = system%particle(i)%contacting_particles(jaddr)%rolresspring_elongation
!xt_accum_proj = xt_accum - DOT_PRODUCT(xt_accum,nij)*nij
!xt_trial = xt_accum_proj + dxt_trial
xij_trial = xijr + dxij_trial

! compute trial rolling resistance moment
mrol_trial = - kr*xij_trial - dr*wrel

! compute trial rolling resistance moment direction
norm_mrol_trial = SQRT(DOT_PRODUCT(mrol_trial,mrol_trial))
IF (norm_mrol_trial/=0.0_real8) THEN
  sij_trial = (one/norm_mrol_trial)*mrol_trial
ELSE
  sij_trial = 0.0_real8
END IF

! check for stick or slip (verify if trial rolling resistance moment violates the maximum possible value)
IF (norm_mrol_trial<=(mir*Req*norm_fcon)) THEN
  mrol = mrol_trial
  xij = xij_trial
  sij = sij_trial
  !"stick region"
ELSE
  sij = sij_trial
  mrol = mir*Req*norm_fcon*sij
  xij = -(one/kr)*(mrol - dr*wrel)
  !"slip region"
END IF

! update pair´s rolling resistance spring elongation
system%particle(i)%contacting_particles(jaddr)%rolresspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%rolresspring_elongation = xij

END SUBROUTINE sph_sph_epsd_torq_roll_resist_w_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_wall_epsd_torq_roll_resist_w_thcond(i,k,wi,ww,overlap,norm_fcon,Ei,system,mrol)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, norm_fcon, Ei
REAL(KIND=real8), DIMENSION(3) :: wi, ww
REAL(KIND=real8), DIMENSION(3), INTENT(OUT) :: mrol

! local variables
INTEGER(KIND=int4) :: imat, s, j, waddr
REAL(KIND=real8)   :: Ri, mi, ji, Req, ni_i, miri, xsiri, mirw, xsirw, Gi, Geq, mir, xsir, &
                      dt, kr, norm_wrel, jr, dr, norm_mrol_trial
REAL(KIND=real8), DIMENSION(3) :: wrel, sij_trial, sij, dxiw_trial, xiwr, xiw_trial, xiw, mrol_trial
REAL(KIND=real8), PARAMETER    :: one=1.0_real8


! get particle´s data
Ri = system%particle(i)%radius
mi = system%particle(i)%mass
ji = system%particle(i)%inertia

! get particle´s material-set number
imat = system%particle(i)%material_set_number

! get particle´s and wall´s material properties
!Ei = system%material_set(imat)%elasticity_modulus  <= this comes as argument (since it´s temperature-dependent)
ni_i = system%material_set(imat)%poisson_coeff
miri = system%material_set(imat)%rolling_resistance_coeff
xsiri = system%material_set(imat)%rolling_resistance_damping_ratio
mirw = system%wall(k)%rolling_resistance_coeff
xsirw = system%wall(k)%rolling_resistance_damping_ratio

! get wall´s address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! get time step size
s = system%solution_cv%current_step
dt = system%solution_cv%step(s)%current_dt

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(1.0_real8+ni_i))

! compute equivalent properties for the pair i-k
Req = Ri
Geq = Gi/(2.0_real8-ni_i)

! compute rolling_resistance parameters for the contacting pair (average)
mir = 0.5_real8*(miri+mirw)
xsir = 0.5_real8*(xsiri+xsirw)

! compute rolling resistance spring stiffness
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!kr = Req*norm_fcon
!kr = 10000_real8*kr
!kr = mir*Req*norm_fcon/0.05_real8
kr = 8.0_real8*Geq*SQRT(Req*overlap)*Req*Req
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute rolling resistance dashpot constant
jr = ji + mi*Ri*Ri
dr = 2.0_real8*xsir*SQRT(jr*kr)

! compute particle-wall relative spin
wrel = wi !- ww   ! CAUTION: think more about this!
norm_wrel = SQRT(DOT_PRODUCT(wrel,wrel))

! compute trial incremental elongation of rolling resistance spring
dxiw_trial = wrel*dt

! compute trial total elongation of rolling resistance spring
xiwr = system%particle(i)%contacting_walls(waddr)%rolresspring_elongation
!xt_accum_proj = xt_accum - DOT_PRODUCT(xt_accum,nij)*nij
!xt_trial = xt_accum_proj + dxt_trial
xiw_trial = xiwr + dxiw_trial

! compute trial rolling resistance moment
mrol_trial = - kr*xiw_trial - dr*wrel

! compute trial rolling resistance moment direction
norm_mrol_trial = SQRT(DOT_PRODUCT(mrol_trial,mrol_trial))
IF (norm_mrol_trial/=0.0_real8) THEN
  sij_trial = (one/norm_mrol_trial)*mrol_trial
ELSE
  sij_trial = 0.0_real8
END IF

! check for stick or slip (verify if trial rolling resistance moment violates the maximum possible value)
IF (norm_mrol_trial<=(mir*Req*norm_fcon)) THEN
  mrol = mrol_trial
  xiw = xiw_trial
  sij = sij_trial
  !"stick region"
ELSE
  sij = sij_trial
  mrol = mir*Req*norm_fcon*sij
  xiw = -(one/kr)*(mrol - dr*wrel)
  !"slip region"
END IF

! update pair´s rolling resistance spring elongation
system%particle(i)%contacting_walls(waddr)%rolresspring_elongation = xiw

END SUBROUTINE sph_wall_epsd_torq_roll_resist_w_thcond
!!---------------------------------------------------------------------------------------------------------------------

! add here subroutines for other rolling resistance models with temperature effects

END MODULE rolling_resistance_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE impulse_based_contact

USE particles_classes
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_contact_by_impulses(i,j,overlap,vi,vj,wi,wj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat
REAL(KIND=real8)   :: delta_t, dt, ei, ej, e, misi, misj, mis, midi, midj, mid, mi, mj, dt1, dt2, vin_t, vjn_t, vin_tpdt, vjn_tpdt, &
                      epsin_compr, epsjn_compr, epsin_recov, epsjn_recov, epsin_bar, epsjn_bar, Deltaij_compr, Deltaij_recov, &
                      In_bar, epsit_bar, epsjt_bar, vit_t, vjt_t, It_bar, norm_vrelctang_t, It_epsi, It_epsj, It_fric, f_It_fric
REAL(KIND=real8), DIMENSION(3) :: vi_t, vj_t, wi_t, wj_t, wi_tpdt, wj_tpdt, rci, rcj, vic_t, vjc_t, vrelc_t, vrelctang_t, &
                                  tauij, fcont, ffric, mifric, mjfric, auxi, auxj, auxji
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8                      


! Note: this subroutine needs to be cleaned and optimized w.r.t dummy arguments and local variables

! get time step info
s = system%solution_cv%current_step
delta_t = system%solution_cv%step(s)%current_dt

! compute contact duration
dt = system%solution_cv%step(s)%collisions_duration_parameter*delta_t
system%particle(i)%collision_duration = dt
system%particle(j)%collision_duration = dt    

! get particles material set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles collision-related material properties
ei = system%material_set(imat)%coefficient_of_restitution
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
ej = system%material_set(jmat)%coefficient_of_restitution
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! compute collision material parameters (average) for the colliding pair
e = half*(ei+ej)
mis = half*(misi+misj)
mid = half*(midi+midj)

! get particles´ masses
mi = system%particle(i)%mass
mj = system%particle(j)%mass

! get particles´ velocity and spin vectors before collision
vi_t = system%particle(i)%velocity
vj_t = system%particle(j)%velocity
IF (system%general_cv%rotational_dofs=="on") THEN
  wi_t = system%particle(i)%spin
  wj_t = system%particle(j)%spin
ELSE
  wi_t = 0.0_real8
  wj_t = 0.0_real8
END IF    

! compute particles´ normal velocities before collision
vin_t = DOT_PRODUCT(vi_t,nij)
vjn_t = DOT_PRODUCT(vj_t,nij)

! compute duration of compressive and recovery stages
dt1 = dt/(one+e)
dt2 = e*dt/(one+e)

! compute external normal impulses during compressive and recovery stages
epsin_compr = DOT_PRODUCT(system%particle(i)%epsilon_force,nij)*dt1
epsjn_compr = DOT_PRODUCT(system%particle(j)%epsilon_force,nij)*dt1
epsin_recov = DOT_PRODUCT(system%particle(i)%epsilon_force,nij)*dt2
epsjn_recov = DOT_PRODUCT(system%particle(j)%epsilon_force,nij)*dt2

! compute external normal impulses during whole collision
epsin_bar = DOT_PRODUCT(system%particle(i)%epsilon_force,nij)*dt
epsjn_bar = DOT_PRODUCT(system%particle(j)%epsilon_force,nij)*dt

! compute auxiliary variables
Deltaij_compr = (one/mi)*epsin_compr - (one/mj)*epsjn_compr
Deltaij_recov = (one/mi)*epsin_recov - (one/mj)*epsjn_recov

! compute normal velocities after collision
vin_tpdt = ( mi*vin_t + mj*(vjn_t - e*(vin_t-vjn_t)) + (epsin_bar + epsjn_bar) - mj*(e*Deltaij_compr - Deltaij_recov) )/(mi+mj)
vjn_tpdt = vin_tpdt - Deltaij_recov + e*(vin_t - vjn_t + Deltaij_compr)

! compute contact force due to collision
In_bar = mi*(vin_tpdt - vin_t)/dt - epsin_bar/dt
fcont = In_bar*nij
system%particle(i)%contact_force = system%particle(i)%contact_force + fcont
system%particle(j)%contact_force = system%particle(j)%contact_force - fcont

! update particles epsilon force vectors (needed if particles have other contact pairs)
system%particle(i)%epsilon_force = system%particle(i)%epsilon_force + fcont
system%particle(j)%epsilon_force = system%particle(j)%epsilon_force - fcont


! build particles´ contact point radial vectors
rci = system%particle(i)%radius*nij
rcj = -system%particle(j)%radius*nij

! compute contact point velocity vectors before collision
vic_t = vi_t + (wi_t .vector. rci)
vjc_t = vj_t + (wj_t .vector. rcj)
vrelc_t = vic_t - vjc_t

! compute tangential relative velocity at the contact point before collision
vrelctang_t = vrelc_t - (DOT_PRODUCT(vrelc_t,nij))*nij
norm_vrelctang_t = SQRT(DOT_PRODUCT(vrelctang_t,vrelctang_t))

! compute friction force when the relative tangential velocity is zero (friction may arise due to ENV forces)
friction_with_zero_relative_vel: IF (norm_vrelctang_t==0.0_real8) THEN

  ! friction force in this case is currently being set to zero (this may be clarified in Berkeley, July 2014  
  ffric = 0.0_real8

  !! if the tangential relative velocity is zero then the direction of ffric will be the opposite direction of the tangent component of ENV;
  !! aqui: supor primeiro que ENV é baixo e que portanto há stick durante toda a colisão, i.e., ffric = -ENV*tauij<=mis*In_bar*tauij;
  !! se ENV>mis*In_bar, então há enduring slip e logo ffric = mid*In_bar*tauij; 
  !! em ambos os casos, a direção tauij é a direção da componente tangente (tangente em relação a nij) de ENV.
  !! Notice: NFs do not generate friction because NFs are central forces, i.e. they always act in the direction of nij!

  !! get environment forces acting on particle i during the collision
  !! check whether to use something like epsilon_force (which is continously updated) or system%part(i)%%env (which is not) as below!!!
  !finfenv = system%particle(i)%environmental_force ! system%particle(i)%epsilon_force
  
  !! compute tangential resultant of the environment forces
  !ftang = finfenv -(DOT_PRODUCT(finfenv,n_out))*n_out
  !norm_ftang = SQRT(DOT_PRODUCT(ftang,ftang))
  
  !! compute friction force
  !IF (norm_ftang==0.0_real8) THEN
  !  ffric = 0.0_real8
  !ELSE IF (norm_ftang<=mis*ABS(In_bar)) THEN   ! stick occurs
  !  ffric = ??? (derive an expression for this force: maybe based on the balance of momentum, analogous to the case where the normal relative velocity is not zero)  !-ftang  !drawback (same as happens in the normal direction when there are NFs and ENV forces): since this force is smeared over delta_t, it doesn´t nullify the action of fienv and finf, which lasts both delta_t; as a consequence, there will be a sliding over delta_t!
  !  system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
  !  !system%wall(w)%friction_force = system%wall(w)%friction_force - ffric 
  !ELSE   ! continuous sliding occurs
  !  tauij = (one/norm_ftang)*ftang
  !  ffric = -mid*ABS(In_bar)*tauij
  !  system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
  !  !system%wall(w)%friction_force = system%wall(w)%friction_force - ffric 
  !END IF
 
  !! update epsilon force vector (needed if particle and wall have other contact pairs)
  !system%particle(i)%epsilon_force = system%particle(i)%epsilon_force + ffric
  !!system%wall(w)%total_force = system%wall(w)%total_force - ffric

END IF friction_with_zero_relative_vel

! compute friction force when the relative tangential velocity is nonzero
friction_with_nonzero_relative_vel: IF (norm_vrelctang_t/=0.0_real8) THEN

  ! compute tangential direction before collision
  tauij = (one/norm_vrelctang_t)*vrelctang_t
  
  ! compute tangential velocities before collision
  vit_t = DOT_PRODUCT(vi_t,tauij)
  vjt_t = DOT_PRODUCT(vj_t,tauij)

  ! compute external tangential impulses during collision
  !epsit_bar = DOT_PRODUCT(system%particle(i)%epsilon_force,tauij)*dt
  !epsjt_bar = DOT_PRODUCT(system%particle(j)%epsilon_force,tauij)*dt
  It_epsi = DOT_PRODUCT(system%particle(i)%epsilon_force,tauij)*dt
  It_epsj = DOT_PRODUCT(system%particle(j)%epsilon_force,tauij)*dt

  ! get iterative (i.e. still not converged) spins at t+dt
  wi_tpdt = wi
  wj_tpdt = wj

  ! compute friction impulse due to collision (in It_fric below, it is assumed first that stick occurs)
  !It_bar = ( (epsjt_bar/mj - epsit_bar/mi) + vjt_t - vit_t ) / ((one/mi + one/mj)*dt)
  auxi = wi_tpdt .vector. rci
  auxj = wj_tpdt .vector. rcj
  auxji = auxj - auxi
  It_fric = ( (It_epsj/mj - It_epsi/mi) + (vjt_t - vit_t) + DOT_PRODUCT(auxji,tauij) ) / (one/mi + one/mj)
  
  ! compute friction force due to collision
  f_It_fric = It_fric/dt
  IF ( ABS(f_It_fric)<=(mis*ABS(In_bar)) ) THEN ! stick occurs
    ffric = f_It_fric*tauij
    system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
    system%particle(j)%friction_force = system%particle(j)%friction_force - ffric
  ELSE
    ffric = -mid*ABS(In_bar)*tauij   ! continuous sliding occurs
    system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
    system%particle(j)%friction_force = system%particle(j)%friction_force - ffric
  END IF                             

  ! update epsilon force vector (needed if particles have other contact pairs)
  system%particle(i)%epsilon_force = system%particle(i)%epsilon_force + ffric
  system%particle(j)%epsilon_force = system%particle(j)%epsilon_force - ffric
  
END IF friction_with_nonzero_relative_vel

! compute moment due to fritcion force
mifric = rci .vector. ffric
mjfric = rcj .vector. (-ffric)
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifric
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfric

END SUBROUTINE sphere_sphere_contact_by_impulses
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_contact_by_impulses(i,k,overlap,vi,vw,wi_tpdt,ww_tpdt,n_out,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi_tpdt, ww_tpdt, n_out

! local variables
INTEGER(KIND=int4) :: s, imat
REAL(KIND=real8)   :: delta_t, dt, ei, ew, e, misi, misw, mis, midi, midw, mid, mi, mj, dt1, dt2, vin_t, vwn_t, &
                      epsin_compr, epsin_recov, epsin_bar, Deltaiw_compr, Deltaiw_recov, vin_tpdt, vwn_tpdt, &
                      In_bar, vit_t, vwt_t, epsit_bar, epsjt_bar, It_bar, norm_vrelctang_t, It_epsi, It_fric, f_It_fric
REAL(KIND=real8), DIMENSION(3) :: vi_t, vw_t, wi_t, rci, vic_t, vrelc_t, vrelctang_t, tauiw, fcont, ffric, mifric, auxi
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8                      


! Note: this subroutine needs to be cleaned and optimized w.r.t dummy arguments and local variables

! get time step info
s = system%solution_cv%current_step
delta_t = system%solution_cv%step(s)%current_dt

! compute contact duration
dt = system%solution_cv%step(s)%collisions_duration_parameter*delta_t
system%particle(i)%collision_duration = dt

! get particle mass and material-set number
mi = system%particle(i)%mass
imat = system%particle(i)%material_set_number

! get particle and wall collision-related material properties
ei = system%material_set(imat)%coefficient_of_restitution
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
ew = system%wall(k)%coefficient_of_restitution
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! compute collision parameters (average) for the pair i-k
e = half*(ei+ew)
mis = half*(misi+misw)
mid = half*(midi+midw)

! get velocity and spin vectors before collision
vi_t = system%particle(i)%velocity
vw_t = system%wall(k)%velocity
IF (system%general_cv%rotational_dofs=="on") THEN
  wi_t = system%particle(i)%spin
  !ww_t = 0.0_real8  !system%wall(k)%spin
ELSE
  wi_t = 0.0_real8
  !ww_t = 0.0_real8
END IF    

! compute normal velocities before collision
vin_t = DOT_PRODUCT(vi_t,n_out)
vwn_t = DOT_PRODUCT(vw_t,n_out)

! compute duration of compressive and recovery stages
dt1 = dt/(one+e)
dt2 = e*dt/(one+e)

! compute external normal impulses during compressive and recovery stages
epsin_compr = DOT_PRODUCT(system%particle(i)%epsilon_force,n_out)*dt1
epsin_recov = DOT_PRODUCT(system%particle(i)%epsilon_force,n_out)*dt2

! compute average external normal impulses
epsin_bar = DOT_PRODUCT(system%particle(i)%epsilon_force,n_out)*dt

! compute auxiliary variables
Deltaiw_compr = (one/mi)*epsin_compr !- (one/mj)*epsjn_compr
Deltaiw_recov = (one/mi)*epsin_recov !- (one/mj)*epsjn_recov

! compute normal velocities after collision
!vin_tpdt = ( mi*vin_t + mj*(vjn_t - e*(vin_t-vjn_t)) + (epsin_bar + epsjn_bar) - mj*(e*Deltaij_compr - Deltaij_recov) )/(mi+mj)
!vwn_tpdt = vin_tpdt - Deltaij_recov + e*(vin_t - vjn_t + Deltaij_compr)
vin_tpdt = Deltaiw_recov - e*(vin_t - vwn_t + Deltaiw_compr)  !+vwn_t
vwn_tpdt = vwn_t

! compute contact force due to collision
In_bar = mi*(vin_tpdt - vin_t)/dt - epsin_bar/dt
fcont = In_bar*n_out
system%particle(i)%contact_force = system%particle(i)%contact_force + fcont
!system%wall(k)%contact_force = system%wall(k)%contact_force - fcont

! update epsilon force vector (needed if particle and wall have other contact pairs)
system%particle(i)%epsilon_force = system%particle(i)%epsilon_force + fcont
!system%wall(k)%total_force = system%wall(k)%total_force - fcont


! build particle´s contact point radial vector
rci = system%particle(i)%radius*n_out
!rcw = Rw*n_out

! compute contact point velocity vectors before collision
vic_t = vi_t + (wi_t .vector. rci)
!vwc_t = vw_t + (ww_t .vector. rcw)
vrelc_t = vic_t - vw_t

! compute tangential relative velocity at the contact point before collision
vrelctang_t = vrelc_t - (DOT_PRODUCT(vrelc_t,n_out))*n_out
norm_vrelctang_t = SQRT(DOT_PRODUCT(vrelctang_t,vrelctang_t))

! compute friction force when the relative tangential velocity is zero (friction may arise due to ENV forces)
friction_with_zero_relative_vel: IF (norm_vrelctang_t==0.0_real8) THEN

  ! friction force is currently being set to zero in this case
  ffric = 0.0_real8

  !! if the tangential relative velocity is zero then the direction of ffric is given by the opposite direction of the tangent component of ENV;
  !! aqui: supor primeiro que ENV é baixo e que portanto há stick durante toda a colisão, i.e., ffric = -ENV*tauij<=mis*In_bar*tauij;
  !! se ENV>mis*In_bar, então há enduring slide e logo ffric = mid*In_bar*tauij; 
  !! em ambos os casos, a direção tauij é a direção da componente tangente (tangente em relação a nij) de ENV.
  ! Note: NFs do not generate friction because they are central forces

  !! get environment forces acting on particle i during the collision
  !! check whether to use something analogous to epsilon_force (which is continously updated) or sys%part(i)%env below!!!
  !finfenv = system%particle(i)%environmental_force
  
  !! compute tangential resultant of the environment forces
  !ftang = finfenv -(DOT_PRODUCT(finfenv,n_out))*n_out
  !norm_ftang = SQRT(DOT_PRODUCT(ftang,ftang))
  
  !! compute friction force
  !IF (norm_ftang==0.0_real8) THEN
  !  ffric = 0.0_real8
  !ELSE IF (norm_ftang<=mis*ABS(In_bar)) THEN   ! stick occurs
  !  ffric = ??? (derive an expression for this force: maybe based on the balance of momentum, analogous to the case where the normal relative velocity is not zero)  !-ftang  !drawback (same as happens in the normal direction when there are NFs and ENV forces): since this force is smeared over delta_t, it doesn´t nullify the action of fienv and finf, which lasts both delta_t; as a consequence, there will be a sliding over delta_t!
  !  system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
  !  !system%wall(w)%friction_force = system%wall(w)%friction_force - ffric 
  !ELSE   ! continuous sliding occurs
  !  tauij = (one/norm_ftang)*ftang
  !  ffric = -mid*ABS(In_bar)*tauij
  !  system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
  !  !system%wall(w)%friction_force = system%wall(w)%friction_force - ffric 
  !END IF
 
  !! update epsilon force vector (needed if particle and wall have other contact pairs)
  !system%particle(i)%epsilon_force = system%particle(i)%epsilon_force + ffric
  !!system%wall(w)%total_force = system%wall(w)%total_force - ffric

END IF friction_with_zero_relative_vel

! compute friction force when the relative tangential velocity is nonzero
friction_with_nonzero_relative_vel: IF (norm_vrelctang_t/=0.0_real8) THEN

  ! compute tangential direction before collision
  tauiw = (one/norm_vrelctang_t)*vrelctang_t
  
  ! compute tangential velocities before collision
  vit_t = DOT_PRODUCT(vi_t,tauiw)
  vwt_t = DOT_PRODUCT(vw_t,tauiw)

  ! compute external tangential impulse on particle i during collision
  !epsit_bar = DOT_PRODUCT(system%particle(i)%epsilon_force,tauiw)*dt
  It_epsi = DOT_PRODUCT(system%particle(i)%epsilon_force,tauiw)*dt

  ! compute friction impulse due to collision (in It_fric below, it is assumed first that stick occurs)
  !It_bar = (mi*(vwt_t - vit_t) - epsit_bar)/dt         
  auxi = wi_tpdt .vector. rci
  It_fric = mi*(vwt_t - vit_t) - It_epsi - mi*DOT_PRODUCT(auxi,tauiw)

  ! compute friction force due to collision
  f_It_fric = It_fric/dt
  IF ( ABS(f_It_fric)<=(mis*ABS(In_bar)) ) THEN  ! stick occurs
    ffric = f_It_fric*tauiw
    system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
    !system%wall(k)%friction_force = system%wall(k)%friction_force - ffric 
  ELSE
    ffric = -mid*ABS(In_bar)*tauiw   ! continuous slide occurs
    system%particle(i)%friction_force = system%particle(i)%friction_force + ffric
    !system%wall(k)%friction_force = system%wall(k)%friction_force - ffric 
  END IF                             

  ! update epsilon force vector (needed if particle and wall have other contact pairs)
  system%particle(i)%epsilon_force = system%particle(i)%epsilon_force + ffric
  !system%wall(k)%total_force = system%wall(k)%total_force - ffric

END IF friction_with_nonzero_relative_vel

! compute moment due to fritcion force
mifric = rci .vector. ffric
!mwfric = rcw .vector. (-ffric)
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifric
!system%wall(k)%friction_moment = system%wall(k)%friction_moment + mwfric

END SUBROUTINE sphere_wall_contact_by_impulses
!!---------------------------------------------------------------------------------------------------------------------


END MODULE impulse_based_contact

!!=============================================================================================================================
!!=============================================================================================================================


MODULE hertz_contact_with_contslide_friction

USE particles_classes
USE array_of_vectors_class
USE rolling_resistance_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_contact_by_hertz_with_contslide(i,j,overlap,vi,vj,wi,wj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, Gi, Gj, Ri, Rj, mi, mj, Eeq, Geq, Req, meq, zetai, zetaj, zeta, &
                      misi, misj, mis, midi, midj, mid, kn, d, cn, vin, vjn, vreln, norm_vrelctang, norm_fcon, kt
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, mjfri, vrelctang, tauij, rci, rcj, vic, vjc, vrelc
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8


! get particles radius and mass
Ri = system%particle(i)%radius
Rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass

! get particles material-set numbers
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
zetai = system%material_set(imat)%contact_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
zetaj = system%material_set(jmat)%contact_damping_ratio
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! compute particles tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute contact parameters (average) for the contacting pair 
zeta = half*(zetai+zetaj)
mis = half*(misi+misj)
mid = half*(midi+midj)

! compute equivalent properties for the pair i-j
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)
Req = Ri*Rj/(Ri+Rj)
meq = mi*mj/(mi+mj)

! compute normal damping variables
kn = 2.0_real8*Eeq*SQRT(Req)*SQRT(overlap)
d = 2.0_real8*SQRT(kn*meq)
cn = zeta*d

! compute normal velocities for normal damping
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! compute normal contact force (includes normal dissipative damping)
fcon = -1.3333333333333333_real8*Eeq*SQRT(Req)*(overlap**1.5_real8)*nij - (cn*vreln)*nij
IF (DOT_PRODUCT(fcon,nij)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! update particles´ contact forces
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%particle(j)%contact_force = system%particle(j)%contact_force - fcon

! compute friction spring stiffness
kt = 8.0_real8*Geq*SQRT(Req)*SQRT(overlap)

! build particles´ contact point radial vectors
rci = Ri*nij
rcj = -Rj*nij

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
vjc = vj + (wj .vector. rcj)
vrelc = vic - vjc

! compute tangential relative velocity at the contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute tangential friction force (continuous slide friction)
IF (norm_vrelctang==0.0_real8) THEN
  ffri = 0.0_real8
ELSE
  tauij = (one/norm_vrelctang)*vrelctang
  norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))
  ffri = -mid*norm_fcon*tauij
END IF

! update particles´ friction forces
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%particle(j)%friction_force = system%particle(j)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri
mjfri = rcj .vector. (-ffri)

! update particles´ friction moments
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfri

END SUBROUTINE sphere_sphere_contact_by_hertz_with_contslide
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_contact_by_hertz_with_contslide(i,k,overlap,vi,vw,wi,ww,n_out,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, n_out

! local variables
INTEGER(KIND=int4) :: imat
REAL(KIND=real8)   :: Ei, ni_i, Gi, Ri, mi, Eeq, Geq, Req, meq, misi, misw, mis, midi, midw, mid, vin, vwn, &
                      vit_t, vwt_t, norm_vrelctang, zetai, zetaw, zeta, kn, d, cn, vreln, norm_fcon, kt
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, vrelctang, tauiw, vi_t, vw_t, rci, vic, vrelc
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8


! get particle data
mi = system%particle(i)%mass
Ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle and wall contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
zetai = system%material_set(imat)%contact_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
zetaw = system%wall(k)%contact_damping_ratio
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))

! compute contact parameters (average) for the pair i-k
zeta = half*(zetai+zetaw)
mis = half*(misi+misw)
mid = half*(midi+midw)

! compute equivalent properties for the pair i-k
Eeq = Ei/(one-ni_i*ni_i)
Geq = Gi/(two-ni_i)
Req = Ri
meq = mi

! compute normal damping variables
kn = 2.0_real8*Eeq*SQRT(Req)*SQRT(overlap)
d = 2.0_real8*SQRT(kn*meq)
cn = zeta*d

! compute normal velocities
vin = DOT_PRODUCT(vi,n_out)
vwn = DOT_PRODUCT(vw,n_out)
vreln = vin - vwn

! compute normal contact force (includes normal dissipative damping)
fcon = -1.3333333333333333_real8*Eeq*SQRT(Req)*(overlap**1.5_real8)*n_out - (cn*vreln)*n_out 
IF (DOT_PRODUCT(fcon,n_out)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! update particle´s contact force
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
!system%wall(k)%contact_force = system%wall(k)%contact_force - fcon

! compute friction spring stiffness
kt = 8.0_real8*Geq*SQRT(Req)*SQRT(overlap)

! build particle´s contact point radial vector
rci = Ri*n_out

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
!vwc = vw + (ww .vector. rcw)
vrelc = vic - vw

! compute tangential relative velocity at the contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,n_out))*n_out
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute tangential friction force
IF (norm_vrelctang==0.0_real8) THEN
  ffri = 0.0_real8
ELSE
  tauiw = (one/norm_vrelctang)*vrelctang
  norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))
  ffri = -mid*norm_fcon*tauiw
END IF

! update particle´s friction force
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
!system%wall(k)%friction_force = system%wall(k)%friction_force - ffri

! compute moment due to fritcion force
mifri = rci .vector. ffri

! update particle´s friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri

END SUBROUTINE sphere_wall_contact_by_hertz_with_contslide
!!---------------------------------------------------------------------------------------------------------------------


END MODULE hertz_contact_with_contslide_friction

!!=============================================================================================================================
!!=============================================================================================================================


MODULE hertz_contact_with_stickslip_friction

USE particles_classes
USE array_of_vectors_class
USE rolling_resistance_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_contact_by_hertz_with_stickslip(i,j,overlap,vi,vj,wi,wj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, Gi, Gj, ri, rj, mi, mj, Eeq, Geq, req, meq, deltat, xsini, xsinj, xsin, &
                      xsiti, xsitj, xsit, misi, misj, mis, midi, midj, mid, kn, dn, vin, vjn, vreln, &
                      norm_fcon, kt, dt, norm_ffri_trial
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, mjfri, vrelctang, tij, rci, rcj, vic, vjc, vrelc, &
                                  xij, xij_trial, dxij_trial, ffri_trial, tij_trial, xijr, xijr_proj, mirol, fadh
REAL(KIND=real8), PARAMETER :: one=1.0_real8, half=0.5_real8, two=2.0_real8, fourthirds=1.3333333333333333_real8                     


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
xsinj = system%material_set(jmat)%contact_damping_ratio
xsitj = system%material_set(jmat)%friction_damping_ratio
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! get particle j address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute contact parameters for the contacting pair (average)
xsin = half*(xsini+xsinj)
xsit = half*(xsiti+xsitj)
mis = half*(misi+misj)
mid = half*(midi+midj)

! compute normal contact spring stiffness and dashpot constant
kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! compute normal contact force
fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nij - dn*vreln*nij
IF (DOT_PRODUCT(fcon,nij)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! udpate particles´ total contact forces
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%particle(j)%contact_force = system%particle(j)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 100.0_real8*kt
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particles´ contact point radial vectors
rci = (ri)*nij
rcj = -(rj)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nij
!rcj = -(rj-overlap)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
vjc = vj + (wj .vector. rcj)
vrelc = vic - vjc

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij

! compute trial incremental elongation of friction spring
dxij_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xijr = system%particle(i)%contacting_particles(jaddr)%fricspring_elongation
xijr_proj = xijr - DOT_PRODUCT(xijr,nij)*nij
xij_trial = xijr_proj + dxij_trial

! compute trial friction force
ffri_trial = - kt*xij_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tij_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tij_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xij = xij_trial
 tij = tij_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tij = tij_trial
  ffri = mid*norm_fcon*tij
  xij = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_particles(jaddr)%fricspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%fricspring_elongation = xij

! update particles´ total friction forces
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%particle(j)%friction_force = system%particle(j)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri
mjfri = rcj .vector. (-ffri)

! update particles´ total friction moments
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfri

! compute moment due to rolling resistance
CALL compute_sphere_sphere_rolling_resistance_moment(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moments
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol
system%particle(j)%rolling_resistance_moment = system%particle(j)%rolling_resistance_moment - mirol

END SUBROUTINE sphere_sphere_contact_by_hertz_with_stickslip
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_contact_by_hertz_with_stickslip(i,k,overlap,vi,vw,wi,ww,nout,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, Gi, ri, mi, Eeq, Geq, req, meq, misi, misw, mis, midi, midw, mid, deltat, xsini, &
                      xsinw, xsin, xsiti, xsitw, xsit, vin, vwn, vreln, kn, dn, norm_fcon, kt, dt, norm_ffri_trial
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, vrelctang, tiw, rci, vic, vrelc, xiw, xiw_trial, dxiw_trial, &
                                  xiwr, xiwr_proj, ffri_trial, tiw_trial, mirol, fadh
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8, fourthirds=1.3333333333333333_real8                    


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s and wall´s contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
xsinw = system%wall(k)%contact_damping_ratio
xsitw = system%wall(k)%friction_damping_ratio
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! get wall k address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))

! compute equivalent properties for the contacting pair
Eeq = Ei/(one-ni_i*ni_i)
Geq = Gi/(two-ni_i)
req = ri
meq = mi

! compute contact parameters for the contacting pair (average) 
xsin = half*(xsini+xsinw)
xsit = half*(xsiti+xsitw)
mis = half*(misi+misw)
mid = half*(midi+midw)

! compute normal contact spring stiffness and dashpot constant
kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nout)
vwn = DOT_PRODUCT(vw,nout)
vreln = vin - vwn

! compute normal contact force
fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nout - dn*vreln*nout 
IF (DOT_PRODUCT(fcon,nout)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! update particle´s and wall´s total contact force
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%wall(k)%contact_force = system%wall(k)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 100.0_real8*kt
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particle´s contact point radial vector
rci = ri*nout
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nout
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
!vwc = vw + (ww .vector. rcw)
vrelc = vic - vw

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout

! compute trial incremental elongation of friction spring
dxiw_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xiwr = system%particle(i)%contacting_walls(waddr)%fricspring_elongation
xiwr_proj = xiwr - DOT_PRODUCT(xiwr,nout)*nout
xiw_trial = xiwr_proj + dxiw_trial

! compute trial friction force
ffri_trial = - kt*xiw_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tiw_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tiw_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xiw = xiw_trial
 tiw = tiw_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tiw = tiw_trial
  ffri = mid*norm_fcon*tiw
  xiw = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_walls(waddr)%fricspring_elongation = xiw

! update particle´s and wall´s total friction force
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%wall(k)%friction_force = system%wall(k)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri

! update particle´s total friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri

! compute moment due to rolling resistance
CALL compute_sphere_wall_rolling_resistance_moment(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moment
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol

END SUBROUTINE sphere_wall_contact_by_hertz_with_stickslip
!!---------------------------------------------------------------------------------------------------------------------


END MODULE hertz_contact_with_stickslip_friction

!!=============================================================================================================================
!!=============================================================================================================================


MODULE hertz_contact_with_thornton_friction

USE particles_classes
USE array_of_vectors_class
USE rolling_resistance_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_contact_by_hertz_with_thornton(i,j,overlap,vi,vj,wi,wj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, Gi, Gj, ri, rj, mi, mj, Eeq, Geq, req, meq, deltat, xsini, xsinj, xsin, &
                      xsiti, xsitj, xsit, misi, misj, mis, midi, midj, mid, kn, dn, vin, vjn, vreln, &
                      norm_vrelctang, norm_fcon, kt, dt, norm_ffri_trial
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, mjfri, vrelctang, tij, rci, rcj, vic, vjc, vrelc, &
                                  xij, xij_trial, dxij_trial, ffri_trial, tij_trial, xijr, xijr_proj, mirol
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8                      


! get particles data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
xsinj = system%material_set(jmat)%contact_damping_ratio
xsitj = system%material_set(jmat)%friction_damping_ratio
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! get particle j address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute contact parameters (average) for the contacting pair 
xsin = half*(xsini+xsinj)
xsit = half*(xsiti+xsitj)
mis = half*(misi+misj)
mid = half*(midi+midj)

! compute equivalent properties for the pair i-j
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute normal contact parameters (spring stiffness and damping constant)
kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! compute normal contact force (includes contact damping)
fcon = - 1.3333333333333333_real8*Eeq*SQRT(req)*(overlap**1.5_real8)*nij - (dn*vreln)*nij
IF (DOT_PRODUCT(fcon,nij)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! udpate particles´ contact forces
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%particle(j)%contact_force = system%particle(j)%contact_force - fcon

! build particles´ contact point radial vectors
rci = (ri)*nij
rcj = -(rj)*nij
!ccccccccccccccccccccccccc
!rci = (ri-overlap)*nij
!rcj = -(rj-overlap)*nij
!ccccccccccccccccccccccccc

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
vjc = vj + (wj .vector. rcj)
vrelc = vic - vjc

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute trial incremental elongation of friction spring
dxij_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xijr = system%particle(i)%contacting_particles(jaddr)%fricspring_elongation
xijr_proj = xijr - DOT_PRODUCT(xijr,nij)*nij
xij_trial = xijr_proj + dxij_trial

! compute friction parameters (spring stiffness and damping constant)
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
dt = 2.0_real8*xsit*SQRT(meq*kt)


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! CAUTION: the particle-particle version of thornton´s model still needs to be implemented (only the particle-wall is finished)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


! compute trial friction force
ffri_trial = - kt*xij_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tij_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tij_trial = 0.0_real8
END IF

! compute norm of contact normal force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xij = xij_trial
 tij = tij_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tij = tij_trial
  ffri = mid*norm_fcon*tij
  xij = -(one/kt)*(ffri - dt*vrelctang)  ! (note: this agrees w/ avci´s model)
  !system%particle(i)%friction_status = "slip region"
END IF

! uddate particles´ total friction forces and pair´s friction spring elongation
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%particle(j)%friction_force = system%particle(j)%friction_force - ffri
system%particle(i)%contacting_particles(jaddr)%fricspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%fricspring_elongation = xij

! compute moment due to friction force
mifri = rci .vector. ffri
mjfri = rcj .vector. (-ffri)

! update particles´ total friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfri

! compute moment due to rolling resistance
CALL compute_sphere_sphere_rolling_resistance_moment(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moment
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol
system%particle(j)%rolling_resistance_moment = system%particle(j)%rolling_resistance_moment - mirol

END SUBROUTINE sphere_sphere_contact_by_hertz_with_thornton
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_contact_by_hertz_with_thornton(i,k,overlap,vi,vw,wi,ww,nout,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, Gi, ri, mi, Eeq, Geq, req, meq, misi, misw, mis, midi, midw, mid, deltat, xsini, xsinw, xsin, &
                      xsiti, xsitw, xsit, vin, vwn, vreln, norm_vrelctang, kn, dn, norm_fcon, kt, dt, norm_ffri_trial, dfcon, norm_fconr, overlap_r, ktr
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, vrelctang, tiw, rci, vic, vrelc, fspring, fspring_trial, dxiw_trial, &
                                  fspringr, fspringr_proj, ffri_trial, tiw_trial, fconr, omegarel, mirol
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8, pi=4.0_real8*ATAN(1.0_real8)                     


! get particle data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle and wall contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
xsinw = system%wall(k)%contact_damping_ratio
xsitw = system%wall(k)%friction_damping_ratio
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! get wall k address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))

! compute contact parameters (average) for the pair i-k
xsin = half*(xsini+xsinw)
xsit = half*(xsiti+xsitw)
mis = half*(misi+misw)
mid = half*(midi+midw)

! compute equivalent properties for the pair i-k
Eeq = Ei/(one-ni_i*ni_i)
Geq = Gi/(two-ni_i)
req = ri
meq = mi

! compute normal contact parameters (spring stiffness and damping constant)
kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal velocities
vin = DOT_PRODUCT(vi,nout)
vwn = DOT_PRODUCT(vw,nout)
vreln = vin - vwn

! compute normal contact force (includes contact damping)
fcon = -1.3333333333333333_real8*Eeq*SQRT(req)*(overlap**1.5_real8)*nout - (dn*vreln)*nout 
IF (DOT_PRODUCT(fcon,nout)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! update particle´s total contact force
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
!system%wall(k)%contact_force = system%wall(k)%contact_force - fcon

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! build particle´s contact point radial vector
rci = ri*nout
!ccccccccccccccccccccccccc
!rci = (ri-overlap)*nout
!ccccccccccccccccccccccccc

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
!vwc = vw + (ww .vector. rcw)
vrelc = vic - vw

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout
norm_vrelctang = SQRT(DOT_PRODUCT(vrelctang,vrelctang))

! compute trial incremental elongation of friction spring
dxiw_trial = vrelctang*deltat

! compute friction parameters (spring stiffness and damping constant)
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
dt = 2.0_real8*xsit*SQRT(meq*kt)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! compute change on normal contact force
fconr = system%particle(i)%contacting_walls(waddr)%normalspring_force
norm_fconr = SQRT(DOT_PRODUCT(fconr,fconr))
dfcon = norm_fcon - norm_fconr
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute trial friction spring force
fspringr = system%particle(i)%contacting_walls(waddr)%fricspring_force
fspringr_proj = fspringr - DOT_PRODUCT(fspringr,nout)*nout
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
IF (dfcon>=0.0_real8) THEN
  fspring_trial = fspringr_proj + kt*dxiw_trial
ELSE
  overlap_r = system%particle(i)%contacting_walls(waddr)%normalspring_elongation
  ktr = 8.0_real8*Geq*SQRT(req)*SQRT(overlap_r)
  fspring_trial = fspringr_proj*(kt/ktr) + kt*dxiw_trial
END IF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute trial friction force
ffri_trial = - fspring_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tiw_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tiw_trial = 0.0_real8
END IF

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 !xiw = xiw_trial
 fspring = fspring_trial
 tiw = tiw_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tiw = tiw_trial
  ffri = mid*norm_fcon*tiw
  fspring = - (ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update particle´s total friction force and pair´s friction spring force
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
!system%wall(k)%friction_force = system%wall(k)%friction_force - ffri
!system%particle(i)%contacting_walls(waddr)%fricspring_elongation = xiw
system%particle(i)%contacting_walls(waddr)%fricspring_force = fspring

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! update pair´s normal spring elongation and contact force
system%particle(i)%contacting_walls(waddr)%normalspring_elongation = overlap
system%particle(i)%contacting_walls(waddr)%normalspring_force = fcon
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! compute moment due to fritcion force
mifri = rci .vector. ffri

! update particle´s total friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri

! compute moment due to rolling resistance
CALL compute_sphere_wall_rolling_resistance_moment(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moment
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol

END SUBROUTINE sphere_wall_contact_by_hertz_with_thornton
!!---------------------------------------------------------------------------------------------------------------------


END MODULE hertz_contact_with_thornton_friction

!!=============================================================================================================================
!!=============================================================================================================================


MODULE linear_contact_with_stickslip_friction

USE particles_classes
USE array_of_vectors_class
USE rolling_resistance_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_sphere_contact_by_linearspr_with_stickslip(i,j,overlap,vi,vj,wi,wj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, Gi, Gj, ri, rj, mi, mj, Eeq, Geq, req, meq, deltat, xsini, xsinj, xsin, &
                      xsiti, xsitj, xsit, misi, misj, mis, midi, midj, mid, kn, dn, vin, vjn, vreln, &
                      norm_fcon, kt, dt, norm_ffri_trial
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, mjfri, vrelctang, tij, rci, rcj, vic, vjc, vrelc, &
                                  xij, xij_trial, dxij_trial, ffri_trial, tij_trial, xijr, xijr_proj, mirol
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
xsinj = system%material_set(jmat)%contact_damping_ratio
xsitj = system%material_set(jmat)%friction_damping_ratio
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! get particle j address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute contact parameters for the contacting pair (average)
xsin = half*(xsini+xsinj)
xsit = half*(xsiti+xsitj)
mis = half*(misi+misj)
mid = half*(midi+midj)

! compute normal contact spring stiffness and dashpot constant
kn = 0.1333333333333333_real8*Eeq*(req**0.25_real8) !based on an average overlap of 0.01Req introduced into Hertz formula
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! compute normal contact force
fcon = - kn*overlap*nij - dn*vreln*nij
IF (DOT_PRODUCT(fcon,nij)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! udpate particles´ total contact forces
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%particle(j)%contact_force = system%particle(j)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = (4.0_real8*Geq/Eeq)*kn  !based on eq. 40 from Thornton, Cummins and Cleary (2013)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particles´ contact point radial vectors
rci = (ri)*nij
rcj = -(rj)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nij
!rcj = -(rj-overlap)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
vjc = vj + (wj .vector. rcj)
vrelc = vic - vjc

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij

! compute trial incremental elongation of friction spring
dxij_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xijr = system%particle(i)%contacting_particles(jaddr)%fricspring_elongation
xijr_proj = xijr - DOT_PRODUCT(xijr,nij)*nij
xij_trial = xijr_proj + dxij_trial

! compute trial friction force
ffri_trial = - kt*xij_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tij_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tij_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xij = xij_trial
 tij = tij_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tij = tij_trial
  ffri = mid*norm_fcon*tij
  xij = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_particles(jaddr)%fricspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%fricspring_elongation = xij

! update particles´ total friction forces
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%particle(j)%friction_force = system%particle(j)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri
mjfri = rcj .vector. (-ffri)

! update particles´ total friction moments
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfri

! compute moment due to rolling resistance
CALL compute_sphere_sphere_rolling_resistance_moment(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moments
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol
system%particle(j)%rolling_resistance_moment = system%particle(j)%rolling_resistance_moment - mirol

END SUBROUTINE sphere_sphere_contact_by_linearspr_with_stickslip
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sphere_wall_contact_by_linearspr_with_stickslip(i,k,overlap,vi,vw,wi,ww,nout,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, Gi, ri, mi, Eeq, Geq, req, meq, misi, misw, mis, midi, midw, mid, deltat, xsini, xsinw, xsin, &
                      xsiti, xsitw, xsit, vin, vwn, vreln, kn, dn, norm_fcon, kt, dt, norm_ffri_trial
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, vrelctang, tiw, rci, vic, vrelc, xiw, xiw_trial, dxiw_trial, &
                                  xiwr, xiwr_proj, ffri_trial, tiw_trial, omegarel, mirol
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8                     


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s and wall´s contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
xsinw = system%wall(k)%contact_damping_ratio
xsitw = system%wall(k)%friction_damping_ratio
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! get wall k address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))

! compute equivalent properties for the contacting pair
Eeq = Ei/(one-ni_i*ni_i)
Geq = Gi/(two-ni_i)
req = ri
meq = mi

! compute contact parameters for the contacting pair (average)
xsin = half*(xsini+xsinw)
xsit = half*(xsiti+xsitw)
mis = half*(misi+misw)
mid = half*(midi+midw)

! compute normal contact spring stiffness and dashpot constant
kn = 0.1333333333333333_real8*Eeq*(req**0.25_real8) !based on an average overlap of 0.01Req introduced into Hertz formula
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nout)
vwn = DOT_PRODUCT(vw,nout)
vreln = vin - vwn

! compute normal contact force
fcon = - kn*overlap*nout -  dn*vreln*nout 
IF (DOT_PRODUCT(fcon,nout)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! update particle´s total contact force
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
!system%wall(k)%contact_force = system%wall(k)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = (4.0_real8*Geq/Eeq)*kn  !based on eq. 40 from Thornton, Cummins and Cleary (2013)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particle´s contact point radial vector
rci = ri*nout
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nout
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
!vwc = vw + (ww .vector. rcw)
vrelc = vic - vw

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout

! compute trial incremental elongation of friction spring
dxiw_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xiwr = system%particle(i)%contacting_walls(waddr)%fricspring_elongation
xiwr_proj = xiwr - DOT_PRODUCT(xiwr,nout)*nout
xiw_trial = xiwr_proj + dxiw_trial

! compute trial friction force
ffri_trial = - kt*xiw_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tiw_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tiw_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xiw = xiw_trial
 tiw = tiw_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tiw = tiw_trial
  ffri = mid*norm_fcon*tiw
  xiw = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_walls(waddr)%fricspring_elongation = xiw

! update particle´s total friction force
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
!system%wall(k)%friction_force = system%wall(k)%friction_force - ffri

! compute moment due to fritcion force
mifri = rci .vector. ffri

! update particle´s total friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri

! compute moment due to rolling resistance
CALL compute_sphere_wall_rolling_resistance_moment(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moment
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol

END SUBROUTINE sphere_wall_contact_by_linearspr_with_stickslip
!!---------------------------------------------------------------------------------------------------------------------


END MODULE linear_contact_with_stickslip_friction

!!=============================================================================================================================
!!=============================================================================================================================


MODULE hertz_contact_with_stsl_friction_and_adhesion

USE particles_classes
USE array_of_vectors_class
USE adhesion_class
USE rolling_resistance_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_sph_cont_by_hertz_w_ss_and_adh(i,j,overlap,vi,vj,wi,wj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, Gi, Gj, ri, rj, mi, mj, Eeq, Geq, req, meq, deltat, xsini, xsinj, xsin, &
                      xsiti, xsitj, xsit, misi, misj, mis, midi, midj, mid, kn, dn, dnnet, vin, vjn, vreln, &
                      norm_fcon, kt, dt, norm_ffri_trial, khertz, kjkr, A, Wadh, equiv_overlap
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, mjfri, vrelctang, tij, rci, rcj, vic, vjc, vrelc, &
                                  xij, xij_trial, dxij_trial, ffri_trial, tij_trial, xijr, xijr_proj, mirol, fadh, fd, vreln_vector, frep, fattr
REAL(KIND=real8), PARAMETER :: one=1.0_real8, half=0.5_real8, two=2.0_real8, fourthirds=1.3333333333333333_real8, twothirds=2.0_real8/3.0_real8, pi=4.0_real8*ATAN(1.0_real8)                     


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
!xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
!xsinj = system%material_set(jmat)%contact_damping_ratio
xsitj = system%material_set(jmat)%friction_damping_ratio
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! get particle j address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute contact parameters for the contacting pair (average)
!xsin = half*(xsini+xsinj)
xsit = half*(xsiti+xsitj)
mis = half*(misi+misj)
mid = half*(midi+midj)

! compute normal contact spring stiffness and dashpot constant
!kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
!dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! compute normal contact force
!fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nij !- dn*vreln*nij
!IF (DOT_PRODUCT(fcon,nij)>0.0_real8) THEN
!  fcon = 0.0_real8
!END IF

!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
! compute equivalent overlap at zero geometric overlap ("apparent overlap" when particles are tangent to each other)
Wadh = system%adhesion_model%model_parameters(1)
equiv_overlap = (1.0_real8/Req)*(two*pi*Req*Req*Wadh/Eeq)**twothirds

! compute total overlap
overlap = overlap + equiv_overlap
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

! compute adhesion effects (if present)
IF (system%general_cv%adhesion_switch=="on") THEN
  CALL compute_sphere_sphere_adhesion_force_and_moment(i,j,overlap,nij,vi,vj,wi,wj,system,fadh,dnnet,khertz,kjkr)
!ELSE
!  fadh = 0.0_real8
END IF

! compute net contact force (includes damping contribution) (think: maybe the IF below should be placed within the adhesion subroutines)
IF (dnnet<0.0_real8) THEN
  !dnnet = 0.0_real8
END IF
!ssssssssssss saturday 19/05/2018 sssssssssssssssssssss
fd = - dnnet*vreln*nij
vreln_vector = vreln*nij
!IF (DOT_PRODUCT(fd,vreln_vector)>0.0_real8) THEN
!  dnnet = 0.0_real8  !thins means fd=0 on fcon below)
!END IF
!!dnnet = -dnnet
A = system%adhesion_model%model_parameters(2)  !!Wadh = system%adhesion_model%model_parameters(1)
frep = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nij - A*vreln*khertz*nij
fattr = fadh - A*vreln*kjkr*nij  !+ A*vreln*kjkr*nij
IF (DOT_PRODUCT(frep,nij)>0.0_real8) THEN
  frep = 0.0_real8
END IF
IF (DOT_PRODUCT(fattr,nij)<0.0_real8) THEN
  fattr = 0.0_real8
  !fattr = fadh
END IF
!ssssssssssssssssssssssssssssssssssssssssssssssssssssss
fcon = frep + fattr
!fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nij + fadh - dnnet*vreln*nij
!IF (dnnet<0.0_real8) THEN
!  fcon = 0.0_real8
!  fadh = 0.0_real8
!END IF

! udpate particles´ total adhesion forces
system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
system%particle(j)%adhesion_force = system%particle(j)%adhesion_force - fadh
!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

! udpate particles´ total contact forces
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%particle(j)%contact_force = system%particle(j)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 100.0_real8*kt
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particles´ contact point radial vectors
rci = (ri)*nij
rcj = -(rj)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nij
!rcj = -(rj-overlap)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
vjc = vj + (wj .vector. rcj)
vrelc = vic - vjc

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij

! compute trial incremental elongation of friction spring
dxij_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xijr = system%particle(i)%contacting_particles(jaddr)%fricspring_elongation
xijr_proj = xijr - DOT_PRODUCT(xijr,nij)*nij
xij_trial = xijr_proj + dxij_trial

! compute trial friction force
ffri_trial = - kt*xij_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tij_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tij_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xij = xij_trial
 tij = tij_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tij = tij_trial
  ffri = mid*norm_fcon*tij
  xij = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_particles(jaddr)%fricspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%fricspring_elongation = xij

! update particles´ total friction forces
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%particle(j)%friction_force = system%particle(j)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri
mjfri = rcj .vector. (-ffri)

! update particles´ total friction moments
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfri

! compute moment due to rolling resistance
CALL compute_sphere_sphere_rolling_resistance_moment(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moments
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol
system%particle(j)%rolling_resistance_moment = system%particle(j)%rolling_resistance_moment - mirol

END SUBROUTINE sph_sph_cont_by_hertz_w_ss_and_adh
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_wall_cont_by_hertz_w_ss_and_adh(i,k,overlap,vi,vw,wi,ww,nout,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, Gi, ri, mi, Eeq, Geq, req, meq, misi, misw, mis, midi, midw, mid, deltat, xsini, &
                      xsinw, xsin, xsiti, xsitw, xsit, vin, vwn, vreln, kn, dn, dnnet, norm_fcon, kt, dt, norm_ffri_trial
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, vrelctang, tiw, rci, vic, vrelc, xiw, xiw_trial, dxiw_trial, &
                                  xiwr, xiwr_proj, ffri_trial, tiw_trial, mirol, fadh
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8, fourthirds=1.3333333333333333_real8                    


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s and wall´s contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
xsinw = system%wall(k)%contact_damping_ratio
xsitw = system%wall(k)%friction_damping_ratio
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! get wall k address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))

! compute equivalent properties for the contacting pair
Eeq = Ei/(one-ni_i*ni_i)
Geq = Gi/(two-ni_i)
req = ri
meq = mi

! compute contact parameters for the contacting pair (average) 
!xsin = half*(xsini+xsinw)
xsit = half*(xsiti+xsitw)
mis = half*(misi+misw)
mid = half*(midi+midw)

! compute normal contact spring stiffness and dashpot constant
!kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
!dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nout)
vwn = DOT_PRODUCT(vw,nout)
vreln = vin - vwn

! compute normal contact force
!fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nout - dn*vreln*nout 
!IF (DOT_PRODUCT(fcon,nout)>0.0_real8) THEN
!  fcon = 0.0_real8
!END IF
!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
IF (system%general_cv%adhesion_switch=="on") THEN
  CALL compute_sphere_wall_adhesion_force_and_moment(i,k,overlap,nout,vi,vw,wi,ww,system,fadh,dnnet)
END IF

! compute net contact force (includes damping contribution)
!IF (dnnet<0.0_real8) dnnet = 0.0_real8
fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nout + fadh - dnnet*vreln*nout
!IF (dnnet<0.0_real8) THEN
  !fcon = 0.0_real8
  !fadh = 0.0_real8
!END IF

! udpate particle´s total adhesion force
system%particle(i)%adhesion_force = system%particle(i)%adhesion_force + fadh
!system%wall(k)%adhesion_force = system%wall(k)%adhesion_force - fadh
!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

! update particle´s total contact force
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
!system%wall(k)%contact_force = system%wall(k)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 100.0_real8*kt
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particle´s contact point radial vector
rci = ri*nout
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nout
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
!vwc = vw + (ww .vector. rcw)
vrelc = vic - vw

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout

! compute trial incremental elongation of friction spring
dxiw_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xiwr = system%particle(i)%contacting_walls(waddr)%fricspring_elongation
xiwr_proj = xiwr - DOT_PRODUCT(xiwr,nout)*nout
xiw_trial = xiwr_proj + dxiw_trial

! compute trial friction force
ffri_trial = - kt*xiw_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tiw_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tiw_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xiw = xiw_trial
 tiw = tiw_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tiw = tiw_trial
  ffri = mid*norm_fcon*tiw
  xiw = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_walls(waddr)%fricspring_elongation = xiw

! update particle´s total friction force
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
!system%wall(k)%friction_force = system%wall(k)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri

! update particle´s total friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri

! compute moment due to rolling resistance
CALL compute_sphere_wall_rolling_resistance_moment(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,system,mirol)

! update particles´ total rolling resistance moment
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol

END SUBROUTINE sph_wall_cont_by_hertz_w_ss_and_adh
!!---------------------------------------------------------------------------------------------------------------------


END MODULE hertz_contact_with_stsl_friction_and_adhesion

!!=============================================================================================================================
!!=============================================================================================================================
    
    
MODULE hertz_contact_with_stsl_friction_and_thermalcond

USE particles_classes
USE array_of_vectors_class
USE rolling_resistance_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_sph_cont_by_hertz_w_ss_and_thcond(i,j,overlap,vi,vj,wi,wj,ti,tj,norm_ximxj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, ti, tj, norm_ximxj
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, nij

! local variables
INTEGER(KIND=int4) :: s, imat, jmat, cp, jaddr, iaddr, ss
REAL(KIND=real8)   :: Ei, Ej, ni_i, ni_j, Gi, Gj, ri, rj, mi, mj, Eeq, Geq, req, meq, deltat, xsini, xsinj, xsin, &
                      xsiti, xsitj, xsit, misi, misj, mis, midi, midj, mid, kn, dn, vin, vjn, vreln, norm_fcon, &
                      kt, dt, norm_ffri_trial, ki, kj, kij, Li, Aijc, qpcond, tidegr, tjdegr
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, mjfri, vrelctang, tij, rci, rcj, vic, vjc, vrelc, &
                                  xij, xij_trial, dxij_trial, ffri_trial, tij_trial, xijr, xijr_proj, mirol, fadh
REAL(KIND=real8), PARAMETER :: one=1.0_real8, half=0.5_real8, two=2.0_real8, fourthirds=1.3333333333333333_real8, pi=4.0_real8*ATAN(1.0_real8)                     


! get particles´ data
ri = system%particle(i)%radius
rj = system%particle(j)%radius
mi = system%particle(i)%mass
mj = system%particle(j)%mass
imat = system%particle(i)%material_set_number
jmat = system%particle(j)%material_set_number

! get particles´ contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
Ej = system%material_set(jmat)%elasticity_modulus
ni_j = system%material_set(jmat)%poisson_coeff
xsinj = system%material_set(jmat)%contact_damping_ratio
xsitj = system%material_set(jmat)%friction_damping_ratio
misj = system%material_set(jmat)%static_friction_coeff
midj = system%material_set(jmat)%dynamic_friction_coeff

! get particles´ thermal properties
ki = system%material_set(imat)%thermal_conductivity
kj = system%material_set(jmat)%thermal_conductivity
tidegr = system%material_set(imat)%degrading_temperature
tjdegr = system%material_set(jmat)%degrading_temperature

! get particle j address at contact list of particle i
DO cp=1,system%particle(i)%no_contacting_particles
  IF (system%particle(i)%contacting_particles(cp)%part_number==j) jaddr = cp
END DO

! get particle i address at contact list of particle j
DO cp=1,system%particle(j)%no_contacting_particles
  IF (system%particle(j)%contacting_particles(cp)%part_number==i) iaddr = cp
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particles´ elasticity moduli (as a function of temperature)
Ei = min(Ei,Ei*exp(one-ti/tidegr))
Ej = min(Ej,Ej*exp(one-tj/tjdegr))

! compute particles´ tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))
Gj = Ej/(2.0_real8*(one+ni_j))

! compute equivalent properties for the contacting pair
Eeq = Ei*Ej/(Ei + Ej - Ei*ni_j*ni_j - Ej*ni_i*ni_i)
Geq = Gi*Gj/(two*Gi + two*Gj - Gi*ni_j - Gj*ni_i)
req = ri*rj/(ri+rj)
meq = mi*mj/(mi+mj)

! compute contact parameters for the contacting pair (average)
xsin = half*(xsini+xsinj)
xsit = half*(xsiti+xsitj)
mis = half*(misi+misj)
mid = half*(midi+midj)

! compute normal contact spring stiffness and dashpot constant
kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nij)
vjn = DOT_PRODUCT(vj,nij)
vreln = vin - vjn

! compute normal contact force
fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nij - dn*vreln*nij
IF (DOT_PRODUCT(fcon,nij)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! udpate particles´ total contact forces
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%particle(j)%contact_force = system%particle(j)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 100.0_real8*kt
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particles´ contact point radial vectors
rci = (ri)*nij
rcj = -(rj)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nij
!rcj = -(rj-overlap)*nij
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
vjc = vj + (wj .vector. rcj)
vrelc = vic - vjc

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nij))*nij

! compute trial incremental elongation of friction spring
dxij_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xijr = system%particle(i)%contacting_particles(jaddr)%fricspring_elongation
xijr_proj = xijr - DOT_PRODUCT(xijr,nij)*nij
xij_trial = xijr_proj + dxij_trial

! compute trial friction force
ffri_trial = - kt*xij_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tij_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tij_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xij = xij_trial
 tij = tij_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tij = tij_trial
  ffri = mid*norm_fcon*tij
  xij = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_particles(jaddr)%fricspring_elongation = xij
system%particle(j)%contacting_particles(iaddr)%fricspring_elongation = xij

! update particles´ total friction forces
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%particle(j)%friction_force = system%particle(j)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri
mjfri = rcj .vector. (-ffri)

! update particles´ total friction moments
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri
system%particle(j)%friction_moment = system%particle(j)%friction_moment + mjfri

! compute moment due to rolling resistance
CALL compute_sph_sph_roll_resist_mom_w_thcond(i,j,vi,vj,wi,wj,overlap,nij,norm_fcon,vrelctang,Ei,Ej,system,mirol)

! update particles´ total rolling resistance moments
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol
system%particle(j)%rolling_resistance_moment = system%particle(j)%rolling_resistance_moment - mirol

! compute particles´ effective conductivity
kij = (ri+rj)/(ri/ki+rj/kj)

! compute distance from particle i´s center to contact point
Li = half*(norm_ximxj-(rj*rj-ri*ri)/norm_ximxj)

! compute particles´ contact area (remark: contact area from Hertz contact theory is not consistent with the geometrical overlap)
Aijc = pi*(ri*ri-Li*Li)

! compute heat power due to conduction
qpcond = kij*(tj-ti)*Aijc/norm_ximxj

! update particles´ total conduction heat power
system%particle(i)%conduction_heat_power = system%particle(i)%conduction_heat_power + qpcond
system%particle(j)%conduction_heat_power = system%particle(j)%conduction_heat_power - qpcond

END SUBROUTINE sph_sph_cont_by_hertz_w_ss_and_thcond
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE sph_wall_cont_by_hertz_w_ss_and_thcond(i,k,overlap,vi,vw,wi,ww,ti,tw,diw,nout,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, ti, tw, diw
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, nout

! local variables
INTEGER(KIND=int4) :: s, imat, waddr, j
REAL(KIND=real8)   :: Ei, ni_i, Gi, ri, mi, Eeq, Geq, req, meq, misi, misw, mis, midi, midw, mid, deltat, xsini, &
                      xsinw, xsin, xsiti, xsitw, xsit, vin, vwn, vreln, kn, dn, norm_fcon, kt, dt, norm_ffri_trial, ki, kiw, Aiwc, qpcond, tidegr
REAL(KIND=real8), DIMENSION(3) :: fcon, ffri, mifri, vrelctang, tiw, rci, vic, vrelc, xiw, xiw_trial, dxiw_trial, &
                                  xiwr, xiwr_proj, ffri_trial, tiw_trial, mirol, fadh
REAL(KIND=real8), PARAMETER    :: one=1.0_real8, half=0.5_real8, two=2.0_real8, fourthirds=1.3333333333333333_real8, pi=4.0_real8*ATAN(1.0_real8)


! get particle´s data
mi = system%particle(i)%mass
ri = system%particle(i)%radius
imat = system%particle(i)%material_set_number

! get particle´s and wall´s contact-related material properties
Ei = system%material_set(imat)%elasticity_modulus
ni_i = system%material_set(imat)%poisson_coeff
xsini = system%material_set(imat)%contact_damping_ratio
xsiti = system%material_set(imat)%friction_damping_ratio
misi = system%material_set(imat)%static_friction_coeff
midi = system%material_set(imat)%dynamic_friction_coeff
xsinw = system%wall(k)%contact_damping_ratio
xsitw = system%wall(k)%friction_damping_ratio
misw = system%wall(k)%static_friction_coeff
midw = system%wall(k)%dynamic_friction_coeff

! get particle´s and wall´s thermal properties
ki = system%material_set(imat)%thermal_conductivity
tidegr = system%material_set(imat)%degrading_temperature
!kw = system%wall(k)%thermal_conductivity

! get wall k address at contact list of particle i
DO j=1,system%particle(i)%no_contacting_walls
  IF (system%particle(i)%contacting_walls(j)%wall_number==k) waddr = j
END DO

! get time step size
s = system%solution_cv%current_step
deltat = system%solution_cv%step(s)%current_dt

! compute particle´s elasticity moduli (as a function of temperature)
Ei = min(Ei,Ei*exp(one-ti/tidegr))

! compute particle´s tangential elastic modulus
Gi = Ei/(2.0_real8*(one+ni_i))

! compute equivalent properties for the contacting pair
Eeq = Ei/(one-ni_i*ni_i)
Geq = Gi/(two-ni_i)
req = ri
meq = mi

! compute contact parameters for the contacting pair (average) 
xsin = half*(xsini+xsinw)
xsit = half*(xsiti+xsitw)
mis = half*(misi+misw)
mid = half*(midi+midw)

! compute normal contact spring stiffness and dashpot constant
kn = 2.0_real8*Eeq*SQRT(req)*SQRT(overlap)
dn = 2.0_real8*xsin*SQRT(kn*meq)

! compute normal relative velocity
vin = DOT_PRODUCT(vi,nout)
vwn = DOT_PRODUCT(vw,nout)
vreln = vin - vwn

! compute normal contact force
fcon = - fourthirds*Eeq*SQRT(req)*(overlap**1.5_real8)*nout - dn*vreln*nout 
IF (DOT_PRODUCT(fcon,nout)>0.0_real8) THEN
  fcon = 0.0_real8
END IF

! update particle´s and wall´s total contact force
system%particle(i)%contact_force = system%particle(i)%contact_force + fcon
system%wall(k)%contact_force = system%wall(k)%contact_force - fcon

! compute friction spring stiffness and dashpot constant
kt = 8.0_real8*Geq*SQRT(req)*SQRT(overlap)
!xxxxxxxxxxxxxxxxxxxxxx
!kt = 100.0_real8*kt
!kt = 1000.0_real8*kt
!xxxxxxxxxxxxxxxxxxxxxx
dt = 2.0_real8*xsit*SQRT(meq*kt)

! build particle´s contact point radial vector
rci = ri*nout
!xxxxxxxxxxxxxxxxxxxxxxxx
!rci = (ri-overlap)*nout
!xxxxxxxxxxxxxxxxxxxxxxxx

! compute contact point velocity vectors
vic = vi + (wi .vector. rci)
!vwc = vw + (ww .vector. rcw)
vrelc = vic - vw

! compute tangential relative velocity at contact point
vrelctang = vrelc - (DOT_PRODUCT(vrelc,nout))*nout

! compute trial incremental elongation of friction spring
dxiw_trial = vrelctang*deltat

! compute trial total elongation of friction spring
xiwr = system%particle(i)%contacting_walls(waddr)%fricspring_elongation
xiwr_proj = xiwr - DOT_PRODUCT(xiwr,nout)*nout
xiw_trial = xiwr_proj + dxiw_trial

! compute trial friction force
ffri_trial = - kt*xiw_trial - dt*vrelctang

! compute trial friction force direction
norm_ffri_trial = SQRT(DOT_PRODUCT(ffri_trial,ffri_trial))
IF (norm_ffri_trial/=0.0_real8) THEN
  tiw_trial = (one/norm_ffri_trial)*ffri_trial
ELSE
  tiw_trial = 0.0_real8
END IF

! compute norm of normal contact force
norm_fcon = SQRT(DOT_PRODUCT(fcon,fcon))

! check for stick or slip
IF (norm_ffri_trial<=(mis*norm_fcon)) THEN
 ffri = ffri_trial
 xiw = xiw_trial
 tiw = tiw_trial
 !system%particle(i)%friction_status = "stick region"
ELSE
  tiw = tiw_trial
  ffri = mid*norm_fcon*tiw
  xiw = -(one/kt)*(ffri - dt*vrelctang)
  !system%particle(i)%friction_status = "slip region"
END IF

! update pair´s friction spring elongation
system%particle(i)%contacting_walls(waddr)%fricspring_elongation = xiw

! update particle´s and wall´s total friction force
system%particle(i)%friction_force = system%particle(i)%friction_force + ffri
system%wall(k)%friction_force = system%wall(k)%friction_force - ffri

! compute moment due to friction force
mifri = rci .vector. ffri

! update particle´s total friction moment
system%particle(i)%friction_moment = system%particle(i)%friction_moment + mifri

! compute moment due to rolling resistance
CALL compute_sph_wall_roll_resist_mom_w_thcond(i,k,vi,vw,wi,ww,overlap,nout,norm_fcon,vrelctang,Ei,system,mirol)

! update particle´ total rolling resistance moment
system%particle(i)%rolling_resistance_moment = system%particle(i)%rolling_resistance_moment + mirol

! compute equivalent conductivity for the contacting pair (assuming no conduction by the wall)
IF (system%wall(k)%thermal_kind=="thermally_active") THEN
  kiw = ki
ELSE IF (system%wall(k)%thermal_kind=="thermally_inactive") THEN
  kiw = 0.0_real8
ELSE
  WRITE(error_unit,*) " ERROR: wall thermal kind not recognized (check for possible mistyping)"
  STOP
END IF

! compute pair´s contact area (remark: contact area from Hertz contact theory is not consistent with the geometrical overlap)
Aiwc = pi*(ri*ri-diw*diw)

! compute heat power due to conduction
qpcond = kiw*(tw-ti)*Aiwc/diw

! update particle´s total conduction heat power
system%particle(i)%conduction_heat_power = system%particle(i)%conduction_heat_power + qpcond
!system%wall(k)%conduction_heat_power = system%wall(k)%conduction_heat_power - qpcond


END SUBROUTINE sph_wall_cont_by_hertz_w_ss_and_thcond
!!---------------------------------------------------------------------------------------------------------------------


END MODULE hertz_contact_with_stsl_friction_and_thermalcond

!!=============================================================================================================================
!!=============================================================================================================================


MODULE contact_models

USE impulse_based_contact
USE hertz_contact_with_contslide_friction
USE hertz_contact_with_stickslip_friction
USE hertz_contact_with_thornton_friction
USE linear_contact_with_stickslip_friction
USE hertz_contact_with_stsl_friction_and_adhesion
USE hertz_contact_with_stsl_friction_and_thermalcond
!USE hertz_contact_with_stsl_friction_and_adhesion_and_thcond

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE resolve_sphere_sphere_contact(i,j,overlap,vi,vj,wi,wj,aidelta,ajdelta,ti,tj,norm_ximxj,nij,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, j
REAL(KIND=real8)   :: overlap, ti, tj, norm_ximxj
REAL(KIND=real8), DIMENSION(3) :: vi, vj, wi, wj, aidelta, ajdelta, nij


! select contact model
SELECT CASE (system%solution_cv%contact_model)
  CASE ("impulse_based_contact")
    CALL sphere_sphere_contact_by_impulses(i,j,overlap,vi,vj,wi,wj,nij,system)
  CASE ("hertz_with_continuous_slide_friction")
    CALL sphere_sphere_contact_by_hertz_with_contslide(i,j,overlap,vi,vj,wi,wj,nij,system)
  CASE ("hertz_with_stickslip_friction")
    CALL sphere_sphere_contact_by_hertz_with_stickslip(i,j,overlap,vi,vj,wi,wj,nij,system)
  CASE ("hertz_with_thornton_friction")
    CALL sphere_sphere_contact_by_hertz_with_thornton(i,j,overlap,vi,vj,wi,wj,nij,system)
  CASE ("linear_with_stickslip_friction")
    CALL sphere_sphere_contact_by_linearspr_with_stickslip(i,j,overlap,vi,vj,wi,wj,nij,system)
  CASE ("hertz_with_stickslip_friction_and_adhesion")
    CALL sph_sph_cont_by_hertz_w_ss_and_adh(i,j,overlap,vi,vj,wi,wj,nij,system)
  CASE ("hertz_with_stickslip_friction_and_thermal_conduction")
    CALL sph_sph_cont_by_hertz_w_ss_and_thcond(i,j,overlap,vi,vj,wi,wj,ti,tj,norm_ximxj,nij,system)
  !CASE ("hertz_with_stickslip_friction_and_adhesion_and_thermal_conduction")
  !  CALL sph_sph_cont_by_hertz_w_ss_and_adh_and_thcond(i,j,overlap,vi,vj,wi,wj,norm_ximxj,nij,system)
  CASE DEFAULT
    WRITE(error_unit,*) " Contact model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE resolve_sphere_sphere_contact
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE resolve_sphere_wall_contact(i,k,overlap,vi,vw,wi,ww,ri,aidelta,ti,tw,diw,n_out,system)

! dymmy arguments
TYPE(system_data)  :: system
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8)   :: overlap, ti, tw, diw
REAL(KIND=real8), DIMENSION(3) :: vi, vw, wi, ww, ri, aidelta, n_out


! select contact model
SELECT CASE (system%solution_cv%contact_model)
  CASE ("impulse_based_contact")
    CALL sphere_wall_contact_by_impulses(i,k,overlap,vi,vw,wi,ww,n_out,system)
  CASE ("hertz_with_continuous_slide_friction")
    CALL sphere_wall_contact_by_hertz_with_contslide(i,k,overlap,vi,vw,wi,ww,n_out,system)
  CASE ("hertz_with_stickslip_friction")
    CALL sphere_wall_contact_by_hertz_with_stickslip(i,k,overlap,vi,vw,wi,ww,n_out,system)
  CASE ("hertz_with_thornton_friction")
    CALL sphere_wall_contact_by_hertz_with_thornton(i,k,overlap,vi,vw,wi,ww,n_out,system)
  CASE ("linear_with_stickslip_friction")
    CALL sphere_wall_contact_by_linearspr_with_stickslip(i,k,overlap,vi,vw,wi,ww,n_out,system)
  CASE ("hertz_with_stickslip_friction_and_adhesion")
    CALL sph_wall_cont_by_hertz_w_ss_and_adh(i,k,overlap,vi,vw,wi,ww,n_out,system)
  CASE ("hertz_with_stickslip_friction_and_thermal_conduction")
    CALL sph_wall_cont_by_hertz_w_ss_and_thcond(i,k,overlap,vi,vw,wi,ww,ti,tw,diw,n_out,system)
  !CASE ("hertz_with_stickslip_friction_and_adhesion_and_thermal_conduction")
  !  CALL sph_wall_cont_by_hertz_w_ss_and_adh_and_thcond(i,k,overlap,vi,vw,wi,ww,diw,n_out,system)
  CASE DEFAULT
    WRITE(error_unit,*) " Contact model not recognized or not yet implemented"
    STOP
END SELECT

END SUBROUTINE resolve_sphere_wall_contact
!!---------------------------------------------------------------------------------------------------------------------

! add here subroutines for other morphologies of contact (e.g. resolve_ellipse_ellipse_contact, etc)

END MODULE contact_models

!!=============================================================================================================================
!!=============================================================================================================================


MODULE spheres_contact_solver_class

USE particles_classes
USE array_of_vectors_class
USE particle_contact_lists_subroutines
USE wall_contact_lists_subroutines
USE contact_models
USE adhesion_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE contact_solver_for_spheres(i,R,V,A,W,T,system)

! dymmy arguments
INTEGER(KIND=int4) :: i
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: A, W, T
TYPE(system_data)  :: system


! select contact detection algorithm
SELECT CASE (system%solution_cv%contact_detection_algorithm)

  CASE ("loop_over_all_particles")
    IF (system%solution_cv%contact_model=="impulse_based_contact") THEN
      CALL detect_and_resolve_sphere_contacts_by_loap_for_hard_spheres(i,R,V,A,W,T,system)
    ELSE
      CALL detect_and_resolve_sphere_contacts_by_looping_over_all_part(i,R,V,A,W,T,system)
    END IF

  CASE ("verlet_list")
    IF (system%solution_cv%contact_model=="impulse_based_contact") THEN
      WRITE(error_unit,*) " Contact model not yet implemented for the selected contact detection algorithm"
      STOP
    ELSE
      CALL detect_and_resolve_sphere_contacts_by_verlet_lists(i,R,V,A,W,T,system)
    END IF

  CASE ("binning")
    IF (system%solution_cv%contact_model=="impulse_based_contact") THEN
      WRITE(error_unit,*) " Contact model not yet implemented for the selected contact detection algorithm"
      STOP
    ELSE
      CALL detect_and_resolve_sphere_contacts_by_binning(i,R,V,A,W,T,system)
    END IF

  CASE DEFAULT
    WRITE(error_unit,*) " Contact detection algorithm not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE contact_solver_for_spheres
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE detect_and_resolve_sphere_contacts_by_loap_for_hard_spheres(i,R,V,A,W,T,system)

! dymmy arguments
INTEGER(KIND=int4) :: i
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: A, W, T
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: j, s, k
REAL(KIND=real8)   :: norm_ri_m_rj, overlap, vin, vjn, normal_relative_vel, dt, diw, vwn, ti_tpdt, tj_tpdt, tw_tpdt
REAL(KIND=real8), DIMENSION(3) :: ri, rj, ri_m_rj, nij, vi, vj, wi_tpdt, wj_tpdt, r1, ri1, n_out, vw, ww_tpdt, aidelta, ajdelta


! check for contacts with other particles and, if occurring, compute contact forces
sphere_sphere_contact: DO j=i+1,system%no_particles
  
  ! compute overlap for pair i-j
  ri = system%particle(i)%position
  rj = system%particle(j)%position
  ri_m_rj = ri - rj
  norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
  overlap = system%particle(i)%radius + system%particle(j)%radius - norm_ri_m_rj

  ! compute contact´s normal direction
  nij = -(1.0_real8/norm_ri_m_rj)*ri_m_rj

  ! compute normal relative velocity for pair i-j
  vi = system%particle(i)%velocity
  vj = system%particle(j)%velocity
  vin = DOT_PRODUCT(vi,nij)
  vjn = DOT_PRODUCT(vj,nij)
  normal_relative_vel = vin - vjn
    
  ! get particles´ incremental rotations
  IF (PRESENT(A)) THEN
    aidelta = A(i)%r
    ajdelta = A(j)%r
  ELSE
    aidelta = 0.0_real8
    ajdelta = 0.0_real8
  END IF    

  ! get particles´ iterative spins (spins at t+dt)
  IF (PRESENT(W)) THEN
    wi_tpdt = W(i)%r
    wj_tpdt = W(j)%r
  ELSE
    wi_tpdt = 0.0_real8
    wj_tpdt = 0.0_real8
  END IF    

  ! get particles´ iterative temperatures (temperatures at t+dt)
  IF (PRESENT(T)) THEN
    ti_tpdt = T(i)%r(1)
    tj_tpdt = T(j)%r(1)
  ELSE
    ti_tpdt = 0.0_real8
    tj_tpdt = 0.0_real8
  END IF    

  ! compute contact forces for pair i-j
  IF (overlap>=0.0_real8 .AND. normal_relative_vel>0.0_real8) THEN
    CALL resolve_sphere_sphere_contact(i,j,overlap,vi,vj,wi_tpdt,wj_tpdt,aidelta,ajdelta,ti_tpdt,tj_tpdt,norm_ri_m_rj,nij,system)
  END IF
    
END DO sphere_sphere_contact

! check for contacts with walls and, if occurring, compute contact forces
sphere_wall_contact: DO k=1,system%no_rigid_walls

  ! compute distance between particle i and wall k
  ri = system%particle(i)%position
  r1 = system%wall(k)%point_position
  n_out = system%wall(k)%outside_normal
  ri1 = r1 - ri
  diw = DOT_PRODUCT(ri1,n_out)
  IF (diw<0.0_real8) WRITE(error_unit,*) " CAUTION: negative distance between particle", i, " and wall", k, " !"
  diw = ABS(diw)   ! CAUTION: check whether this may spoil something...

  ! compute overlap for pair i-k
  overlap = system%particle(i)%radius - diw

  ! compute normal relative velocity for pair i-k
  vi = system%particle(i)%velocity
  vw = system%wall(k)%velocity
  vin = DOT_PRODUCT(vi,n_out)
  vwn = DOT_PRODUCT(vw,n_out)
  normal_relative_vel = vin - vwn
    
  ! get incremental rotations for particle i
  IF (PRESENT(A)) THEN
    aidelta = A(i)%r
  ELSE
    aidelta = 0.0_real8
  END IF    

  ! get pair´s iterative spins (spins at t+dt)
  IF (PRESENT(W)) THEN
    wi_tpdt = W(i)%r
    ww_tpdt = 0.0_real8
  ELSE
    wi_tpdt = 0.0_real8
    ww_tpdt = 0.0_real8
  END IF    

  ! get pair´s iterative temperatures (temperatures at t+dt)
  IF (PRESENT(T)) THEN
    ti_tpdt = T(i)%r(1)
    tw_tpdt = system%wall(k)%temperature
  ELSE
    ti_tpdt = 0.0_real8
    tw_tpdt = 0.0_real8
  END IF    

  ! compute contact forces for pair i-k
  IF (overlap>=0.0_real8 .AND. normal_relative_vel>0.0_real8 ) THEN
    CALL resolve_sphere_wall_contact(i,k,overlap,vi,vw,wi_tpdt,ww_tpdt,ri,aidelta,ti_tpdt,tw_tpdt,diw,n_out,system)
  END IF

END DO sphere_wall_contact

END SUBROUTINE detect_and_resolve_sphere_contacts_by_loap_for_hard_spheres
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE detect_and_resolve_sphere_contacts_by_looping_over_all_part(i,R,V,A,W,T,system)

! dymmy arguments
INTEGER(KIND=int4) :: i
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: A, W, T
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: j, k, ncp, ncw
REAL(KIND=real8) :: norm_ri_m_rj, overlap, vin, vjn, normal_relative_vel, diw, vwn, time, &
                    xi_wall, xf_wall, yi_wall, yf_wall, zi_wall, zf_wall, ti, tj, tw
REAL(KIND=real8), DIMENSION(3) :: ri, rj, ri_m_rj, nij, vi, vj, wi, wj, r1, ri1, n_out, vw, ww, aidelta, ajdelta, xciw
LOGICAL :: contact_point_falls_within_wall_limits

! check for contacts with other particles and, if occurring, compute contact forces
sphere_sphere_contacts: DO j=i+1,system%no_particles

  ! compute overlap for pair i-j
  ri = R(i)%r 
  rj = R(j)%r
  ri_m_rj = ri - rj
  norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
  overlap = system%particle(i)%radius + system%particle(j)%radius - norm_ri_m_rj

  ! compute normal direction for pair i-j
  nij = -(1.0_real8/norm_ri_m_rj)*ri_m_rj

  ! get particles´ incremental rotations
  IF (PRESENT(A)) THEN
    aidelta = A(i)%r
    ajdelta = A(j)%r
  ELSE
    aidelta = 0.0_real8
    ajdelta = 0.0_real8
  END IF    

  ! get particles´ velocities and spins
  vi = V(i)%r
  vj = V(j)%r
  IF (PRESENT(W)) THEN
    wi = W(i)%r
    wj = W(j)%r
  ELSE
    wi = 0.0_real8
    wj = 0.0_real8
  END IF

  ! get particles´ temperatures
  IF (PRESENT(T)) THEN
    ti = T(i)%r(1)
    tj = T(j)%r(1)
  ELSE
    ti = 0.0_real8
    tj = 0.0_real8
  END IF

  ! check for contact release between pair i-j
  ncp = system%particle(i)%no_contacting_particles
  IF (ANY(system%particle(i)%contacting_particles(1:ncp)%part_number==j) .AND. overlap<=0.0_real8) THEN 
    CALL remove_particle_from_particle_contact_list(i,j,system)
  END IF

  ! compute contact forces for pair i-j
  IF (overlap>0.0_real8) THEN
    IF (overlap>=(0.5_real8*MIN(system%particle(i)%radius,system%particle(j)%radius))) THEN
      WRITE(error_unit,*) "WARNING: particles", i, " and", j, " overlap more than 50% of their radii!"
    END IF
    ncp = system%particle(i)%no_contacting_particles
    IF (ALL(system%particle(i)%contacting_particles(1:ncp)%part_number/=j)) THEN
      CALL add_particle_to_particle_contact_list(i,j,system)
    END IF
    CALL resolve_sphere_sphere_contact(i,j,overlap,vi,vj,wi,wj,aidelta,ajdelta,ti,tj,norm_ri_m_rj,nij,system)
  END IF

END DO sphere_sphere_contacts

! check for contacts with walls and, if occurring, compute contact forces
sphere_wall_contacts: DO k=1,system%no_rigid_walls

  ! compute distance between particle i and wall k
  ri = R(i)%r
  r1 = system%wall(k)%point_position
  n_out = system%wall(k)%outside_normal
  ri1 = r1 - ri
  diw = DOT_PRODUCT(ri1,n_out)
  IF (diw<0.0_real8) THEN
    WRITE(error_unit,*) " CAUTION: center of particle", i, " is beyond wall", k, "; check for possible problem"  !WRITE(error_unit,*) " CAUTION: negative distance between particle", i, " and wall", k, " !"
  END IF
  ! Think: maybe the warning above is not needed: only the warning for excessive overlap (further down) is needed...
  diw = ABS(diw)   ! CAUTION: check whether this may spoil something...

  ! compute overlap for pair i-k
  overlap = system%particle(i)%radius - diw

  ! get incremental rotations for particle i
  IF (PRESENT(A)) THEN
    aidelta = A(i)%r
  ELSE
    aidelta = 0.0_real8
  END IF    

  ! get velocities and spins for pair i-k
  vi = V(i)%r
  vw = system%wall(k)%velocity
  IF (PRESENT(W)) THEN
    wi = W(i)%r
    ww = 0.0_real8 !system%wall(k)%spin
  ELSE
    wi = 0.0_real8
    ww = 0.0_real8
  END IF

  ! get temperatures for pair i-k
  IF (PRESENT(T)) THEN
    ti = T(i)%r(1)
    tw = system%wall(k)%temperature
  ELSE
    ti = 0.0_real8
    tw = 0.0_real8
  END IF
  
  ! check for contact release between pair i-k
  ncw = system%particle(i)%no_contacting_walls
  IF (ANY(system%particle(i)%contacting_walls(1:ncw)%wall_number==k) .AND. overlap<=0.0_real8) THEN 
    CALL remove_wall_from_particle_contact_list(i,k,system)
  END IF

  ! compute contact forces for pair i-k
  IF (overlap>0.0_real8) THEN
    IF (overlap>=(0.5_real8*system%particle(i)%radius)) THEN
      WRITE(error_unit,*) "WARNING: particle", i, " and wall", k, " overlap more than 50% of the particle´s radius!"
    END IF    
    ncw = system%particle(i)%no_contacting_walls
    IF (ALL(system%particle(i)%contacting_walls(1:ncw)%wall_number/=k)) THEN
      CALL add_wall_to_particle_contact_list(i,k,system)
    END IF
    CALL resolve_sphere_wall_contact(i,k,overlap,vi,vw,wi,ww,ri,aidelta,ti,tw,diw,n_out,system)
  END IF

END DO sphere_wall_contacts

END SUBROUTINE detect_and_resolve_sphere_contacts_by_looping_over_all_part
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE detect_and_resolve_sphere_contacts_by_verlet_lists(i,R,V,A,W,T,system)

! dymmy arguments
INTEGER(KIND=int4) :: i
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: A, W, T
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: iv, j, k, ncn, ncp, ncw
REAL(KIND=real8)   :: norm_ri_m_rj, overlap, vin, vjn, normal_relative_vel, diw, vwn, time, &
                    xi_wall, xf_wall, yi_wall, yf_wall, zi_wall, zf_wall, ti, tj, tw
REAL(KIND=real8), DIMENSION(3) :: ri, rj, ri_m_rj, nij, vi, vj, wi, wj, r1, ri1, n_out, vw, ww, aidelta, ajdelta, xciw
LOGICAL :: contact_point_falls_within_wall_limits


! get number of close neighbors
ncn = system%particle(i)%no_particles_into_verlet_list

! check for contacts with close neighbors and, if occurring, compute contact forces
sphere_sphere_contacts: DO iv=1,ncn

  ! get close neighbor´s number
  j = system%particle(i)%verlet_list(iv)
  
  ! check for contacts only for j>i (to avoid double check)
  IF (j>i) THEN
  
    ! compute overlap for pair i-j
    ri = R(i)%r 
    rj = R(j)%r
    ri_m_rj = ri - rj
    norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
    overlap = system%particle(i)%radius + system%particle(j)%radius - norm_ri_m_rj

    ! compute normal direction for pair i-j
    nij = -(1.0_real8/norm_ri_m_rj)*ri_m_rj

    ! get particles´ incremental rotations
    IF (PRESENT(A)) THEN
      aidelta = A(i)%r
      ajdelta = A(j)%r
    ELSE
      aidelta = 0.0_real8
      ajdelta = 0.0_real8
    END IF    

    ! get particles´ velocities and spins
    vi = V(i)%r
    vj = V(j)%r
    IF (PRESENT(W)) THEN
      wi = W(i)%r
      wj = W(j)%r
    ELSE
      wi = 0.0_real8
      wj = 0.0_real8
    END IF
    
    ! get particles´ temperatures
    IF (PRESENT(T)) THEN
      ti = T(i)%r(1)
      tj = T(j)%r(1)
    ELSE
      ti = 0.0_real8
      tj = 0.0_real8
    END IF

    ! check for contact release between pair i-j
    ncp = system%particle(i)%no_contacting_particles
    IF (ANY(system%particle(i)%contacting_particles(1:ncp)%part_number==j) .AND. overlap<=0.0_real8) THEN 
      CALL remove_particle_from_particle_contact_list(i,j,system)
    END IF
    
    ! compute contact forces for pair i-j
    IF (overlap>0.0_real8) THEN
      IF (overlap>=(0.5_real8*MIN(system%particle(i)%radius,system%particle(j)%radius))) THEN
      WRITE(error_unit,*) "WARNING: particles", i, " and", j, " overlap more than 50% of the smaller radius!"
      END IF
      ncp = system%particle(i)%no_contacting_particles
      IF (ALL(system%particle(i)%contacting_particles(1:ncp)%part_number/=j)) THEN
        CALL add_particle_to_particle_contact_list(i,j,system)
      END IF
      CALL resolve_sphere_sphere_contact(i,j,overlap,vi,vj,wi,wj,aidelta,ajdelta,ti,tj,norm_ri_m_rj,nij,system)
    END IF  

  END IF

END DO sphere_sphere_contacts

! check for contacts with walls and, if occurring, compute contact forces
sphere_wall_contacts: DO k=1,system%no_rigid_walls

  ! compute distance between particle i and wall k
  ri = R(i)%r
  r1 = system%wall(k)%point_position
  n_out = system%wall(k)%outside_normal
  ri1 = r1 - ri
  diw = DOT_PRODUCT(ri1,n_out)
  !xxxxxxxxxxxxxxxxxxxxxxxxxx 15/08/2017: comment lines below if wall limits are to be enforced further down xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  IF (diw<0.0_real8) THEN
    WRITE(error_unit,*) " CAUTION: center of particle", i, " is beyond wall", k, "; check for possible problem"  !WRITE(error_unit,*) " CAUTION: negative distance between particle", i, " and wall", k, " !"
  END IF
  !Think: maybe the warning above is not needed: only the warning for excessive overlap (further down) is needed...
  !xxxxxxxxxxxxxxxxxxxxxxxxxx 15/08/2017 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  diw = ABS(diw)   ! CAUTION: check whether this may spoil something...

  ! compute overlap for pair i-k
  overlap = system%particle(i)%radius - diw

  ! get incremental rotations for particle i
  IF (PRESENT(A)) THEN
    aidelta = A(i)%r
  ELSE
    aidelta = 0.0_real8
  END IF    

  ! get velocities and spins for pair i-k
  vi = V(i)%r
  vw = system%wall(k)%velocity
  IF (PRESENT(W)) THEN
    wi = W(i)%r
    ww = 0.0_real8  !system%wall(k)%spin
  ELSE
    wi = 0.0_real8
    ww = 0.0_real8
  END IF    

  ! get temperatures for pair i-k
  IF (PRESENT(T)) THEN
    ti = T(i)%r(1)
    tw = system%wall(k)%temperature
  ELSE
    ti = 0.0_real8
    tw = 0.0_real8
  END IF
  
  !xxxxxxxxxxxxxx 15/08/2017 xxxxxxxxxxxxxxxxxxx define wall limits and check if contact point falls within xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ! define wall limits (Ai et al.´s example 7.5)
  !IF (k==1) THEN
  !  xi_wall = -0.3775_real8 !system%wall(k)%xi
  !  xf_wall = -0.0275_real8 !system%wall(k)%xf
  !  yi_wall = 0.57_real8 !system%wall(k)%yi
  !  yf_wall = 1.01_real8 !system%wall(k)%yf
  !  zi_wall = -0.1_real8 !system%wall(k)%zi
  !  zf_wall = 0.1_real8 !system%wall(k)%zf
  !ELSE IF (k==2) THEN
  !  xi_wall = -0.0280_real8
  !  xf_wall = -0.0275_real8
  !  yi_wall = 0.35_real8
  !  yf_wall = 0.57_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !ELSE IF (k==3) THEN
  !  xi_wall = 0.0275_real8
  !  xf_wall = 0.0280_real8
  !  yi_wall = 0.35_real8
  !  yf_wall = 0.57_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !ELSE IF (k==4) THEN
  !  xi_wall = 0.0275_real8
  !  xf_wall = 0.3775_real8
  !  yi_wall = 0.57_real8
  !  yf_wall = 1.01_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !ELSE IF (k==5) THEN
  !  xi_wall = -5.0_real8
  !  xf_wall = 5.0_real8
  !  yi_wall = -0.01_real8
  !  yf_wall = 0.0_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !ELSE IF (k==6) THEN   ! bottom of hopper
  !  xi_wall = -0.0275_real8
  !  xf_wall = 0.0275_real8
  !  yi_wall = 0.56_real8
  !  yf_wall = 0.57_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !ELSE IF (k==7) THEN   ! left lateral wall
  !  xi_wall = -0.3780_real8
  !  xf_wall = -0.3775_real8
  !  yi_wall = 1.01_real8
  !  yf_wall = 1.61_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !ELSE IF (k==8) THEN   ! right lateral wall
  !  xi_wall = 0.3775_real8
  !  xf_wall = 0.3780_real8
  !  yi_wall = 1.01_real8
  !  yf_wall = 1.61_real8
  !  zi_wall = -0.1_real8
  !  zf_wall = 0.1_real8
  !END IF
  !contact_point_falls_within_wall_limits = .FALSE.
  !IF (k==9 .AND. overlap>0.0_real8) contact_point_falls_within_wall_limits = .TRUE.
  !IF (overlap>0.0_real8) THEN
  !  ! compute position vector of contact point
  !  xciw = ri + diw*n_out
  !  ! check if contact point falls within wall limits
  !  IF (xciw(1)>=xi_wall .AND. xciw(1)<=xf_wall) THEN
  !     IF (xciw(2)>=yi_wall .AND. xciw(2)<=yf_wall) THEN
  !       IF (xciw(3)>=zi_wall .AND. xciw(3)<=zf_wall) THEN
  !         contact_point_falls_within_wall_limits = .TRUE.
  !       END IF
  !     END IF
  !  END IF
  !END IF
  !! check if particle is beyond wall (MAYBE THIS IS NOT NECESSARY... only the check for excessive overlap is needed...
  ! Note: the IF below needs to be corrected: it may not be over diw, since diw is positive here due to ABS(diw) further above
  !IF (diw<0.0_real8 .AND. contact_point_falls_within_wall_limits==.TRUE.) THEN
  !  WRITE(error_unit,*) " CAUTION: center of particle", i, " is beyond wall", k, "; check if this is ok"
  !END IF
  !diw = ABS(diw)   ! CAUTION: check whether this may spoil something...
   
  !! check for contact release between pair i-k
  !ncw = system%particle(i)%no_contacting_walls
  !!IF (ANY(system%particle(i)%contacting_walls(1:ncw)%wall_number==k) .AND. overlap<=0.0_real8) THEN 
  !!xxxxxx15/08Ixxxxxx
  !IF (ANY(system%particle(i)%contacting_walls(1:ncw)%wall_number==k) .AND. (overlap<=0.0_real8 .OR. contact_point_falls_within_wall_limits==.FALSE.)) THEN 
  !!xxxxxx15/08Fxxxxxx
  !  CALL remove_wall_from_particle_contact_list(i,k,system)
  !END IF

  !! compute contact forces for pair i-k
  !!IF (overlap>0.0_real8) THEN
  !!xxxxx15/08Ixxxxx
  !IF (overlap>0.0_real8 .AND. contact_point_falls_within_wall_limits==.TRUE.) THEN  
  !!xxxxx15/08Fxxxxx
  !  !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
  !  IF (overlap>=(0.5_real8*system%particle(i)%radius)) THEN
  !    WRITE(error_unit,*) "WARNING: particle", i, " and wall", k, " overlap more than 50% of the particle´s radius!"
  !  END IF    
  !  !wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
  !  ncw = system%particle(i)%no_contacting_walls
  !  IF (ALL(system%particle(i)%contacting_walls(1:ncw)%wall_number/=k)) THEN
  !    CALL add_wall_to_particle_contact_list(i,k,system)
  !  END IF
  !  CALL resolve_sphere_wall_contact(i,k,overlap,vi,vw,wi,ww,ri,aidelta,n_out,system)
  !  !IF (system%general_cv%adhesion_switch=="on") THEN
  !  !  CALL compute_sphere_wall_adhesion_force_and_moment(i,k,overlap,n_out,vi,vw,wi,ww,system)
  !  !END IF     
  !END IF
  !xxxxxxx 15/08/2017 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
  ! check for contact release between pair i-k
  ncw = system%particle(i)%no_contacting_walls
  IF (ANY(system%particle(i)%contacting_walls(1:ncw)%wall_number==k) .AND. overlap<=0.0_real8) THEN 
    CALL remove_wall_from_particle_contact_list(i,k,system)
    CALL remove_particle_from_wall_contact_list(k,i,system)
  END IF
  
  ! compute contact forces for pair i-k
  IF (overlap>0.0_real8) THEN
    IF (overlap>=(0.5_real8*system%particle(i)%radius)) THEN
      WRITE(error_unit,*) "WARNING: particle", i, " and wall", k, " overlap more than 50% of the particle´s radius!"
    END IF    
    ncw = system%particle(i)%no_contacting_walls
    IF (ALL(system%particle(i)%contacting_walls(1:ncw)%wall_number/=k)) THEN
      CALL add_wall_to_particle_contact_list(i,k,system)
      CALL add_particle_to_wall_contact_list(k,i,system)
    END IF
    CALL resolve_sphere_wall_contact(i,k,overlap,vi,vw,wi,ww,ri,aidelta,ti,tw,diw,n_out,system)
  END IF

END DO sphere_wall_contacts

END SUBROUTINE detect_and_resolve_sphere_contacts_by_verlet_lists
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE detect_and_resolve_sphere_contacts_by_binning(i,R,V,A,W,T,system)

! dymmy arguments
INTEGER(KIND=int4) :: i
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: A, W, T
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: cellx, celly, cellz, iix, iiy, iiz, no_neighbors_into_current_cell, ic, j, k, ncp, ncw
REAL(KIND=real8)   :: norm_ri_m_rj, overlap, vin, vjn, normal_relative_vel, diw, vwn, time, ti, tj, tw
REAL(KIND=real8), DIMENSION(3) :: ri, rj, ri_m_rj, nij, vi, vj, wi, wj, r1, ri1, n_out, vw, ww, aidelta, ajdelta


! get particle´s cell address
cellx = system%particle(i)%cell_address(1)
celly = system%particle(i)%cell_address(2)
cellz = system%particle(i)%cell_address(3)

! check for contacts with particles of neighboring bins and, if occurring, compute contact forces
sphere_sphere_contacts: DO iix=cellx-1,cellx+1
  DO iiy=celly-1,celly+1
    DO iiz=cellz-1,cellz+1
      
      ! get number of neighbors into current cell
      no_neighbors_into_current_cell = system%grid%cell(iix,iiy,iiz)%no_particles_into_cell
      
      ! loop over neighbors of current cell
      DO ic=1,no_neighbors_into_current_cell
        
        ! get neighbor´s number
        j = system%grid%cell(iix,iiy,iiz)%particle_list(ic)
  
        ! check for contacts only for j>i (to avoid double check)
        IF (j>i) THEN
       
          ! compute overlap for pair i-j
          ri = R(i)%r 
          rj = R(j)%r
          ri_m_rj = ri - rj
          norm_ri_m_rj = SQRT(DOT_PRODUCT(ri_m_rj,ri_m_rj))
          overlap = system%particle(i)%radius + system%particle(j)%radius - norm_ri_m_rj

          ! compute normal direction for pair i-j
          nij = -(1.0_real8/norm_ri_m_rj)*ri_m_rj

          ! get particles´ incremental rotations
          IF (PRESENT(A)) THEN
            aidelta = A(i)%r
            ajdelta = A(j)%r
          ELSE
            aidelta = 0.0_real8
            ajdelta = 0.0_real8
          END IF    

          ! get particles´ velocities and spins
          vi = V(i)%r
          vj = V(j)%r
          IF (PRESENT(W)) THEN
            wi = W(i)%r
            wj = W(j)%r
          ELSE
            wi = 0.0_real8
            wj = 0.0_real8
          END IF    
    
          ! get particles´ temperatures
          IF (PRESENT(T)) THEN
          ti = T(i)%r(1)
          tj = T(j)%r(1)
          ELSE
            ti = 0.0_real8
            tj = 0.0_real8
          END IF

          ! check for contact release between pair i-j
          ncp = system%particle(i)%no_contacting_particles
          IF (ANY(system%particle(i)%contacting_particles(1:ncp)%part_number==j) .AND. overlap<=0.0_real8) THEN 
            CALL remove_particle_from_particle_contact_list(i,j,system)
          END IF

          ! compute contact forces for pair i-j
          IF (overlap>0.0_real8) THEN
            IF (overlap>=(0.5_real8*MIN(system%particle(i)%radius,system%particle(j)%radius))) THEN
              WRITE(error_unit,*) "WARNING: particles", i, " and", j, " overlap more than 50% of the smaller radius!"
            END IF
            ncp = system%particle(i)%no_contacting_particles
            IF (ALL(system%particle(i)%contacting_particles(1:ncp)%part_number/=j)) THEN
              CALL add_particle_to_particle_contact_list(i,j,system)
            END IF
            CALL resolve_sphere_sphere_contact(i,j,overlap,vi,vj,wi,wj,aidelta,ajdelta,ti,tj,norm_ri_m_rj,nij,system)
          END IF
        
        END IF            
      
      END DO
    
    END DO
  END DO
END DO sphere_sphere_contacts

! check for contacts with walls and, if occurring, compute contact forces
sphere_wall_contacts: DO k=1,system%no_rigid_walls

  ! compute distance between particle i and wall k
  ri = R(i)%r
  r1 = system%wall(k)%point_position
  n_out = system%wall(k)%outside_normal
  ri1 = r1 - ri
  diw = DOT_PRODUCT(ri1,n_out)
  IF (diw<0.0_real8) THEN
    WRITE(error_unit,*) " CAUTION: center of particle", i, " is beyond wall", k, "; check for possible problem"  !WRITE(error_unit,*) " CAUTION: negative distance between particle", i, " and wall", k, " !"
  END IF
  diw = ABS(diw)   ! CAUTION: check whether this may spoil something...

  ! compute overlap for pair i-k
  overlap = system%particle(i)%radius - diw

  ! get incremental rotations for particle i
  IF (PRESENT(A)) THEN
    aidelta = A(i)%r
  ELSE
    aidelta = 0.0_real8
  END IF    

  ! get velocities and spins for pair i-k
  vi = V(i)%r
  vw = system%wall(k)%velocity
  IF (PRESENT(W)) THEN
    wi = W(i)%r
    ww = 0.0_real8  !system%wall(k)%spin
  ELSE
    wi = 0.0_real8
    ww = 0.0_real8
  END IF    

  ! get temperatures for pair i-k
  IF (PRESENT(T)) THEN
    ti = T(i)%r(1)
    tw = system%wall(k)%temperature
  ELSE
    ti = 0.0_real8
    tw = 0.0_real8
  END IF

  ! check for contact release between pair i-k
  ncw = system%particle(i)%no_contacting_walls
  IF (ANY(system%particle(i)%contacting_walls(1:ncw)%wall_number==k) .AND. overlap<=0.0_real8) THEN 
    CALL remove_wall_from_particle_contact_list(i,k,system)
  END IF

  ! compute contact forces for pair i-k
  IF (overlap>0.0_real8) THEN
    IF (overlap>=(0.5_real8*system%particle(i)%radius)) THEN
      WRITE(error_unit,*) "WARNING: particle", i, " and wall", k, " overlap more than 50% of the particle´s radius!"
    END IF    
    ncw = system%particle(i)%no_contacting_walls
    IF (ALL(system%particle(i)%contacting_walls(1:ncw)%wall_number/=k)) THEN
      CALL add_wall_to_particle_contact_list(i,k,system)
    END IF
    CALL resolve_sphere_wall_contact(i,k,overlap,vi,vw,wi,ww,ri,aidelta,ti,tw,diw,n_out,system)
  END IF

END DO sphere_wall_contacts

END SUBROUTINE detect_and_resolve_sphere_contacts_by_binning
!!---------------------------------------------------------------------------------------------------------------------


END MODULE spheres_contact_solver_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE contacts_solver

USE particles_classes
USE array_of_vectors_class
USE spheres_contact_solver_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE detect_and_resolve_particle_contacts(i,R,V,A,W,T,sys)

! dymmy arguments
INTEGER(KIND=int4) :: i
TYPE(array_of_vectors), DIMENSION(:), POINTER :: R, V
TYPE(array_of_vectors), DIMENSION(:), POINTER, OPTIONAL :: A, W, T
TYPE(system_data) :: sys


! select contact solver
SELECT CASE (sys%solution_cv%particles_shape)
  CASE ("spherical")
    CALL contact_solver_for_spheres(i,R,V,A,W,T,sys)
  CASE ("superellipsoidal")
    WRITE(error_unit,*) " Contact solver not yet implemented for this particle shape"
    STOP
  CASE ("general")
    WRITE(error_unit,*) " Contact solver not yet implemented for this particle shape"
    STOP
  CASE DEFAULT
    WRITE(error_unit,*) " Contact solver not yet implemented for this particle shape"
    STOP
END SELECT

END SUBROUTINE detect_and_resolve_particle_contacts
!!---------------------------------------------------------------------------------------------------------------------


END MODULE contacts_solver

!!=============================================================================================================================
!!=============================================================================================================================


MODULE dof_constraints

USE system_data_types

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE enforce_velocity_constraints(i,system,vtpdt)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: i
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(3), INTENT(INOUT) :: vtpdt

! local variables
INTEGER(KIND=int4) :: j

! loop over velocity dofs
DO j=4,6

  ! constant velocity
  IF (system%particle(i)%translational_dof_codes(j)==1 .AND. system%particle(i)%translational_dof_codes(j-3)==0) THEN
    vtpdt(j-3) = system%particle(i)%initial_velocity(j-3)
    !WRITE(error_unit,*) " Warning: velocity constraints were applied for particle = ", i, " dof = ", j
  
  ! other dof codes
  ELSE
    WRITE(error_unit,*) " Warning: velocity constraint type not yet implemented; check particle = ", i, " dof_code = ", j
    STOP

  END IF

END DO

END SUBROUTINE enforce_velocity_constraints
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE enforce_position_constraints(i,system,rtpdt,vtpdt)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: i
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(3), INTENT(INOUT) :: rtpdt, vtpdt

! local variables
INTEGER(KIND=int4) :: j, s
REAL(KIND=real8) :: a, w, t


! loop over position dofs
DO j=1,3

  ! fixed positions (no motion)
  IF (system%particle(i)%translational_dof_codes(j)==1) THEN
    rtpdt(j) = system%particle(i)%xyz_coordinates(j)
    vtpdt(j) = 0.0_real8
    !WRITE(error_unit,*) " Warning: position constraints were applied for particle = ", i, " dof = ", j
  
  ! harmonic motion
  ELSE IF (system%particle(i)%translational_dof_codes(j)==2) THEN
    a = system%particle(i)%harmonic_constraints_data(j)
    w = system%particle(i)%harmonic_constraints_data(j+3)
    s = system%solution_cv%current_step
    t = system%solution_cv%step(s)%current_time
    rtpdt(j) = system%particle(i)%xyz_coordinates(j) + a*sin(w*t)
    vtpdt(j) = a*w*cos(w*t)
    !WRITE(error_unit,*) " Warning: position constraints were applied for particle = ", i, " dof = ", j
  
  ! other dof codes
  ELSE
    WRITE(error_unit,*) " Warning: position constraint type not yet implemented; check particle = ", i, " dof_code = ", j
    STOP

  END IF

END DO

END SUBROUTINE enforce_position_constraints
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE enforce_spin_constraints(i,system,wtpdt)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: i
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(3), INTENT(INOUT) :: wtpdt

! local variables
INTEGER(KIND=int4) :: j

! loop over spin dofs
DO j=4,6

  ! constant spin
  IF (system%particle(i)%rotational_dof_codes(j)==1 .AND. system%particle(i)%rotational_dof_codes(j-3)==0) THEN
    wtpdt(j-3) = system%particle(i)%initial_spin(j-3)
    !WRITE(error_unit,*) " Warning: spin constraints were applied for particle = ", i, " dof = ", j

  ! other dof codes
  ELSE
    WRITE(error_unit,*) " Warning: spin constraint type not yet implemented; check particle = ", i, " dof_code = ", j
    STOP

  END IF

END DO

END SUBROUTINE enforce_spin_constraints
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE enforce_angle_constraints(i,system,adelta,wtpdt)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: i
TYPE(system_data) :: system
REAL(KIND=real8), DIMENSION(3), INTENT(INOUT) :: adelta, wtpdt

! local variables
INTEGER(KIND=int4) :: j

! loop over angle dofs
DO j=1,3

  ! fixed angles (no motion)
  IF (system%particle(i)%rotational_dof_codes(j)==1) THEN
    adelta(j) = 0.0_real8  !system%particle(i)%initial_angles(j)
    wtpdt(j) = 0.0_real8
    !WRITE(error_unit,*) " Warning: angle constraints were applied for particle = ", i, " dof = ", j

  ! harmonic angular motion
  ELSE IF (system%particle(i)%rotational_dof_codes(j)==2) THEN
    WRITE(error_unit,*) " Warning: harmonic angular constraints not yet implemented; check particle = ", i, " dof_code = ", j
    STOP

  ! other dof codes
  ELSE
    WRITE(error_unit,*) " Warning: angular constraint type not yet implemented; check particle = ", i, " dof_code = ", j
    STOP
  
  END IF

END DO

END SUBROUTINE enforce_angle_constraints
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE enforce_temperature_constraints(i,system,ttpdt)

! dummy arguments
INTEGER(KIND=int4), INTENT(IN) :: i
TYPE(system_data) :: system
REAL(KIND=real8), INTENT(INOUT) :: ttpdt

! local variables
INTEGER(KIND=int4) :: j, s
REAL(KIND=real8) :: a, w, t


! Note: currently there is only one temperature dof (the particles´ 13th dof), so there is no need to loop over temperature dofs

! fixed temperature (no heating or cooling)
IF (system%particle(i)%temperature_dof_codes==1) THEN
  ttpdt = system%particle(i)%initial_temperature
  !WRITE(error_unit,*) " Warning: temperature constraints were applied for particle = ", i

! harmonic heating or cooling
ELSE IF (system%particle(i)%temperature_dof_codes==2) THEN
    a = system%particle(i)%harmonic_constraints_data(1)
    w = system%particle(i)%harmonic_constraints_data(4)
    s = system%solution_cv%current_step
    t = system%solution_cv%step(s)%current_time
    ttpdt = system%particle(i)%initial_temperature + a*sin(w*t)
    !WRITE(error_unit,*) " Warning: temperature constraints were applied for particle = ", i
  
! other dof codes
ELSE
  WRITE(error_unit,*) " Warning: temperature constraint type not yet implemented; check particle = ", i, " dof_code(13)"
  STOP

END IF


END SUBROUTINE enforce_temperature_constraints
!!---------------------------------------------------------------------------------------------------------------------


END MODULE dof_constraints

!!=============================================================================================================================
!!=============================================================================================================================


MODULE system_properties

USE particles_classes
USE array_of_vectors_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_system_properties(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, k
REAL(KIND=real8) :: total_mass, mi, ji
REAL(KIND=real8), DIMENSION(3) :: ri, vi, wi

! clear system values
total_mass = 0.0_real8
system%center_of_mass = 0.0_real8
system%translational_kinetic_energy = 0.0_real8
system%rotational_kinetic_energy = 0.0_real8
system%linear_mom = 0.0_real8
system%angular_mom = 0.0_real8

! compute system energy, momenta and center of mass
DO i=1,system%no_particles
  mi = system%particle(i)%mass
  ji = system%particle(i)%inertia
  ri = system%particle(i)%position
  vi = system%particle(i)%velocity
  wi = system%particle(i)%spin
  total_mass = total_mass + mi
  system%center_of_mass = system%center_of_mass + mi*ri
  system%translational_kinetic_energy = system%translational_kinetic_energy + 0.5_real8*mi*DOT_PRODUCT(vi,vi)
  system%rotational_kinetic_energy = system%rotational_kinetic_energy + 0.5_real8*ji*DOT_PRODUCT(wi,wi)
  system%linear_mom = system%linear_mom + mi*vi
  system%angular_mom = system%angular_mom + ji*wi
END DO
system%center_of_mass = system%center_of_mass/total_mass
system%kinetic_energy = system%translational_kinetic_energy + system%rotational_kinetic_energy

! add here any other desired system properties (e.g. coordination number, fabric tensor, etc)


END SUBROUTINE compute_system_properties
!!---------------------------------------------------------------------------------------------------------------------


END MODULE system_properties

!!=============================================================================================================================
!!=============================================================================================================================


MODULE print_results

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_system_results(system)

! dummy arguments
TYPE(system_data) :: system

! print results according to desired file format
SELECT CASE (system%general_cv%results_file_format)
  CASE ("psy_format")
    CALL print_results_in_psy_format(system)
  CASE ("gid_format")
    CALL print_results_in_gid_format(system)
  CASE ("vtk_format")
    CALL print_vtk(system)
  CASE DEFAULT
    WRITE(error_unit,*) "results file format not recognized (please check caps: only small caps are allowed for the results file format type)"
    STOP
END SELECT

! print extra results in extra results file (if desired)
IF (system%general_cv%print_extra_results_file=="yes") CALL print_extra_results(system)


END SUBROUTINE print_system_results
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_results_in_psy_format(system)

! dummy argument
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, s, ss
REAL(KIND=real8)   :: time
REAL(KIND=real4)   :: time_real4, first_tprint_real4


! get step, substep and time info
s = system%solution_cv%current_step
ss = system%solution_cv%step(s)%current_substep
time = system%solution_cv%step(s)%current_time

! transform time values to single precision (to avoid round-off errors)
time_real4 = REAL(time,real4)
first_tprint_real4 = REAL(system%general_cv%dt_for_results_printing,real4)

! print system name and particles initial positions (only when printing results for the first time)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*) "System name: ", system%name
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*) "Particles initial positions (rx,ry,rz):"
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,"(6X,A,T21,A,T36,A,T51,A)") "Particle", "rx", "ry", "rz"
  DO i=1,system%no_particles
    WRITE(psy_results_unit,"(I10,T14,3ES15.6)") i,system%particle(i)%xyz_coordinates
  END DO
  WRITE(psy_results_unit,*)
END IF

! print step, substep and time info
WRITE(psy_results_unit,*) "==============================================================================="
WRITE(psy_results_unit,"(A,T9,I2,T16,A,T26,I10,T46,A,T53,ES15.6)") " Step =", s, "Substep =", ss, "Time =", time 
WRITE(psy_results_unit,*)

! print particles positions, velocities, angles and spins
IF (system%general_cv%rotational_dofs=="off" .OR. system%general_cv%print_rotational_dofs=="no") THEN
  WRITE(psy_results_unit,*) "Particles positions (rx,ry,rz) and velocities (vx,vy,vz):"
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,"(6X,A,T22,A,T37,A,T52,A,T67,A,T82,A,T97,A)") "Particle", "rx", "ry", "rz", "vx", "vy", "vz"
  DO i=1,system%no_particles
    WRITE(psy_results_unit,"(I10,T14,6ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity
  END DO
  WRITE(psy_results_unit,*)
ELSE
  WRITE(psy_results_unit,*) "Particles positions (r), velocities (v), angles (a) and spins (w):"
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,"(6X,A,T22,A,T37,A,T52,A,T67,A,T82,A,T97,A,T112,A,T127,A,T142,A,T157,A,T172,A,T187,A)") &
                      "Particle", "rx", "ry", "rz", "vx", "vy", "vz", "ax", "ay", "az", "wx", "wy", "wz"
  DO i=1,system%no_particles
    WRITE(psy_results_unit,"(I10,T14,12ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity, &
                                                 system%particle(i)%angle, system%particle(i)%spin
  END DO
  WRITE(psy_results_unit,*)
END IF

! print system properties (energy, momenta, center of mass, etc)
IF (system%general_cv%compute_system_properties=="yes") THEN
  IF (system%general_cv%rotational_dofs=="off") THEN
    system%translational_kinetic_energy = system%kinetic_energy
    system%rotational_kinetic_energy = 0.0_real8
  END IF
  WRITE(psy_results_unit,"(6X,A,T45,3ES15.6)") "System kinetic energy (tot,trans,rot):", system%kinetic_energy, &
                                                system%translational_kinetic_energy, system%rotational_kinetic_energy
  WRITE(psy_results_unit,"(6X,A,T45,3ES15.6)") "System linear momentum:", system%linear_mom
  WRITE(psy_results_unit,"(6X,A,T45,3ES15.6)") "System center of mass:", system%center_of_mass
  WRITE(psy_results_unit,*)
END IF

! add here any additional desired results (ex: particles´ forces, accelerations, etc)

END SUBROUTINE print_results_in_psy_format
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_results_in_gid_format(system)

! dummy argument
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, s, ss
REAL(KIND=real8)   :: time
REAL(KIND=real4)   :: time_real4, first_tprint_real4, vrelct


! get step, substep and time info
s = system%solution_cv%current_step
ss = system%solution_cv%step(s)%current_substep
time = system%solution_cv%step(s)%current_time

! transform time values to single precision (to avoid round-off errors)
time_real4 = REAL(time,real4)
first_tprint_real4 = REAL(system%general_cv%dt_for_results_printing,real4)

! print file header (only when printing results for the first time)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(GiD_results_unit,*) "GiD Post Results File 1.0"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*)
END IF

! print substep header and time info (as comment lines)
WRITE(GiD_results_unit,*) "#=============================================================================="
WRITE(GiD_results_unit,"(A,T9,I2,T16,A,T26,I10,T45,A,T52,ES15.6)") " # Step =", s, "Substep =", ss, "Time =", time 
WRITE(GiD_results_unit,*)

! print particles positions, velocities, angles and spins
IF (system%general_cv%rotational_dofs=="off" .OR. system%general_cv%print_rotational_dofs=="no") THEN
  WRITE(GiD_results_unit,"(A,A,A40,A,I10,A)") " ResultGroup", ' "', system%name, '" ', ss, "   OnNodes"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "positions"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "rx"', '  "ry"', '  "rz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "velocities"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "vx"', '  "vy"', '  "vz"'
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "Values"
  DO i=1,system%no_particles
    WRITE(GiD_results_unit,"(I10,T13,6ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity
  END DO
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
ELSE
  WRITE(GiD_results_unit,"(A,A,A40,A,I10,A)") " ResultGroup", ' "', system%name, '" ', ss, "   OnNodes"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "positions"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "rx"', '  "ry"', '  "rz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "velocities"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "vx"', '  "vy"', '  "vz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "angles"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "ax"', '  "ay"', '  "az"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "spins"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "wx"', '  "wy"', '  "wz"'
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "Values"
  DO i=1,system%no_particles
    WRITE(GiD_results_unit,"(I10,T13,12ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity, &
                                                   system%particle(i)%angle, system%particle(i)%spin 
  END DO
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
END IF

! print system kinetic energy
IF (system%general_cv%compute_system_properties=="yes") THEN
  IF (system%general_cv%rotational_dofs=="off") THEN
    WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "system kinetic energy"', ' "', system%name, '" ', ss, "   Scalar", "   OnNodes"
    WRITE(GiD_results_unit,*) "Values"
    WRITE(GiD_results_unit,"(I10,T13,ES15.6)") 1, system%kinetic_energy
    WRITE(GiD_results_unit,*) "End Values"
    WRITE(GiD_results_unit,*)
  ELSE
    WRITE(GiD_results_unit,"(A,A,A40,A,I10,A)") " ResultGroup", ' "', system%name, '" ', ss, "   OnNodes"
    WRITE(GiD_results_unit,*)
    WRITE(GiD_results_unit,*) "  ResultDescription", '  "system total kinetic energy"  ', "Scalar"
    WRITE(GiD_results_unit,*) "  ResultDescription", '  "system translational kinetic energy"  ', "Scalar"
    WRITE(GiD_results_unit,*) "  ResultDescription", '  "system rotational kinetic energy"  ', "Scalar"
    WRITE(GiD_results_unit,*)
    WRITE(GiD_results_unit,*) "Values"
    WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") 1, system%kinetic_energy, system%translational_kinetic_energy, system%rotational_kinetic_energy 
    WRITE(GiD_results_unit,*) "End Values"
    WRITE(GiD_results_unit,*)
  END IF
END IF

! print system linear momentum
IF (system%general_cv%compute_system_properties=="yes") THEN
  WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "system linear momentum"', ' "', system%name, '" ', ss, "   Vector", "   OnNodes"
  WRITE(GiD_results_unit,*) "Values"
  WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") 1, system%linear_mom
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
END IF

! print system center of mass
IF (system%general_cv%compute_system_properties=="yes") THEN
  WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "system center of mass"', ' "', system%name, '" ', ss, "   Vector", "   OnNodes"
  WRITE(GiD_results_unit,*) "Values"
  WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") 1, system%center_of_mass
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
END IF

! print post-process mesh file (particles´ numbers, radii and initial coordinates)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(GiD_mesh_unit,"(A,A,A40,A,A,A,A)") " MESH", ' "', system%name, '" ', " Dimension 3", "  ElemType Sphere", "  Nnode 1"  !CAUTION: it would be better if ElemType comes from the particle´s attributes or parameters...!
  WRITE(GiD_mesh_unit,*)
  WRITE(GiD_mesh_unit,*) " Coordinates"
  WRITE(GiD_mesh_unit,"(5X,A,TR6,A,TR13,A,TR13,A)") "#Particle", "rx", "ry", "rz"
  DO i=1,system%no_particles
    WRITE(GiD_mesh_unit,"(I10,T13,3ES15.6)") i, system%particle(i)%xyz_coordinates
  END DO
  WRITE(GiD_mesh_unit,*) " End Coordinates"
  WRITE(GiD_mesh_unit,*)
  WRITE(GiD_mesh_unit,*) " Elements"
  WRITE(GiD_mesh_unit,"(5X,A,TR5,A,TR7,A)") "#Particle", "node", "radius"
  DO i=1,system%no_particles
    WRITE(GiD_mesh_unit,"(I10,TR2,I10,TR3,ES15.6)") i, i, system%particle(i)%radius
  END DO
  WRITE(GiD_mesh_unit,*) " End Elements"
END IF


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! add below any additional desired results (ex: particles´ forces, accelerations, etc)

! print contact force for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(ES15.6,3ES15.6)") time, system%particle(1)%contact_force

! print contact and friction forces for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(A,T12,I5,T18,A,T25,ES15.6,3ES15.6,3ES15.6)") " Substep =", ss, "Time =", time, system%particle(1)%contact_force, system%particle(1)%friction_force

! print dof results and forces for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(A,TR3,A,ES15.6,TR3,A,15ES15.6)") " particle=1", "Time =", time, "xi,vi,wi,ficon,fifri=", system%particle(1)%position, system%particle(1)%velocity, system%particle(1)%spin, system%particle(1)%contact_force, system%particle(1)%friction_force


!below is the print done while fixing friction and rolling resistance (jan-feb-march 2017)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! print desired dof results and forces for a specific particle (print on error file for easiness)
!vrelct = system%particle(1)%velocity(1) + system%particle(1)%spin(3)*system%particle(1)%radius
!WRITE(error_unit,"(A,8ES15.6)") "t,xx,vx,wz,ffrix,mrolz,elong,vrelct=", time, system%particle(1)%position(1), system%particle(1)%velocity(1), system%particle(1)%spin(3), system%particle(1)%friction_force(1), system%particle(1)%rolling_resistance_moment(3), system%particle(1)%contacting_walls(1)%fricspring_elongation(1), vrelct
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


! print desired dof results and forces for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(A,TR3,A,ES15.6,TR3,A,3ES17.8)") " part=2", "t =", time, "x2x,v2x,w2z=", system%particle(2)%position(1), system%particle(2)%velocity(1), system%particle(2)%spin(3)


! print particles´ contact forces
!WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "particles contact force"', ' "', system%name, '" ', ss, "   Vector", "   OnNodes"
!WRITE(GiD_results_unit,*) "Values"
!DO i=1,system%no_particles
!  WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") i, system%particle(i)%contact_force
!END DO
!WRITE(GiD_results_unit,*) "End Values"
!WRITE(GiD_results_unit,*)

! print system kinetic energy (on error file for easiness)
!IF (system%general_cv%compute_system_properties=="yes") THEN
!  IF (system%general_cv%rotational_dofs=="off") THEN
!    system%translational_kinetic_energy = system%kinetic_energy
!    system%rotational_kinetic_energy = 0.0_real8
!  END IF
!  WRITE(error_unit,"(A,T10,ES15.6,T28,A,T66,3ES15.6)") " Time =", time, "System kinetic energy (tot,trans,rot):", system%kinetic_energy, system%translational_kinetic_energy, system%rotational_kinetic_energy
!END IF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END SUBROUTINE print_results_in_gid_format
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_vtk(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, unit_vtk=0
CHARACTER(len=256) :: character_buffer

! Henrique Gropelli´s code below
WRITE (character_buffer,'(A8,A1,I9.9,A)') input_file,'_',system%solution_cv%step(system%solution_cv%current_step)%current_substep, ".vtk"
OPEN(newunit=unit_vtk, access='stream', form='unformatted', FILE=character_buffer, STATUS='replace', convert='big_endian', BUFFERED='yes')
WRITE(unit_vtk) "# vtk DataFile Version 3.0",char(10),system%name,char(10),"BINARY",char(10),"DATASET UNSTRUCTURED_GRID",char(10)
WRITE(character_buffer,fmt="('POINTS',I,1X,'DOUBLE')")system%no_particles
IF (REAL(system%solution_cv%step(system%solution_cv%current_step)%current_time,real4)==REAL(system%general_cv%dt_for_results_printing,real4)) THEN
    WRITE(unit_vtk)character_buffer,char(10),(real(system%particle(i)%xyz_coordinates,real8),i=1,system%no_particles)
ELSE
    WRITE(unit_vtk)character_buffer,char(10),(real(system%particle(i)%position,real8),i=1,system%no_particles)
END IF
WRITE(unit_vtk)'FIELD FieldData 1',char(10),'Time 1 1 double',char(10),real(system%solution_cv%step(system%solution_cv%current_step)%current_time,real8)
WRITE(character_buffer,fmt="('POINT_DATA',I)")system%no_particles
WRITE(unit_vtk)character_buffer,char(10),'SCALARS Particle_radius double',char(10),'LOOKUP_TABLE default',char(10),(real(system%particle(i)%radius,real8),i=1,system%no_particles),'VECTORS Particle_velocity double',char(10),(real(system%particle(i)%velocity,real8),i=1,system%no_particles)
IF (system%general_cv%rotational_dofs=="on" .OR. system%general_cv%print_rotational_dofs=="yes") THEN
    WRITE(unit_vtk)'VECTORS Particle_spin double',char(10),(real(system%particle(i)%spin,real8),i=1,system%no_particles),'VECTORS Particle_angle double',char(10),(real(system%particle(i)%angle,real8),i=1,system%no_particles)
END IF
CLOSE(unit=unit_vtk)


END SUBROUTINE print_vtk
!!---------------------------------------------------------------------------------------------------------------------

 
SUBROUTINE print_extra_results(system)

! dummy argument
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, k, s, ss, no_infiltrated_particles, no_particles_that_crossed_midsection
REAL(KIND=real8)   :: time, vol_infiltrated_particles, vol_particles_that_crossed_midsection
REAL(KIND=real4)   :: time_real4, first_tprint_real4
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get step, substep and time info
s = system%solution_cv%current_step
ss = system%solution_cv%step(s)%current_substep
time = system%solution_cv%step(s)%current_time

! transform time values to single precision (to avoid round-off errors)
time_real4 = REAL(time,real4)
first_tprint_real4 = REAL(system%general_cv%dt_for_results_printing,real4)

! print system name (only at the time of first print)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(extra_results_unit,*)
  WRITE(extra_results_unit,*) "System name: ", system%name
  WRITE(extra_results_unit,*)
END IF

! print step, substep and time info
WRITE(extra_results_unit,*) "==============================================================================="
WRITE(extra_results_unit,"(A,T9,I2,T16,A,T26,I10,T46,A,T53,ES15.6)") " Step =", s, "Substep =", ss, "Time =", time 
WRITE(extra_results_unit,*)

! print extra results (add below any desired result upon demand)

! print walls´ contact-related quantities (number of contacting particles and total contact force)
WRITE(extra_results_unit,*) "Rigid walls contact-related quantities:" 
WRITE(extra_results_unit,"(3X,A,T10,A,T22,A,T37,A,T52,A,T67,A,T82,A,T97,A)") "Wall", "ncpart", "fconx", "fcony", "fconz", "ffricx", "ffricy", "ffricz"
DO k=1,system%no_rigid_walls
  WRITE(extra_results_unit,"(I6,T8,I6,T16,6ES15.6)") k, system%wall(k)%no_contacting_particles, system%wall(k)%contact_force, system%wall(k)%friction_force
END DO
WRITE(extra_results_unit,*)

! compute number and volume of infiltrated particles into porous solid (infiltration model, joint paper with Tarek, Feb/2018)
no_infiltrated_particles = 0
vol_infiltrated_particles = 0.0_real8
no_particles_that_crossed_midsection = 0
vol_particles_that_crossed_midsection = 0.0_real8
DO i=1,system%no_particles
  IF (system%particle(i)%position(1)>=0.0545_real8 .AND. system%particle(i)%parameters%kind=="granular_sphere_in_given_fluid_flow") THEN
    no_infiltrated_particles = no_infiltrated_particles + 1
    vol_infiltrated_particles = vol_infiltrated_particles + 1.333333333333333_real8*pi*system%particle(i)%radius**3
  END IF
  IF (system%particle(i)%position(1)>=0.1_real8 .AND. system%particle(i)%parameters%kind=="granular_sphere_in_given_fluid_flow") THEN
    no_particles_that_crossed_midsection = no_particles_that_crossed_midsection + 1
    vol_particles_that_crossed_midsection = vol_particles_that_crossed_midsection + 1.333333333333333_real8*pi*system%particle(i)%radius**3
  END IF
END DO

! print number and volume of particles that crossed porous solid´s front face (infiltration model, joint paper with Tarek, Feb/2018)
WRITE(extra_results_unit,*) "Infiltration results: number and volume of particles for which x>=0.0545m" 
WRITE(extra_results_unit,"(3X,A,T27,A)") "number of particles", "volume of particles"
WRITE(extra_results_unit,"(3X,I10,T27,ES15.6)") no_infiltrated_particles, vol_infiltrated_particles
WRITE(extra_results_unit,*)

! print number and volume of particles that crossed porous solid´s midsection (infiltration model, joint paper with Tarek, Feb/2018)
WRITE(extra_results_unit,*) "Infiltration results: number and volume of particles for which x>=0.1m" 
WRITE(extra_results_unit,"(3X,A,T27,A)") "number of particles", "volume of particles"
WRITE(extra_results_unit,"(3X,I10,T27,ES15.6)") no_particles_that_crossed_midsection, vol_particles_that_crossed_midsection
WRITE(extra_results_unit,*)

! print forced walls positions and velocities(testing the implementation of forced walls)
WRITE(extra_results_unit,*) "Forced walls positions and velocities" 
WRITE(extra_results_unit,"(3X,A,T16,A,T31,A,T46,A,T61,A,T76,A,T91,A)") "Wall", "x", "y", "z", "vx", "vy", "vz"
DO k=1,system%no_rigid_walls
  IF (system%wall(k)%kind=="forced_flat_rigid_wall") THEN
    WRITE(extra_results_unit,"(I6,T8,6ES15.6)") k, system%wall(k)%point_position, system%wall(k)%velocity
  END IF
END DO
WRITE(extra_results_unit,*)

! print system properties if this is more convenient here than in the .res file (energy, momenta, center of mass, etc)
!IF (system%general_cv%compute_system_properties=="yes") THEN
!  IF (system%general_cv%rotational_dofs=="off") THEN
!    system%translational_kinetic_energy = system%kinetic_energy
!    system%rotational_kinetic_energy = 0.0_real8
!  END IF
!  WRITE(extra_results_unit,"(6X,A,T45,3ES15.6)") "System kinetic energy (tot,trans,rot):", system%kinetic_energy, &
!                                                system%translational_kinetic_energy, system%rotational_kinetic_energy
!  WRITE(extra_results_unit,"(6X,A,T45,3ES15.6)") "System linear momentum:", system%linear_mom
!  WRITE(extra_results_unit,"(6X,A,T45,3ES15.6)") "System center of mass:", system%center_of_mass
!  WRITE(extra_results_unit,*)
!END IF

! add here any additional desired results (ex: particles´ forces, accelerations, etc)

END SUBROUTINE print_extra_results
!!---------------------------------------------------------------------------------------------------------------------


END MODULE print_results

!!=============================================================================================================================
!!=============================================================================================================================


MODULE print_results_for_thermomechanics_solver

USE particles_classes

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_system_thermomechanics_results(system)

! dummy arguments
TYPE(system_data) :: system

! print results according to desired file format
SELECT CASE (system%general_cv%results_file_format)
  CASE ("psy_format")
    CALL print_thermomechanics_results_in_psy_format(system)
  CASE ("gid_format")
    CALL print_thermomechanics_results_in_gid_format(system)
  CASE ("vtk_format")
    CALL print_vtk_temp(system)
  CASE DEFAULT
    WRITE(error_unit,*) "results file format not recognized (please check caps: only small caps are allowed for the results file format type)"
    STOP
END SELECT

! print extra results in extra results file (if desired)
IF (system%general_cv%print_extra_results_file=="yes") CALL print_thermomechanics_extra_results(system)


END SUBROUTINE print_system_thermomechanics_results
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_thermomechanics_results_in_psy_format(system)

! dummy argument
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, s, ss
REAL(KIND=real8)   :: time
REAL(KIND=real4)   :: time_real4, first_tprint_real4


! get step, substep and time info
s = system%solution_cv%current_step
ss = system%solution_cv%step(s)%current_substep
time = system%solution_cv%step(s)%current_time

! transform time values to single precision (to avoid round-off errors)
time_real4 = REAL(time,real4)
first_tprint_real4 = REAL(system%general_cv%dt_for_results_printing,real4)

! print system name and particles initial positions (only when printing results for the first time)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*) "System name: ", system%name
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,*) "Particles initial positions (rx,ry,rz):"
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,"(6X,A,T21,A,T36,A,T51,A)") "Particle", "rx", "ry", "rz"
  DO i=1,system%no_particles
    WRITE(psy_results_unit,"(I10,T14,3ES15.6)") i,system%particle(i)%xyz_coordinates
  END DO
  WRITE(psy_results_unit,*)
END IF

! print step, substep and time info
WRITE(psy_results_unit,*) "==============================================================================="
WRITE(psy_results_unit,"(A,T9,I2,T16,A,T26,I10,T46,A,T53,ES15.6)") " Step =", s, "Substep =", ss, "Time =", time 
WRITE(psy_results_unit,*)

! print particles positions, velocities, angles and spins
IF (system%general_cv%rotational_dofs=="off" .OR. system%general_cv%print_rotational_dofs=="no") THEN
  WRITE(psy_results_unit,*) "Particles positions (rx,ry,rz) and velocities (vx,vy,vz):"
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,"(6X,A,T22,A,T37,A,T52,A,T67,A,T82,A,T97,A)") "Particle", "rx", "ry", "rz", "vx", "vy", "vz"
  DO i=1,system%no_particles
    WRITE(psy_results_unit,"(I10,T14,6ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity
  END DO
  WRITE(psy_results_unit,*)
ELSE
  WRITE(psy_results_unit,*) "Particles positions (r), velocities (v), angles (a) and spins (w):"
  WRITE(psy_results_unit,*)
  WRITE(psy_results_unit,"(6X,A,T22,A,T37,A,T52,A,T67,A,T82,A,T97,A,T112,A,T127,A,T142,A,T157,A,T172,A,T187,A)") &
                      "Particle", "rx", "ry", "rz", "vx", "vy", "vz", "ax", "ay", "az", "wx", "wy", "wz"
  DO i=1,system%no_particles
    WRITE(psy_results_unit,"(I10,T14,12ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity, &
                                                 system%particle(i)%angle, system%particle(i)%spin
  END DO
  WRITE(psy_results_unit,*)
END IF

! print system properties (energy, momenta, center of mass, etc)
IF (system%general_cv%compute_system_properties=="yes") THEN
  IF (system%general_cv%rotational_dofs=="off") THEN
    system%translational_kinetic_energy = system%kinetic_energy
    system%rotational_kinetic_energy = 0.0_real8
  END IF
  WRITE(psy_results_unit,"(6X,A,T45,3ES15.6)") "System kinetic energy (tot,trans,rot):", system%kinetic_energy, &
                                                system%translational_kinetic_energy, system%rotational_kinetic_energy
  WRITE(psy_results_unit,"(6X,A,T45,3ES15.6)") "System linear momentum:", system%linear_mom
  WRITE(psy_results_unit,"(6X,A,T45,3ES15.6)") "System center of mass:", system%center_of_mass
  WRITE(psy_results_unit,*)
END IF

! add here any additional desired results (ex: particles´ forces, accelerations, etc)

END SUBROUTINE print_thermomechanics_results_in_psy_format
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_thermomechanics_results_in_gid_format(system)

! dummy argument
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, s, ss
REAL(KIND=real8)   :: time
REAL(KIND=real4)   :: time_real4, first_tprint_real4, vrelct


! get step, substep and time info
s = system%solution_cv%current_step
ss = system%solution_cv%step(s)%current_substep
time = system%solution_cv%step(s)%current_time

! transform time values to single precision (to avoid round-off errors)
time_real4 = REAL(time,real4)
first_tprint_real4 = REAL(system%general_cv%dt_for_results_printing,real4)

! print file header (only when printing results for the first time)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(GiD_results_unit,*) "GiD Post Results File 1.0"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*)
END IF

!-----------------------------------------------------
! print substep header and time info (as comment lines)
WRITE(GiD_results_unit,*) "#=============================================================================="
WRITE(GiD_results_unit,"(A,T9,I2,T16,A,T26,I10,T45,A,T52,ES15.6)") " # Step =", s, "Substep =", ss, "Time =", time 
WRITE(GiD_results_unit,*)

! print particles positions, velocities, angles and spins
IF (system%general_cv%rotational_dofs=="off" .OR. system%general_cv%print_rotational_dofs=="no") THEN
  WRITE(GiD_results_unit,"(A,A,A40,A,I10,A)") " ResultGroup", ' "', system%name, '" ', ss, "   OnNodes"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "positions"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "rx"', '  "ry"', '  "rz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "velocities"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "vx"', '  "vy"', '  "vz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "temperatures"  ', "Scalar"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "Values"
  DO i=1,system%no_particles
    WRITE(GiD_results_unit,"(I10,T13,7ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity, &
                                                   system%particle(i)%temperature
  END DO
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
ELSE
  WRITE(GiD_results_unit,"(A,A,A40,A,I10,A)") " ResultGroup", ' "', system%name, '" ', ss, "   OnNodes"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "positions"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "rx"', '  "ry"', '  "rz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "velocities"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "vx"', '  "vy"', '  "vz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "angles"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "ax"', '  "ay"', '  "az"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "spins"  ', "Vector"
  WRITE(GiD_results_unit,*) "    ComponentNames", '  "wx"', '  "wy"', '  "wz"'
  WRITE(GiD_results_unit,*) "  ResultDescription", '  "temperatures"  ', "Scalar"
  WRITE(GiD_results_unit,*)
  WRITE(GiD_results_unit,*) "Values"
  DO i=1,system%no_particles
    WRITE(GiD_results_unit,"(I10,T13,13ES15.6)") i, system%particle(i)%position, system%particle(i)%velocity, &
                                                   system%particle(i)%angle, system%particle(i)%spin, system%particle(i)%temperature
  END DO
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
END IF
!-----------------------------------------------------

!DO i=1,1
!    WRITE(GiD_results_unit,"(ES15.6,3ES15.6)") system%particle(i)%velocity(1)
!END DO


! print system kinetic energy
IF (system%general_cv%compute_system_properties=="yes") THEN
  IF (system%general_cv%rotational_dofs=="off") THEN
    WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "system kinetic energy"', ' "', system%name, '" ', ss, "   Scalar", "   OnNodes"
    WRITE(GiD_results_unit,*) "Values"
    WRITE(GiD_results_unit,"(I10,T13,ES15.6)") 1, system%kinetic_energy
    WRITE(GiD_results_unit,*) "End Values"
    WRITE(GiD_results_unit,*)
  ELSE
    WRITE(GiD_results_unit,"(A,A,A40,A,I10,A)") " ResultGroup", ' "', system%name, '" ', ss, "   OnNodes"
    WRITE(GiD_results_unit,*)
    WRITE(GiD_results_unit,*) "  ResultDescription", '  "system total kinetic energy"  ', "Scalar"
    WRITE(GiD_results_unit,*) "  ResultDescription", '  "system translational kinetic energy"  ', "Scalar"
    WRITE(GiD_results_unit,*) "  ResultDescription", '  "system rotational kinetic energy"  ', "Scalar"
    WRITE(GiD_results_unit,*)
    WRITE(GiD_results_unit,*) "Values"
    WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") 1, system%kinetic_energy, system%translational_kinetic_energy, system%rotational_kinetic_energy 
    WRITE(GiD_results_unit,*) "End Values"
    WRITE(GiD_results_unit,*)
  END IF
END IF

! print system linear momentum
IF (system%general_cv%compute_system_properties=="yes") THEN
  WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "system linear momentum"', ' "', system%name, '" ', ss, "   Vector", "   OnNodes"
  WRITE(GiD_results_unit,*) "Values"
  WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") 1, system%linear_mom
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
END IF

! print system center of mass
IF (system%general_cv%compute_system_properties=="yes") THEN
  WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "system center of mass"', ' "', system%name, '" ', ss, "   Vector", "   OnNodes"
  WRITE(GiD_results_unit,*) "Values"
  WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") 1, system%center_of_mass
  WRITE(GiD_results_unit,*) "End Values"
  WRITE(GiD_results_unit,*)
END IF

! print post-process mesh file (particles´ numbers, radii and initial coordinates)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(GiD_mesh_unit,"(A,A,A40,A,A,A,A)") " MESH", ' "', system%name, '" ', " Dimension 3", "  ElemType Sphere", "  Nnode 1"  !CAUTION: it would be better if ElemType comes from the particle´s attributes or parameters...!
  WRITE(GiD_mesh_unit,*)
  WRITE(GiD_mesh_unit,*) " Coordinates"
  WRITE(GiD_mesh_unit,"(5X,A,TR6,A,TR13,A,TR13,A)") "#Particle", "rx", "ry", "rz"
  DO i=1,system%no_particles
    WRITE(GiD_mesh_unit,"(I10,T13,3ES15.6)") i, system%particle(i)%xyz_coordinates
  END DO
  WRITE(GiD_mesh_unit,*) " End Coordinates"
  WRITE(GiD_mesh_unit,*)
  WRITE(GiD_mesh_unit,*) " Elements"
  WRITE(GiD_mesh_unit,"(5X,A,TR5,A,TR7,A)") "#Particle", "node", "radius"
  DO i=1,system%no_particles
    WRITE(GiD_mesh_unit,"(I10,TR2,I10,TR3,ES15.6)") i, i, system%particle(i)%radius
  END DO
  WRITE(GiD_mesh_unit,*) " End Elements"
END IF


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! add below any additional desired results (ex: particles´ forces, accelerations, etc)

! print contact force for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(ES15.6,3ES15.6)") time, system%particle(1)%contact_force

! print contact and friction forces for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(A,T12,I5,T18,A,T25,ES15.6,3ES15.6,3ES15.6)") " Substep =", ss, "Time =", time, system%particle(1)%contact_force, system%particle(1)%friction_force

! print dof results and forces for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(A,TR3,A,ES15.6,TR3,A,15ES15.6)") " particle=1", "Time =", time, "xi,vi,wi,ficon,fifri=", system%particle(1)%position, system%particle(1)%velocity, system%particle(1)%spin, system%particle(1)%contact_force, system%particle(1)%friction_force


!below is the print done while fixing friction and rolling resistance (jan-feb-march 2017)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! print desired dof results and forces for a specific particle (print on error file for easiness)
!vrelct = system%particle(1)%velocity(1) + system%particle(1)%spin(3)*system%particle(1)%radius
!WRITE(error_unit,"(A,8ES15.6)") "t,xx,vx,wz,ffrix,mrolz,elong,vrelct=", time, system%particle(1)%position(1), system%particle(1)%velocity(1), system%particle(1)%spin(3), system%particle(1)%friction_force(1), system%particle(1)%rolling_resistance_moment(3), system%particle(1)%contacting_walls(1)%fricspring_elongation(1), vrelct
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


! print desired dof results and forces for a specific particle (print on error file for easiness)
!WRITE(error_unit,"(A,TR3,A,ES15.6,TR3,A,3ES17.8)") " part=2", "t =", time, "x2x,v2x,w2z=", system%particle(2)%position(1), system%particle(2)%velocity(1), system%particle(2)%spin(3)


! print particles´ contact forces
!WRITE(GiD_results_unit,"(A,A,A,A40,A,I10,A,A)") " Result", '  "particles contact force"', ' "', system%name, '" ', ss, "   Vector", "   OnNodes"
!WRITE(GiD_results_unit,*) "Values"
!DO i=1,system%no_particles
!  WRITE(GiD_results_unit,"(I10,T13,3ES15.6)") i, system%particle(i)%contact_force
!END DO
!WRITE(GiD_results_unit,*) "End Values"
!WRITE(GiD_results_unit,*)

! print system kinetic energy (on error file for easiness)
!IF (system%general_cv%compute_system_properties=="yes") THEN
!  IF (system%general_cv%rotational_dofs=="off") THEN
!    system%translational_kinetic_energy = system%kinetic_energy
!    system%rotational_kinetic_energy = 0.0_real8
!  END IF
!  WRITE(error_unit,"(A,T10,ES15.6,T28,A,T66,3ES15.6)") " Time =", time, "System kinetic energy (tot,trans,rot):", system%kinetic_energy, system%translational_kinetic_energy, system%rotational_kinetic_energy
!END IF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END SUBROUTINE print_thermomechanics_results_in_gid_format
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE print_vtk_temp(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, unit_vtk=0
CHARACTER(len=256) :: character_buffer

! Henrique Gropelli´s code below
WRITE (character_buffer,'(A8,A1,I9.9,A)') input_file,'_',system%solution_cv%step(system%solution_cv%current_step)%current_substep, ".vtk"
OPEN(newunit=unit_vtk, access='stream', form='unformatted', FILE=character_buffer, STATUS='replace', convert='big_endian', BUFFERED='yes')
WRITE(unit_vtk) "# vtk DataFile Version 3.0",char(10),system%name,char(10),"BINARY",char(10),"DATASET UNSTRUCTURED_GRID",char(10)
WRITE(character_buffer,fmt="('POINTS',I,1X,'DOUBLE')")system%no_particles
IF (REAL(system%solution_cv%step(system%solution_cv%current_step)%current_time,real4)==REAL(system%general_cv%dt_for_results_printing,real4)) THEN
    WRITE(unit_vtk)character_buffer,char(10),(real(system%particle(i)%xyz_coordinates,real8),i=1,system%no_particles)
ELSE
    WRITE(unit_vtk)character_buffer,char(10),(real(system%particle(i)%position,real8),i=1,system%no_particles)
END IF
WRITE(unit_vtk)'FIELD FieldData 1',char(10),'Time 1 1 double',char(10),real(system%solution_cv%step(system%solution_cv%current_step)%current_time,real8)
WRITE(character_buffer,fmt="('POINT_DATA',I)")system%no_particles
WRITE(unit_vtk)character_buffer,char(10),'SCALARS Particle_radius double',char(10),'LOOKUP_TABLE default',char(10),(real(system%particle(i)%radius,real8),i=1,system%no_particles),'SCALARS Particle_temperature double',char(10),'LOOKUP_TABLE default',char(10),(real(system%particle(i)%temperature,real8),i=1,system%no_particles),'VECTORS Particle_velocity double',char(10),(real(system%particle(i)%velocity,real8),i=1,system%no_particles)
IF (system%general_cv%rotational_dofs=="on" .OR. system%general_cv%print_rotational_dofs=="yes") THEN
    WRITE(unit_vtk)'VECTORS Particle_spin double',char(10),(real(system%particle(i)%spin,real8),i=1,system%no_particles),'VECTORS Particle_angle double',char(10),(real(system%particle(i)%angle,real8),i=1,system%no_particles)
END IF
CLOSE(unit=unit_vtk)


END SUBROUTINE print_vtk_temp
!!---------------------------------------------------------------------------------------------------------------------

 
SUBROUTINE print_thermomechanics_extra_results(system)

! dummy argument
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i, k, s, ss, no_infiltrated_particles, no_particles_that_crossed_midsection
REAL(KIND=real8)   :: time, vol_infiltrated_particles, vol_particles_that_crossed_midsection
REAL(KIND=real4)   :: time_real4, first_tprint_real4
REAL(KIND=real8), PARAMETER :: pi=4.0_real8*ATAN(1.0_real8)


! get step, substep and time info
s = system%solution_cv%current_step
ss = system%solution_cv%step(s)%current_substep
time = system%solution_cv%step(s)%current_time

! transform time values to single precision (to avoid round-off errors)
time_real4 = REAL(time,real4)
first_tprint_real4 = REAL(system%general_cv%dt_for_results_printing,real4)

! print system name (only at the time of first print)
IF (time_real4==first_tprint_real4) THEN 
  WRITE(extra_results_unit,*)
  WRITE(extra_results_unit,*) "System name: ", system%name
  WRITE(extra_results_unit,*)
END IF

! print step, substep and time info
WRITE(extra_results_unit,*) "==============================================================================="
WRITE(extra_results_unit,"(A,T9,I2,T16,A,T26,I10,T46,A,T53,ES15.6)") " Step =", s, "Substep =", ss, "Time =", time 
WRITE(extra_results_unit,*)

! print extra results (add below any desired result upon demand)

! print walls´ contact-related quantities (number of contacting particles and total contact force)
WRITE(extra_results_unit,*) "Rigid walls contact-related quantities:" 
WRITE(extra_results_unit,"(3X,A,T10,A,T22,A,T37,A,T52,A,T67,A,T82,A,T97,A)") "Wall", "ncpart", "fconx", "fcony", "fconz", "ffricx", "ffricy", "ffricz"
DO k=1,system%no_rigid_walls
  WRITE(extra_results_unit,"(I6,T8,I6,T16,6ES15.6)") k, system%wall(k)%no_contacting_particles, system%wall(k)%contact_force, system%wall(k)%friction_force
END DO
WRITE(extra_results_unit,*)

! compute number and volume of infiltrated particles into porous solid (infiltration model, joint paper with Tarek, Feb/2018)
no_infiltrated_particles = 0
vol_infiltrated_particles = 0.0_real8
no_particles_that_crossed_midsection = 0
vol_particles_that_crossed_midsection = 0.0_real8
DO i=1,system%no_particles
  IF (system%particle(i)%position(1)>=0.0545_real8 .AND. system%particle(i)%parameters%kind=="granular_sphere_in_given_fluid_flow") THEN
    no_infiltrated_particles = no_infiltrated_particles + 1
    vol_infiltrated_particles = vol_infiltrated_particles + 1.333333333333333_real8*pi*system%particle(i)%radius**3
  END IF
  IF (system%particle(i)%position(1)>=0.1_real8 .AND. system%particle(i)%parameters%kind=="granular_sphere_in_given_fluid_flow") THEN
    no_particles_that_crossed_midsection = no_particles_that_crossed_midsection + 1
    vol_particles_that_crossed_midsection = vol_particles_that_crossed_midsection + 1.333333333333333_real8*pi*system%particle(i)%radius**3
  END IF
END DO

! print number and volume of particles that crossed porous solid´s front face (infiltration model, joint paper with Tarek, Feb/2018)
WRITE(extra_results_unit,*) "Infiltration results: number and volume of particles for which x>=0.0545m" 
WRITE(extra_results_unit,"(3X,A,T27,A)") "number of particles", "volume of particles"
WRITE(extra_results_unit,"(3X,I10,T27,ES15.6)") no_infiltrated_particles, vol_infiltrated_particles
WRITE(extra_results_unit,*)

! print number and volume of particles that crossed porous solid´s midsection (infiltration model, joint paper with Tarek, Feb/2018)
WRITE(extra_results_unit,*) "Infiltration results: number and volume of particles for which x>=0.1m" 
WRITE(extra_results_unit,"(3X,A,T27,A)") "number of particles", "volume of particles"
WRITE(extra_results_unit,"(3X,I10,T27,ES15.6)") no_particles_that_crossed_midsection, vol_particles_that_crossed_midsection
WRITE(extra_results_unit,*)

! print forced walls positions and velocities(testing the implementation of forced walls)
WRITE(extra_results_unit,*) "Forced walls positions and velocities" 
WRITE(extra_results_unit,"(3X,A,T16,A,T31,A,T46,A,T61,A,T76,A,T91,A)") "Wall", "x", "y", "z", "vx", "vy", "vz"
DO k=1,system%no_rigid_walls
  IF (system%wall(k)%kind=="forced_flat_rigid_wall") THEN
    WRITE(extra_results_unit,"(I6,T8,6ES15.6)") k, system%wall(k)%point_position, system%wall(k)%velocity
  END IF
END DO
WRITE(extra_results_unit,*)

! print system properties if this is more convenient here than in the .res file (energy, momenta, center of mass, etc)
!IF (system%general_cv%compute_system_properties=="yes") THEN
!  IF (system%general_cv%rotational_dofs=="off") THEN
!    system%translational_kinetic_energy = system%kinetic_energy
!    system%rotational_kinetic_energy = 0.0_real8
!  END IF
!  WRITE(extra_results_unit,"(6X,A,T45,3ES15.6)") "System kinetic energy (tot,trans,rot):", system%kinetic_energy, &
!                                                system%translational_kinetic_energy, system%rotational_kinetic_energy
!  WRITE(extra_results_unit,"(6X,A,T45,3ES15.6)") "System linear momentum:", system%linear_mom
!  WRITE(extra_results_unit,"(6X,A,T45,3ES15.6)") "System center of mass:", system%center_of_mass
!  WRITE(extra_results_unit,*)
!END IF

! add here any additional desired results (ex: particles´ forces, accelerations, etc)

END SUBROUTINE print_thermomechanics_extra_results
!!---------------------------------------------------------------------------------------------------------------------


END MODULE print_results_for_thermomechanics_solver

!!=============================================================================================================================
!!=============================================================================================================================


MODULE solvers_utilities

USE particles_classes
USE array_of_vectors_class
USE contact_history_array_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)

! dummy arguments
INTEGER(KIND=int4) :: niter, desired_niter
REAL(KIND=real8) :: dt_min, dt_max
REAL(KIND=real8), INTENT(INOUT) :: dt


! compute new dt
dt = (DSQRT(DBLE(desired_niter)/DBLE(niter)))*dt

! avoid being trapped
!!IF (dt==dt_min .AND. (niter==desired_niter .OR. niter==2)) dt=1.2_real8*dt
!IF ((niter==desired_niter .OR. niter==2) .OR. dt==dt_min) dt=1.2_real8*dt

! avoid being trapped
!IF (dt<dt_max .AND. niter<=desired_niter) THEN
! dt=1.2_real8*dt
!END IF

! enforce limits
IF (dt>dt_max) dt=dt_max  !dt = MIN(dt,system%solution_cv%step(s)%dt_max)
IF (dt<dt_min) dt=dt_min

END SUBROUTINE adapt_time_step_size
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_iteration_errors(Rresidual,Vresidual,Aresidual,Wresidual,Rnormalizer,Vnormalizer,Anormalizer,Wnormalizer, &
                                    Rtpdt,Vtpdt,Adelta,Wtpdt,Rerror,Verror,Aerror,Werror,tol1,tol2,tolR,tolV,tolA,tolW)

! dummy arguments
REAL(KIND=real8), INTENT(IN)  :: tol1, tol2
REAL(KIND=real8), INTENT(OUT) :: Rerror, Verror, Aerror, Werror, tolR, tolV, tolA, tolW
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rresidual, Vresidual, Aresidual, Wresidual, Rnormalizer, &
                                                 Vnormalizer, Anormalizer, Wnormalizer, Rtpdt, Vtpdt, Adelta, Wtpdt

! local variables
REAL(KIND=real8) :: norm_Rnormalizer, norm_Vnormalizer, norm_Anormalizer, norm_Wnormalizer, &
                    norm_Rresidual, norm_Vresidual, norm_Aresidual, norm_Wresidual, smallest_normalizer


! compute residual norms
norm_Rresidual = norm(Rresidual)
norm_Vresidual = norm(Vresidual)
norm_Aresidual = norm(Aresidual)
norm_Wresidual = norm(Wresidual)

! compute normalizer norms
norm_Rnormalizer = norm(Rnormalizer)
norm_Vnormalizer = norm(Vnormalizer)
norm_Anormalizer = norm(Anormalizer)
norm_Wnormalizer = norm(Wnormalizer)

! initialize tolerances
tolR = tol1
tolV = tol2
tolA = tol1
tolW = tol2

! set normalizers and tolerances for positions
smallest_normalizer = 1.0D-14*norm(Rtpdt)
IF (norm_Rnormalizer < smallest_normalizer) THEN
  norm_Rnormalizer = norm(Rtpdt)
  tolR = 0.01_real8*tolR
  IF (norm_Rnormalizer < 1.0D-14) THEN
    norm_Rnormalizer = 1.0_real8
    tolR = 1.0D-15
  END IF  
END IF

! set normalizers and tolerances for velocities
smallest_normalizer = 1.0D-14*norm(Vtpdt)
IF (norm_Vnormalizer < smallest_normalizer) THEN
  norm_Vnormalizer = norm(Vtpdt)
  tolV = 0.01_real8*tolV
  IF (norm_Vnormalizer < 1.0D-14) THEN
    norm_Vnormalizer = 1.0_real8
    tolV = 1.0D-15
  END IF  
END IF

! set normalizers and tolerances for angles
smallest_normalizer = 1.0D-14
IF (norm_Anormalizer < smallest_normalizer) THEN
  norm_Anormalizer = 1.0_real8
  tolA = 1.0D-15
END IF

! set normalizers and tolerances for spins
smallest_normalizer = 1.0D-14*norm(Wtpdt)
IF (norm_Wnormalizer < smallest_normalizer) THEN
  norm_Wnormalizer = norm(Wtpdt)  
  tolW = 0.01_real8*tolW
  IF (norm_Wnormalizer < 1.0D-14) THEN
    norm_Wnormalizer = 1.0_real8
    tolW = 1.0D-15
  END IF  
END IF

! double-check for small normalizers
IF (norm_Rnormalizer<=TINY(1.0_real8)) norm_Rnormalizer = 1.0_real8
IF (norm_Vnormalizer<=TINY(1.0_real8)) norm_Vnormalizer = 1.0_real8
IF (norm_Anormalizer<=TINY(1.0_real8)) norm_Anormalizer = 1.0_real8
IF (norm_Wnormalizer<=TINY(1.0_real8)) norm_Wnormalizer = 1.0_real8

! compute iteration errors
Rerror = norm_Rresidual/norm_Rnormalizer
Verror = norm_Vresidual/norm_Vnormalizer
Aerror = norm_Aresidual/norm_Anormalizer
Werror = norm_Wresidual/norm_Wnormalizer

END SUBROUTINE compute_iteration_errors
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE compute_iteration_errors_without_rotdofs(Rresidual,Vresidual,Rnormalizer,Vnormalizer,&
                                                    Rtpdt,Vtpdt,Rerror,Verror,tol1,tol2,tolR,tolV)

! dummy arguments
REAL(KIND=real8), INTENT(IN)  :: tol1, tol2
REAL(KIND=real8), INTENT(OUT) :: Rerror, Verror, tolR, tolV
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rresidual, Vresidual, Rnormalizer, Vnormalizer, Rtpdt, Vtpdt

! local variables
REAL(KIND=real8) :: norm_Rnormalizer, norm_Vnormalizer, norm_Rresidual, norm_Vresidual, smallest_normalizer


! compute residual norms
norm_Rresidual = norm(Rresidual)
norm_Vresidual = norm(Vresidual)

! compute normalizer norms
norm_Rnormalizer = norm(Rnormalizer)
norm_Vnormalizer = norm(Vnormalizer)

! initialize tolerances
tolR = tol1
tolV = tol2

! set normalizers and tolerances for positions
smallest_normalizer = 1.0D-14*norm(Rtpdt)
IF (norm_Rnormalizer < smallest_normalizer) THEN
  norm_Rnormalizer = norm(Rtpdt)
  tolR = 0.01_real8*tolR
  IF (norm_Rnormalizer < 1.0D-14) THEN
    norm_Rnormalizer = 1.0_real8
    tolR = 1.0D-15
  END IF  
END IF

! set normalizers and tolerances for velocities
smallest_normalizer = 1.0D-14*norm(Vtpdt)
IF (norm_Vnormalizer < smallest_normalizer) THEN
  norm_Vnormalizer = norm(Vtpdt)
  tolV = 0.01_real8*tolV
  IF (norm_Vnormalizer < 1.0D-14) THEN
    norm_Vnormalizer = 1.0_real8
    tolV = 1.0D-15
  END IF  
END IF

! double-check for small normalizers
IF (norm_Rnormalizer<=TINY(1.0_real8)) norm_Rnormalizer = 1.0_real8
IF (norm_Vnormalizer<=TINY(1.0_real8)) norm_Vnormalizer = 1.0_real8

! compute iteration errors
Rerror = norm_Rresidual/norm_Rnormalizer
Verror = norm_Vresidual/norm_Vnormalizer

END SUBROUTINE compute_iteration_errors_without_rotdofs
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE get_particle_contact_history_arrays(i,system,PPCHt,PWCHt)

! dummy arguments
INTEGER(KIND=int4) :: i
TYPE(system_data)  :: system
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHt, PWCHt

! local variables
INTEGER(KIND=int4) :: j, ncp, ncw


! deallocate old arrays
!IF (ASSOCIATED(PPCHt(i)%pair)) DEALLOCATE(PPCHt(i)%pair)
!IF (ASSOCIATED(PWCHt(i)%pair)) DEALLOCATE(PWCHt(i)%pair)

! deallocate old arrays
IF (ALLOCATED(PPCHt(i)%pair)) DEALLOCATE(PPCHt(i)%pair)
IF (ALLOCATED(PWCHt(i)%pair)) DEALLOCATE(PWCHt(i)%pair)

! allocate new arrays
ncp = system%particle(i)%no_contacting_particles
ncw = system%particle(i)%no_contacting_walls
ALLOCATE(PPCHt(i)%pair(1:ncp))
ALLOCATE(PWCHt(i)%pair(1:ncw))

! pass contact history of particle i to new arrays
DO j=1,ncp
  PPCHt(i)%pair(j)%number = system%particle(i)%contacting_particles(j)%part_number
  PPCHt(i)%pair(j)%normalspring_elongation = system%particle(i)%contacting_particles(j)%normalspring_elongation
  PPCHt(i)%pair(j)%normalspring_force = system%particle(i)%contacting_particles(j)%normalspring_force
  PPCHt(i)%pair(j)%fricspring_elongation = system%particle(i)%contacting_particles(j)%fricspring_elongation
  PPCHt(i)%pair(j)%fricspring_force = system%particle(i)%contacting_particles(j)%fricspring_force
  PPCHt(i)%pair(j)%rolresspring_elongation = system%particle(i)%contacting_particles(j)%rolresspring_elongation
END DO
DO j=1,ncw
  PWCHt(i)%pair(j)%number = system%particle(i)%contacting_walls(j)%wall_number
  PWCHt(i)%pair(j)%normalspring_elongation = system%particle(i)%contacting_walls(j)%normalspring_elongation
  PWCHt(i)%pair(j)%normalspring_force = system%particle(i)%contacting_walls(j)%normalspring_force
  PWCHt(i)%pair(j)%fricspring_elongation = system%particle(i)%contacting_walls(j)%fricspring_elongation
  PWCHt(i)%pair(j)%fricspring_force = system%particle(i)%contacting_walls(j)%fricspring_force
  PWCHt(i)%pair(j)%rolresspring_elongation = system%particle(i)%contacting_walls(j)%rolresspring_elongation
END DO

END SUBROUTINE get_particle_contact_history_arrays
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE update_particle_contact_history_arrays(i,system,PPCH,PWCH)

! dummy arguments
INTEGER(KIND=int4) :: i
TYPE(system_data)  :: system
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCH, PWCH

! local variables
INTEGER(KIND=int4) :: j, ncp, ncw


! deallocate old arrays
!IF (ASSOCIATED(PPCHtpdt(i)%pair)) DEALLOCATE(PPCHtpdt(i)%pair)
!IF (ASSOCIATED(PWCHtpdt(i)%pair)) DEALLOCATE(PWCHtpdt(i)%pair)

! deallocate old arrays
IF (ALLOCATED(PPCH(i)%pair)) DEALLOCATE(PPCH(i)%pair)
IF (ALLOCATED(PWCH(i)%pair)) DEALLOCATE(PWCH(i)%pair)

! allocate new arrays
ncp = system%particle(i)%no_contacting_particles
ncw = system%particle(i)%no_contacting_walls
ALLOCATE(PPCH(i)%pair(1:ncp))
ALLOCATE(PWCH(i)%pair(1:ncw))

! save contact history of particle i on new arrays
DO j=1,ncp
  PPCH(i)%pair(j)%number = system%particle(i)%contacting_particles(j)%part_number
  PPCH(i)%pair(j)%normalspring_elongation = system%particle(i)%contacting_particles(j)%normalspring_elongation
  PPCH(i)%pair(j)%normalspring_force = system%particle(i)%contacting_particles(j)%normalspring_force
  PPCH(i)%pair(j)%fricspring_elongation = system%particle(i)%contacting_particles(j)%fricspring_elongation
  PPCH(i)%pair(j)%fricspring_force = system%particle(i)%contacting_particles(j)%fricspring_force
  PPCH(i)%pair(j)%rolresspring_elongation = system%particle(i)%contacting_particles(j)%rolresspring_elongation
END DO
DO j=1,ncw
  PWCH(i)%pair(j)%number = system%particle(i)%contacting_walls(j)%wall_number
  PWCH(i)%pair(j)%normalspring_elongation = system%particle(i)%contacting_walls(j)%normalspring_elongation
  PWCH(i)%pair(j)%normalspring_force = system%particle(i)%contacting_walls(j)%normalspring_force
  PWCH(i)%pair(j)%fricspring_elongation = system%particle(i)%contacting_walls(j)%fricspring_elongation
  PWCH(i)%pair(j)%fricspring_force = system%particle(i)%contacting_walls(j)%fricspring_force
  PWCH(i)%pair(j)%rolresspring_elongation = system%particle(i)%contacting_walls(j)%rolresspring_elongation
END DO

END SUBROUTINE update_particle_contact_history_arrays
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE deallocate_contact_history_arrays(np,PPCHt,PWCHt,PPCHtpdt,PWCHtpdt)

! dummy arguments
INTEGER(KIND=int4) :: np
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHt, PWCHt
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE, OPTIONAL :: PPCHtpdt, PWCHtpdt

! local variables
INTEGER(KIND=int4) :: i, derror1, derror2, derror3, derror4, derror5, derror6


! deallocate arrays of time t
DO i=1,np
  DEALLOCATE(PPCHt(i)%pair,STAT=derror1)
  DEALLOCATE(PWCHt(i)%pair,STAT=derror2)
END DO
DEALLOCATE(PPCHt,PWCHt,STAT=derror3)

! deallocate arrays of time t+dt
derror4 = 0; derror5 = 0; derror6 = 0
IF (PRESENT(PPCHtpdt) .AND. PRESENT(PWCHtpdt)) THEN
  DO i=1,np
    DEALLOCATE(PPCHtpdt(i)%pair,STAT=derror4) !IF (ASSOCIATED(PPCHtpdt(i)%pair)) DEALLOCATE(PPCHtpdt(i)%pair)
    DEALLOCATE(PWCHtpdt(i)%pair,STAT=derror5) !IF (ASSOCIATED(PWCHtpdt(i)%pair)) DEALLOCATE(PWCHtpdt(i)%pair)
  END DO
DEALLOCATE(PPCHtpdt,PWCHtpdt,STAT=derror6)
END IF

! check for deallocation error
derror6 = derror1 + derror2 + derror3 + derror4 + derror5 + derror6
IF (derror6/=0) WRITE(error_unit,*) "deallocation failure of contact history arrays at the end of solver"

END SUBROUTINE deallocate_contact_history_arrays
!!---------------------------------------------------------------------------------------------------------------------


END MODULE solvers_utilities

!!=============================================================================================================================
!!=============================================================================================================================


MODULE mechanics_solver_class

USE particles_classes
USE rigid_walls_classes
USE periodic_boundary_conditions
USE grid_subroutines
USE verlet_lists_subroutines
USE contact_history_array_class
USE contacts_solver
USE dof_constraints
USE system_properties
USE print_results
USE solvers_utilities

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE mechanics_problem_solver(system)

! dummy arguments
TYPE(system_data) :: system


! select time integration solver
SELECT CASE (system%solution_cv%solver_type)

  CASE("midpoint_implicit_solver_for_hard_spheres_worot")
    CALL midpoint_implicit_solver_for_hard_spheres_worot(system)
    
  CASE("midpoint_implicit_solver_for_soft_spheres_worot")
    CALL midpoint_implicit_solver_for_soft_spheres_worot(system)

  CASE("euler_explicit_solver_for_soft_spheres_worot")
    CALL euler_explicit_solver_for_soft_spheres_worot(system)

  CASE("midpoint_implicit_solver_for_hard_spheres")
    CALL midpoint_implicit_solver_for_hard_spheres(system)

  CASE("euler_explicit_solver_for_soft_spheres")
    CALL euler_explicit_solver_for_soft_spheres(system)
  
  CASE("midpoint_implicit_solver_for_soft_spheres")
    CALL midpoint_implicit_solver_for_soft_spheres(system)

  CASE("newmark_implicit_solver")
    CALL newmark_implicit_solver(system)

  CASE DEFAULT
    WRITE(error_unit,*) "Solver type not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE mechanics_problem_solver
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE midpoint_implicit_solver_for_hard_spheres_worot(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind
INTEGER(KIND=int4) :: s, ss, i, j, w, iter, niter, max_niter, desired_niter, &
                      dealloc_error_R, dealloc_error_V, dealloc_error_F, dealloc_error_Rres_Vres_and_normalizers
REAL(KIND=real8)   :: t, tpdt, tprint, dt, dt_min, dt_max, dtcol, dt_print, phi, phi2, mi, tol1, tol2, tolR, tolV, Rerror, Verror 
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, finft, fienvt, finftpdt, fienvtpdt, &
                                  fiapplt, ri, rj, vi, ficont, fifrit, ficontpdt, fifritpdt
REAL(KIND=real8), PARAMETER    :: one = 1.0_real8, half=0.5_real8   
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, Fnft, Fnftpdt, Fenvt, Fenvtpdt, Fcont, Fcontpdt, &
                                                 Rtpdt_old, Vtpdt_old, Rresidual, Vresidual, Rnormalizer, Vnormalizer

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals refer to a single particle vector)
! Rt = vector of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = vector of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = vector of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = vector of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnftpdt = vector of near-field force vectors at time t+dt (contains the external force vectors of all particles at time t+dt)
! rit = position vector of particle i at time t
! ritpdt = position vector of particle i at time t+dt
! vit = velocity vector of particle i at time t
! vitpdt = velocity vector of particle i at time t+dt
! finft = near-field force vector of particle i at time t

! check for for proper solver type
IF (system%general_cv%rotational_dofs=="on") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match rotational dofs switch. Program aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match rotational dofs switch. Program aborted. "
  STOP
END IF

! allocate vectors
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles),Rtpdt_old(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles),Vtpdt_old(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fnftpdt(1:system%no_particles))
ALLOCATE(Fenvt(1:system%no_particles),Fenvtpdt(1:system%no_particles))
ALLOCATE(Fcont(1:system%no_particles),Fcontpdt(1:system%no_particles))
ALLOCATE(Rresidual(1:system%no_particles),Vresidual(1:system%no_particles))
ALLOCATE(Rnormalizer(1:system%no_particles),Vnormalizer(1:system%no_particles))

! get solution control variables
max_niter = system%solution_cv%max_no_iterations
desired_niter = system%solution_cv%desired_no_iterations
tol1 = system%solution_cv%positions_tolerance
tol2 = system%solution_cv%velocities_tolerance

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8

! get particles´ initial conditions
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
END DO

! compute particles´ mass and inertia
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
END DO

! clear contact-related quantities
DO i=1,system%no_particles
  system%particle(i)%contact_force = 0.0_real8
  system%particle(i)%friction_force = 0.0_real8
  system%particle(i)%friction_moment = 0.0_real8
  system%particle(i)%rolling_resistance_moment = 0.0_real8
END DO

! initialize particles´ force and moment vectors (forces and moments at t=0)
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL detect_and_resolve_particle_contacts(i,Rt,Vt,sys=system)  ! caution: if there are particles in contact at t=0, they must have nonzero relative velocity for the contact not to be missed!
  CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
  CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt)
  system%particle(i)%nearfield_force = finft
  system%particle(i)%environment_force = fienvt
END DO

! set phi-parameter (generalized mid-point rule)
phi = 0.5_real8
phi2 = 0.25_real8

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_max = system%solution_cv%step(s)%dt_max
  dt_min = system%solution_cv%step(s)%dt_min
  dt_print = system%general_cv%dt_for_results_printing

  ! initialize time values (tpdt is t+dt, i.e., is the time at the end of the substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)

  ! set step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt
  
  ! initialize substep counter
  ss = 0
  
  ! perform time stepping (one time step is named here a "substep")
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

     ! increment substep counter
     ss = ss + 1
     system%solution_cv%step(s)%current_substep = ss
          
     ! get particles positions, velocities and forces at the beginning of the substep
     DO i=1,system%no_particles
       Rt(i)%r = system%particle(i)%position
       Vt(i)%r = system%particle(i)%velocity
       Fnft(i)%r = system%particle(i)%nearfield_force
       Fenvt(i)%r = system%particle(i)%environment_force
     END DO

     ! print substep information on screen
  10 WRITE(default_output_unit,*)
     WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

     ! initialize collision-related quantities
     DO i=1,system%no_particles
       system%particle(i)%contact_force = 0.0_real8
       system%particle(i)%friction_force = 0.0_real8
       system%particle(i)%epsilon_force = Fnft(i)%r + Fenvt(i)%r
       system%particle(i)%collision_duration = 0.0_real8
     END DO
     
     ! set predictor for positions and velocities at t+dt
     DO i=1,system%no_particles
       Rtpdt(i)%r = system%particle(i)%position
       Vtpdt(i)%r = system%particle(i)%velocity
     END DO

     ! compute rigid walls positions, velocities and outside normals at time t+dt
     DO w=1,system%no_rigid_walls
       wkind = system%wall(w)%kind
       CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
       CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
     END DO

     ! perform iterations
     iterations: DO iter=1,max_niter
       
       ! save number of iterations performed
       niter = iter

       ! NOTE: in the IB-collisions solver, collision-related quantities cannot be cleared within the iterations, because
       ! it would clear the contact influence of one particle over the others in collisions involving more than 2 particles!

       ! loop over particles (compute iterative velocities and positions at t+dt)
       DO i=1,system%no_particles
         mi = system%particle(i)%mass
         i_kind = system%particle(i)%parameters%kind
         rit = Rt(i)%r
         vit = Vt(i)%r
         finft = Fnft(i)%r
         fienvt = Fenvt(i)%r
         CALL detect_and_resolve_particle_contacts(i,Rt,Vt,sys=system)
         ficont = system%particle(i)%contact_force
         fifrit = system%particle(i)%friction_force
         dtcol = system%particle(i)%collision_duration
         CALL compute_particle_nearfield_forces(i,i_kind,Rtpdt,Vtpdt,system,finftpdt)
         CALL compute_particle_environment_forces(i,i_kind,Rtpdt,Vtpdt,system,fienvtpdt) 
         vitpdt = vit + (dt/mi)*(phi*(finftpdt+fienvtpdt) + (one-phi)*(finft+fienvt)) + (dtcol/mi)*(ficont+fifrit)
         IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
         ritpdt = rit + (phi*vitpdt + (one-phi)*vit)*dt
         IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
         Rtpdt_old(i)%r = Rtpdt(i)%r
         Rtpdt(i)%r = ritpdt
         Vtpdt_old(i)%r = Vtpdt(i)%r
         Vtpdt(i)%r = vitpdt
         Fnftpdt(i)%r = finftpdt
         Fenvtpdt(i)%r = fienvtpdt
         Rresidual(i)%r = Rtpdt(i)%r - Rtpdt_old(i)%r
         Rnormalizer(i)%r = Rtpdt(i)%r - Rt(i)%r
         Vresidual(i)%r = Vtpdt(i)%r - Vtpdt_old(i)%r
         Vnormalizer(i)%r = Vtpdt(i)%r - Vt(i)%r
       END DO

       ! compute iteration errors and tolerances
       CALL compute_iteration_errors_without_rotdofs(Rresidual,Vresidual,Rnormalizer,Vnormalizer,&
                                                     Rtpdt,Vtpdt,Rerror,Verror,tol1,tol2,tolR,tolV)
       
  	   ! print iteration information
  	   WRITE(default_output_unit,*) "  iter =", INT(iter,int2), " error_R =", REAL(Rerror,real4), " error_V =", REAL(Verror,real4)

       ! check for convergence
       IF (Rerror<=tolR .AND. Verror<=tolV) EXIT iterations       

       ! check for divergence 
	   IF (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) EXIT iterations

     END DO iterations

     ! stop solution or adapt time ste in case of nonconvergence or divergence
	 IF (iter>max_niter .OR. (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) ) THEN
	   WRITE(error_unit,*) " nonconvergence or divergence at step =", s, ", t+dt =", REAL(tpdt,real4)
	   WRITE(error_unit,*) " time step size will be reduced if adaptive time stepping is ON"
	   WRITE(default_output_unit,*) " nonconvergence or divergence at step =", s, ", t+dt =", REAL(tpdt,real4)
	   WRITE(default_output_unit,*) " time step size will be reduced if adaptive time stepping is ON"
	   IF (system%solution_cv%step(s)%adaptive_time_stepping == "off") EXIT steps
	   IF (system%solution_cv%step(s)%adaptive_time_stepping == "on") THEN
         WRITE(error_unit,*) " time adaptivity still not implemented for this solver; program will be terminated"
         STOP
         !CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
	     !tpdt = t + dt
	     !system%solution_cv%step(s)%current_time = tpdt
         !system%solution_cv%step(s)%current_dt = dt
	     !GOTO 10
	   END IF  
	 END IF

     ! update particles positions, velocities and force vectors
     DO i=1,system%no_particles
       system%particle(i)%position = Rtpdt(i)%r 
       system%particle(i)%velocity = Vtpdt(i)%r
       system%particle(i)%nearfield_force = Fnftpdt(i)%r
       system%particle(i)%environment_force = Fenvtpdt(i)%r
     END DO 

     ! print substep results in results file
     IF (tpdt_real4==tprint_real4) THEN 
       IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
       CALL print_system_results(system)
       tprint = tprint + dt_print
       tprint_real4 = REAL(tprint,real4)
     END IF

     ! compute time step size for the next substep (adaptive time stepping)
     IF (system%solution_cv%step(s)%adaptive_time_stepping=="on" .AND. tpdt<system%solution_cv%step(s)%final_time) THEN
       WRITE(error_unit,*) " time adaptivity still not implemented for this solver; program will be terminated"
       STOP
       !CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
       !system%solution_cv%step(s)%current_dt = dt
     END IF

     ! increment time for next substep
     t = tpdt
     tpdt = tpdt + dt
     tpdt_real4 = REAL(tpdt,real4)

     ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)
     IF (tpdt_real4==tprint_real4) THEN
       tpdt = tprint
     END IF

     ! update time control variable
     system%solution_cv%step(s)%current_time = tpdt

  END DO substeps

END DO steps

! deallocate vectors
DEALLOCATE(Rt,Rtpdt,Rtpdt_old,STAT=dealloc_error_R)
DEALLOCATE(Vt,Vtpdt,Vtpdt_old,STAT=dealloc_error_V)
DEALLOCATE(Fcont,Fcontpdt,Fnft,Fnftpdt,Fenvt,Fenvtpdt,STAT=dealloc_error_F)
DEALLOCATE(Rresidual,Vresidual,Rnormalizer,Vnormalizer,STAT=dealloc_error_Rres_Vres_and_normalizers)

! check for deallocation error
dealloc_error_R = dealloc_error_R + dealloc_error_V + dealloc_error_F + dealloc_error_Rres_Vres_and_normalizers
IF (dealloc_error_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) at the end of the implicit solver"

END SUBROUTINE midpoint_implicit_solver_for_hard_spheres_worot
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE midpoint_implicit_solver_for_soft_spheres_worot(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind
INTEGER(KIND=int4) :: s, ss, i, j, w, iter, niter, max_niter, desired_niter, dealloc_error_R, dealloc_error_V, &
                      dealloc_error_F, dealloc_error_Rres_Vres_and_normalizers
REAL(KIND=real8)   :: t, tpdt, tprint, tcellupdate, tverletupdate, dt, dt_min, dt_max, dt_print, dt_cell_list_update, &
                      dt_verlet_list_update, phi, phi2, mi, tol1, tol2, tolR, tolV, Rerror, Verror
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4, tcellupdate_real4, tverletupdate_real4
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, finft, fienvt, finftpdt, fienvtpdt, fiapplt, & 
                                  ri, rj, vi, ficont, fifrit, ficontpdt, fifritpdt
REAL(KIND=real8), PARAMETER    :: one = 1.0_real8, half=0.5_real8   
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, Fnft, Fnftpdt, Fenvt, Fenvtpdt, Fcont, Fcontpdt, Ffrict, &
                                                 Ffrictpdt, Rtpdt_old, Vtpdt_old, Rresidual, Vresidual, Rnormalizer, Vnormalizer

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals refer to a single particle vector)
! Rt = vector of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = vector of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = vector of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = vector of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnftpdt = vector of near-field force vectors at time t+dt (contains the external force vectors of all particles at time t+dt)
! rit = position vector of particle i at time t
! ritpdt = position vector of particle i at time t+dt
! vit = velocity vector of particle i at time t
! vitpdt = velocity vector of particle i at time t+dt
! finft = near-field force vector of particle i at time t

! check for for proper solver type
IF (system%general_cv%rotational_dofs=="on") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match rotational dofs switch. Program aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match rotational dofs switch. Program aborted. "
  STOP
END IF

! allocate vectors
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles),Rtpdt_old(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles),Vtpdt_old(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fnftpdt(1:system%no_particles))
ALLOCATE(Fenvt(1:system%no_particles),Fenvtpdt(1:system%no_particles))
ALLOCATE(Fcont(1:system%no_particles),Fcontpdt(1:system%no_particles))
ALLOCATE(Ffrict(1:system%no_particles),Ffrictpdt(1:system%no_particles))
ALLOCATE(Rresidual(1:system%no_particles),Vresidual(1:system%no_particles))
ALLOCATE(Rnormalizer(1:system%no_particles),Vnormalizer(1:system%no_particles))

! get solution control variables
max_niter = system%solution_cv%max_no_iterations
desired_niter = system%solution_cv%max_no_iterations
tol1 = system%solution_cv%positions_tolerance
tol2 = system%solution_cv%velocities_tolerance

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8
tcellupdate = 0.0_real8
tverletupdate = 0.0_real8

! get particles´ initial conditions
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
END DO

! compute particles´ mass and inertia
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
END DO

! initialize grid and verlet lists
IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
  CALL sort_particles_into_grid_cells(system)
END IF
IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  CALL sort_particles_into_grid_cells(system)
  CALL build_particles_verlet_lists(system)
END IF

! clear contact-related quantities
DO i=1,system%no_particles
  system%particle(i)%contact_force = 0.0_real8
  system%particle(i)%friction_force = 0.0_real8
  system%particle(i)%friction_moment = 0.0_real8
  system%particle(i)%rolling_resistance_moment = 0.0_real8
END DO

! initialize particles´ force and moment vectors (forces and moments at t=0)
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL detect_and_resolve_particle_contacts(i,Rt,Vt,sys=system)
  CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
  CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt)
  system%particle(i)%nearfield_force = finft
  system%particle(i)%environment_force = fienvt
END DO

! set phi-parameter (generalized mid-point rule)
phi = 0.5_real8
phi2 = 0.25_real8

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_max = system%solution_cv%step(s)%dt_max
  dt_min = system%solution_cv%step(s)%dt_min
  dt_print = system%general_cv%dt_for_results_printing
  dt_cell_list_update = system%grid%dt_for_cell_list_update
  dt_verlet_list_update = system%general_cv%dt_for_verlet_list_update

  ! initialize time values (tpdt is t+dt, i.e., the time at the end of the substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  tcellupdate = tcellupdate + dt_cell_list_update
  tverletupdate = tverletupdate + dt_verlet_list_update
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tcellupdate_real4 = REAL(tverletupdate,real4)
  tverletupdate_real4 = REAL(tverletupdate,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)

  ! set step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt
  
  ! initialize substep counter
  ss = 0
  
  ! perform time stepping (one time step is named here a "substep")
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

     ! increment substep counter
     ss = ss + 1
     system%solution_cv%step(s)%current_substep = ss
          
     ! get particles positions, velocities and forces at time t
     DO i=1,system%no_particles
       Rt(i)%r = system%particle(i)%position
       Vt(i)%r = system%particle(i)%velocity
       Fcont(i)%r = system%particle(i)%contact_force
       Ffrict(i)%r = system%particle(i)%friction_force
       Fnft(i)%r  = system%particle(i)%nearfield_force
       Fenvt(i)%r = system%particle(i)%environment_force
     END DO

     ! print substep information on screen
  10 WRITE(default_output_unit,*)
     WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

     ! set predictor for positions and velocities at time t+dt
     DO i=1,system%no_particles
       Rtpdt(i)%r = system%particle(i)%position
       Vtpdt(i)%r = system%particle(i)%velocity
     END DO

     ! compute rigid walls positions, velocities and outside normals at time t+dt
     DO w=1,system%no_rigid_walls
       wkind = system%wall(w)%kind
       CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
       CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
     END DO

     ! perform iterations
     iterations: DO iter=1,max_niter
       
       ! save number of iterations performed
       niter = iter
       
       ! clear iterative contact and friction force vectors
       DO i=1,system%no_particles
         system%particle(i)%contact_force = 0.0_real8
         system%particle(i)%friction_force = 0.0_real8
       END DO  
       
       ! loop over particles (compute iterative velocities and positions at time t+dt)
       DO i=1,system%no_particles
         mi = system%particle(i)%mass
         i_kind = system%particle(i)%parameters%kind
         rit = Rt(i)%r
         vit = Vt(i)%r
         finft = Fnft(i)%r
         fienvt = Fenvt(i)%r
         ficont = Fcont(i)%r
         fifrit = Ffrict(i)%r
         CALL detect_and_resolve_particle_contacts(i,Rtpdt,Vtpdt,sys=system)
         ficontpdt = system%particle(i)%contact_force
         fifritpdt = system%particle(i)%friction_force
         CALL compute_particle_nearfield_forces(i,i_kind,Rtpdt,Vtpdt,system,finftpdt)
         CALL compute_particle_environment_forces(i,i_kind,Rtpdt,Vtpdt,system,fienvtpdt) 
         vitpdt = vit + (dt/mi)*(phi*(finftpdt+fienvtpdt) + (one-phi)*(finft+fienvt)) & 
                      + (dt/mi)*(phi*(ficontpdt+fifritpdt) + (one-phi)*(ficont+fifrit))
         IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
         ritpdt = rit + (phi*vitpdt + (one-phi)*vit)*dt
         IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
         Rtpdt_old(i)%r = Rtpdt(i)%r
         Rtpdt(i)%r = ritpdt
         Vtpdt_old(i)%r = Vtpdt(i)%r
         Vtpdt(i)%r = vitpdt
         Fcontpdt(i)%r = ficontpdt
         Ffrictpdt(i)%r = fifritpdt
         Fnftpdt(i)%r = finftpdt
         Fenvtpdt(i)%r = fienvtpdt
         Rresidual(i)%r = Rtpdt(i)%r - Rtpdt_old(i)%r
         Rnormalizer(i)%r = Rtpdt(i)%r - Rt(i)%r
         Vresidual(i)%r = Vtpdt(i)%r - Vtpdt_old(i)%r
         Vnormalizer(i)%r = Vtpdt(i)%r - Vt(i)%r
       END DO

       ! compute iteration errors and tolarances
       CALL compute_iteration_errors_without_rotdofs(Rresidual,Vresidual,Rnormalizer,Vnormalizer,&
                                                     Rtpdt,Vtpdt,Rerror,Verror,tol1,tol2,tolR,tolV)

  	   ! print iteration information
  	   WRITE(default_output_unit,*) "  iter =", INT(iter,int2), " error_R =", REAL(Rerror,real4), " error_V =", REAL(Verror,real4)
       
       ! check for convergence
       IF (Rerror<=tolR .AND. Verror<=tolV) EXIT iterations

       ! check for divergence 
	   IF (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) EXIT iterations

     END DO iterations 

     ! stop solution or adapt time ste in case of nonconvergence or divergence 
	 IF (iter>max_niter .OR. (Rerror>=1.0E+06 .OR. Verror>=1.0E+06)) THEN
	   WRITE(error_unit,*) " Failure to converge at substep =", ss, ", t+dt =", REAL(tpdt,real4)
	   WRITE(error_unit,*) " Program will be terminated"
	   WRITE(default_output_unit,*) " failure to converge at substep =", ss, ", t+dt =", REAL(tpdt,real4)
	   WRITE(default_output_unit,*) " Program will be terminated"
	   EXIT steps
	 END IF

     ! update particles positions, velocities and force vectors
     DO i=1,system%no_particles
       system%particle(i)%position = Rtpdt(i)%r 
       system%particle(i)%velocity = Vtpdt(i)%r
       system%particle(i)%nearfield_force = Fnftpdt(i)%r
       system%particle(i)%environment_force = Fenvtpdt(i)%r
       system%particle(i)%contact_force = Fcontpdt(i)%r
       system%particle(i)%friction_force = Ffrictpdt(i)%r
     END DO 

     ! print substep results in results file
     IF (tpdt_real4==tprint_real4) THEN 
       IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
       CALL print_system_results(system)
       tprint = tprint + dt_print
       tprint_real4 = REAL(tprint,real4)
     END IF

     ! update cell and verlet lists
     IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
       IF (tpdt_real4==tcellupdate_real4) THEN 
         CALL update_grid_cells_lists(system)
         tcellupdate = tcellupdate + dt_cell_list_update
         tcellupdate_real4 = REAL(tcellupdate,real4)
       END IF
     END IF
     IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
       IF (tpdt_real4==tverletupdate_real4) THEN
         CALL update_grid_cells_lists(system)
         CALL update_particles_verlet_lists(system)
         tverletupdate = tverletupdate + dt_verlet_list_update
         tverletupdate_real4 = REAL(tverletupdate,real4)
       END IF
     END IF

     ! compute dt for next substep (adaptive time stepping)
     IF (system%solution_cv%step(s)%adaptive_time_stepping=="on" .AND. tpdt<system%solution_cv%step(s)%final_time) THEN
       CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
     END IF

     ! increment time for next substep
     t = tpdt
     tpdt = tpdt + dt
     tpdt_real4 = REAL(tpdt,real4)

     ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)
     IF (tpdt_real4==tprint_real4) THEN
       tpdt = tprint
     END IF

     ! update time control variables
     system%solution_cv%step(s)%current_time = tpdt
     system%solution_cv%step(s)%current_dt = dt

  END DO substeps

END DO steps

! deallocate vectors
DEALLOCATE(Rt,Rtpdt,Rtpdt_old,STAT=dealloc_error_R)
DEALLOCATE(Vt,Vtpdt,Vtpdt_old,STAT=dealloc_error_V)
DEALLOCATE(Fcont,Fcontpdt,Ffrict,Ffrictpdt,Fnft,Fnftpdt,Fenvt,Fenvtpdt,STAT=dealloc_error_F)
DEALLOCATE(Rresidual,Vresidual,Rnormalizer,Vnormalizer,STAT=dealloc_error_Rres_Vres_and_normalizers)

! check for deallocation error
dealloc_error_R = dealloc_error_R + dealloc_error_V + dealloc_error_F + dealloc_error_Rres_Vres_and_normalizers
IF (dealloc_error_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) at the end of midpoint_implicit solver"

END SUBROUTINE midpoint_implicit_solver_for_soft_spheres_worot
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE euler_explicit_solver_for_soft_spheres_worot(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind
INTEGER(KIND=int4) :: s, ss, i, j, w, dealloc_error_R, dealloc_error_V, dealloc_error_F
REAL(KIND=real8)   :: t, tpdt, tprint, tcellupdate, tverletupdate, dt, dt_print, dt_cell_list_update, dt_verlet_list_update, & 
                      phi, phi2, mi, Force_norm, Ref_norm, xmin, xmax, ymin, ymax, zmin, zmax, xi, yi, zi 
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4, tcellupdate_real4, tverletupdate_real4
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, ri, rj, ficont, fifrit, finft, fienvt,  &
                                  finftpdt, fienvtpdt, ficontpdt, fifritpdt, figivt, vi
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, Fnft, Fenvt, Fcont, Ffrict
REAL(KIND=real8), PARAMETER    :: one = 1.0_real8   

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals refer to a single particle vector)
! Rt = vector of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = vector of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = vector of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = vector of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnft = vector of near-field force vectors at time t (contains the external force vectors of all particles at time t)
! rit = position vector of particle "i" at time t
! ritpdt = position vector of particle "i" at time t+dt
! vit = velocity vector of particle "i" at time t
! vitpdt = velocity vector of particle "i" at time t+dt
! finft = near-field force vector of particle "i" at time t


! check for for proper solver type
IF (system%general_cv%rotational_dofs=="on") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match rotational dofs switch. Program aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match rotational dofs switch. Program aborted. "
  STOP
END IF

! allocate vectors
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fenvt(1:system%no_particles))
ALLOCATE(Fcont(1:system%no_particles),Ffrict(1:system%no_particles))

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8
tcellupdate = 0.0_real8
tverletupdate = 0.0_real8

! get particles´ initial conditions
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
END DO

! compute particles´ mass and inertia
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
END DO

! initialize grid and verlet lists
IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
  CALL sort_particles_into_grid_cells(system)
END IF
IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  CALL sort_particles_into_grid_cells(system)
  CALL build_particles_verlet_lists(system)
END IF

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_print = system%general_cv%dt_for_results_printing
  dt_cell_list_update = system%grid%dt_for_cell_list_update
  dt_verlet_list_update = system%general_cv%dt_for_verlet_list_update  

  ! initialize time values (tpdt is t+dt, i.e., the time at the end of the substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  tcellupdate = tcellupdate + dt_cell_list_update
  tverletupdate = tverletupdate + dt_verlet_list_update
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tcellupdate_real4 = REAL(tverletupdate,real4)
  tverletupdate_real4 = REAL(tverletupdate,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)

  ! initialize step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt
  
  ! initialize substep counter
  ss = 0

  ! perform time stepping (one time step is named here a "substep")
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

    ! increment substep counter
    ss = ss + 1
    system%solution_cv%step(s)%current_substep = ss
     
    ! get particles positions and velocities at time t
    DO i=1,system%no_particles
      Rt(i)%r = system%particle(i)%position
      Vt(i)%r = system%particle(i)%velocity
    END DO 
     
     ! clear contact-related quantities
     DO i=1,system%no_particles
       system%particle(i)%contact_force = 0.0_real8
       system%particle(i)%friction_force = 0.0_real8
     END DO

     ! compute particles´ contact and friction force vectors at time t
     DO i=1,system%no_particles
       CALL detect_and_resolve_particle_contacts(i,Rt,Vt,sys=system)
       Fcont(i)%r = system%particle(i)%contact_force
       Ffrict(i)%r = system%particle(i)%friction_force
     END DO

     ! compute particles´ nearfield and environment force vectors at time t
     DO i=1,system%no_particles
       i_kind = system%particle(i)%parameters%kind
       CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
       CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt)
       Fnft(i)%r  = finft
       Fenvt(i)%r = fienvt
     END DO

     ! enforce periodic boundary conditions
     IF (system%general_cv%periodic_bc=="on") THEN
       CALL enforce_periodic_boundary_conditions(Rt,system)
     END IF

    ! print substep information on screen
    WRITE(default_output_unit,*)
    WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

    ! loop over particles (compute particles´ velocities and positions at time t+dt)
    DO i=1,system%no_particles
      mi = system%particle(i)%mass
      i_kind = system%particle(i)%parameters%kind
      rit = Rt(i)%r
      vit = Vt(i)%r
      finft  = Fnft(i)%r
      fienvt = Fenvt(i)%r
      ficont = Fcont(i)%r
      fifrit = Ffrict(i)%r
      figivt = 0.0_real8
      IF (s==1 .AND. ss==1) figivt = system%particle(i)%initial_given_force
      vitpdt = vit + (dt/mi)*(finft+fienvt+figivt) + (dt/mi)*(ficont+fifrit)
      IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
      ritpdt = rit + dt*vitpdt
      IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
      Rtpdt(i)%r = ritpdt
      Vtpdt(i)%r = vitpdt
    END DO
       
    ! check for numerical instability (not a rigorous check)
	!V_norm = norm(Vtpdt)
	!Ref_norm = norm(Vt)
	!IF (Ref_norm<=1.0E-10) Ref_norm = 1.0_real8
	!Relative_velocity = V_norm/Ref_norm
	!IF (Relative_velocity>=1.0E+08) THEN
	!   WRITE(error_unit,*) " CAUTION: large increase on velocities at substep =", ss, ", time =", REAL(time,real4)
	!   WRITE(error_unit,*) " Be aware of possible instabilities within the explicit solver"
	!   WRITE(default_output_unit,*) " CAUTION: large increase on velocities at substep =", ss, ", time =", REAL(time,real4)
	!   WRITE(default_output_unit,*) " Be aware of possible instabilities within the explicit solver"
	!END IF

    ! update particles´ positions and velocities
    DO i=1,system%no_particles
      system%particle(i)%position = Rtpdt(i)%r 
      system%particle(i)%velocity = Vtpdt(i)%r
    END DO

    ! print substep results in results file
    IF (tpdt_real4==tprint_real4) THEN 
      IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
      CALL print_system_results(system)
      tprint = tprint + dt_print
      tprint_real4 = REAL(tprint,real4)
    END IF

    ! update cell and verlet lists
    IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
      IF (tpdt_real4==tcellupdate_real4) THEN 
        CALL update_grid_cells_lists(system)
        tcellupdate = tcellupdate + dt_cell_list_update
        tcellupdate_real4 = REAL(tcellupdate,real4)
      END IF
    END IF
    IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
      IF (tpdt_real4==tverletupdate_real4) THEN
        CALL update_grid_cells_lists(system)
        CALL update_particles_verlet_lists(system)
        tverletupdate = tverletupdate + dt_verlet_list_update
        tverletupdate_real4 = REAL(tverletupdate,real4)
      END IF
    END IF

    ! update rigid walls positions, velocities and outside normals (values at t+dt, i.e., for next substep)
    DO w=1,system%no_rigid_walls
      wkind = system%wall(w)%kind
      CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
      CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
    END DO

    ! increment time for next substep
    t = tpdt
    tpdt = tpdt + dt
    tpdt_real4 = REAL(tpdt,real4)

    ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)
    IF (tpdt_real4==tprint_real4) THEN
      tpdt = tprint
    END IF

    ! update time control variable
    system%solution_cv%step(s)%current_time = tpdt

  END DO substeps

END DO steps

! deallocate vectors
DEALLOCATE(Rt,Rtpdt,STAT=dealloc_error_R)
DEALLOCATE(Vt,Vtpdt,STAT=dealloc_error_V)
DEALLOCATE(Fnft,Fenvt,Fcont,Ffrict,STAT=dealloc_error_F)

! check for deallocation error
dealloc_error_R = dealloc_error_R + dealloc_error_V + dealloc_error_F
IF (dealloc_error_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) at the end of the explicit solver"

END SUBROUTINE euler_explicit_solver_for_soft_spheres_worot
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE midpoint_implicit_solver_for_hard_spheres(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind
INTEGER(KIND=int4) :: s, ss, i, j, w, iter, niter, max_niter, desired_niter, &
                      dealloc_error_R, dealloc_error_V, dealloc_error_A, dealloc_error_W, dealloc_error_F, &
                      dealloc_error_Rres_Vres_and_normalizers, dealloc_error_Ares_Wres_and_normalizers
REAL(KIND=real8)   :: t, tpdt, tprint, dt, dt_min, dt_max, dtcol, dt_print, phi, phi2, mi, ji, &
                      tol1, tol2, tolR, tolV, tolA, tolW, Rerror, Verror, Aerror, Werror
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, ait, aidelta, wit, witpdt, finft, fienvt, finftpdt, fienvtpdt, &
                                  fiapplt, ficont, fifrit, ficontpdt, fifritpdt, mifrit, ri, rj, vi
REAL(KIND=real8), PARAMETER    :: one = 1.0_real8, half=0.5_real8, four=4.0_real8   
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, At, Adelta, Wt, Wtpdt, Fnft, Fnftpdt, Fenvt, Fenvtpdt, &
                                                 Rtpdt_old, Vtpdt_old, Adelta_old, Wtpdt_old, Rresidual, Vresidual, Aresidual, Wresidual, &
                                                 Rnormalizer, Vnormalizer, Anormalizer, Wnormalizer

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals refer to a single particle vector)
! Rt = vector of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = vector of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = vector of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = vector of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnftpdt = vector of near-field force vectors at time t+dt (contains the external force vectors of all particles at time t+dt)
! rit = position vector of particle i at time t
! ritpdt = position vector of particle i at time t+dt
! vit = velocity vector of particle i at time t
! vitpdt = velocity vector of particle i at time t+dt
! finft = near-field force vector of particle i at time t


! allocate vectors
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles),Rtpdt_old(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles),Vtpdt_old(1:system%no_particles))
ALLOCATE(At(1:system%no_particles),Adelta(1:system%no_particles),Adelta_old(1:system%no_particles))
ALLOCATE(Wt(1:system%no_particles),Wtpdt(1:system%no_particles),Wtpdt_old(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fnftpdt(1:system%no_particles))
ALLOCATE(Fenvt(1:system%no_particles),Fenvtpdt(1:system%no_particles))
ALLOCATE(Rresidual(1:system%no_particles),Vresidual(1:system%no_particles))
ALLOCATE(Aresidual(1:system%no_particles),Wresidual(1:system%no_particles))
ALLOCATE(Rnormalizer(1:system%no_particles),Vnormalizer(1:system%no_particles))
ALLOCATE(Anormalizer(1:system%no_particles),Wnormalizer(1:system%no_particles))

! get solution control variables
max_niter = system%solution_cv%max_no_iterations
desired_niter = system%solution_cv%desired_no_iterations
tol1 = system%solution_cv%positions_tolerance
tol2 = system%solution_cv%velocities_tolerance

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8

! get particles´ initial conditions
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  system%particle(i)%angle = system%particle(i)%initial_angles
  system%particle(i)%spin = system%particle(i)%initial_spin
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
  At(i)%r = system%particle(i)%angle
  Wt(i)%r = system%particle(i)%spin
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
END DO

! compute particles´ mass and inertia
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
END DO

! clear contact-related quantities
DO i=1,system%no_particles
  system%particle(i)%contact_force = 0.0_real8
  system%particle(i)%friction_force = 0.0_real8
  system%particle(i)%friction_moment = 0.0_real8
  system%particle(i)%rolling_resistance_moment = 0.0_real8
END DO

! initialize particles´ force and moment vectors (forces and moments at t=0)
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL detect_and_resolve_particle_contacts(i,Rt,Vt,Adelta,Wt,sys=system)  ! caution: if there are particles in contact at t=0, they must have nonzero relative velocity for the contact not to be missed!
  CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
  CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt)
  system%particle(i)%nearfield_force = finft
  system%particle(i)%environment_force = fienvt
END DO

! set phi-parameter (generalized mid-point rule)
phi = 0.5_real8
phi2 = 0.25_real8

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_max = system%solution_cv%step(s)%dt_max
  dt_min = system%solution_cv%step(s)%dt_min
  dt_print = system%general_cv%dt_for_results_printing

  ! initialize time values (tpdt is t+dt, i.e., is the time at the end of the substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)
    
  ! set step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt

  ! initialize substep counter
  ss = 0
  
  ! perform time stepping (one time step is named here a "substep")
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

     ! increment substep counter
     ss = ss + 1
     system%solution_cv%step(s)%current_substep = ss
          
     ! get particles positions, velocities, angles, spins and forces at time t
     DO i=1,system%no_particles
       Rt(i)%r = system%particle(i)%position
       Vt(i)%r = system%particle(i)%velocity
       At(i)%r = system%particle(i)%angle
       Wt(i)%r = system%particle(i)%spin
       Fnft(i)%r = system%particle(i)%nearfield_force
       Fenvt(i)%r = system%particle(i)%environment_force
       ! caution: maybe contact and friction should enter here!
     END DO

     ! print substep information
  10 WRITE(default_output_unit,*)
     WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

     ! initialize contact-related quantities
     DO i=1,system%no_particles
       system%particle(i)%contact_force = 0.0_real8
       system%particle(i)%friction_force = 0.0_real8
       system%particle(i)%epsilon_force = Fnft(i)%r + Fenvt(i)%r
       system%particle(i)%friction_moment = 0.0_real8
       system%particle(i)%rolling_resistance_moment = 0.0_real8
       system%particle(i)%collision_duration = 0.0_real8
     END DO
     
     ! set predictor for positions, velocities, angles and spins at time t+dt
     DO i=1,system%no_particles
       Rtpdt(i)%r = system%particle(i)%position
       Vtpdt(i)%r = system%particle(i)%velocity
       Adelta(i)%r = 0.0_real8
       Wtpdt(i)%r = system%particle(i)%spin
     END DO

     ! compute rigid walls positions, velocities and outside normals at time t+dt
     DO w=1,system%no_rigid_walls
       wkind = system%wall(w)%kind
       CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
       CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
     END DO

     ! perform iterations
     iterations: DO iter=1,max_niter
       
       ! save number of iterations performed
       niter = iter

       ! NOTE: in the IB-contacts solver, contact-related quantities cannot be cleared within the iterations, otherwise
       ! it would clear the contact influence of one particle over the others in contacts involving more than 2 particles!

       ! loop over particles (compute particles velocities, positions, spins and angles at t+dt)
       DO i=1,system%no_particles
         mi = system%particle(i)%mass
         ji = system%particle(i)%inertia
         i_kind = system%particle(i)%parameters%kind
         rit = Rt(i)%r
         vit = Vt(i)%r
         ait = At(i)%r
         wit = Wt(i)%r
         finft = Fnft(i)%r
         fienvt = Fenvt(i)%r
         CALL detect_and_resolve_particle_contacts(i,Rt,Vt,Adelta,Wtpdt,sys=system)
         ficont = system%particle(i)%contact_force
         fifrit = system%particle(i)%friction_force
         mifrit = system%particle(i)%friction_moment
         dtcol = system%particle(i)%collision_duration
         CALL compute_particle_nearfield_forces(i,i_kind,Rtpdt,Vtpdt,system,finftpdt)
         CALL compute_particle_environment_forces(i,i_kind,Rtpdt,Vtpdt,system,fienvtpdt) 
         vitpdt = vit + (dt/mi)*(phi*(finftpdt+fienvtpdt) + (one-phi)*(finft+fienvt)) + (dtcol/mi)*(ficont+fifrit)
         IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
         ritpdt = rit + (phi*vitpdt + (one-phi)*vit)*dt
         IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
         witpdt = wit + (dtcol/ji)*mifrit
         IF (ANY(system%particle(i)%rotational_dof_codes(4:6)/=0)) CALL enforce_spin_constraints(i,system,witpdt)
         aidelta = (phi*witpdt + (one-phi)*wit)*dt
         IF (ANY(system%particle(i)%rotational_dof_codes(1:3)/=0)) CALL enforce_angle_constraints(i,system,aidelta,witpdt)
         Rtpdt_old(i)%r = Rtpdt(i)%r
         Rtpdt(i)%r = ritpdt
         Vtpdt_old(i)%r = Vtpdt(i)%r
         Vtpdt(i)%r = vitpdt
         Adelta_old(i)%r = Adelta(i)%r
         Adelta(i)%r = aidelta
         Wtpdt_old(i)%r = Wtpdt(i)%r
         Wtpdt(i)%r = witpdt
         Fnftpdt(i)%r = finftpdt
         Fenvtpdt(i)%r = fienvtpdt
         Rresidual(i)%r = Rtpdt(i)%r - Rtpdt_old(i)%r
         Rnormalizer(i)%r = Rtpdt(i)%r - Rt(i)%r
         Vresidual(i)%r = Vtpdt(i)%r - Vtpdt_old(i)%r
         Vnormalizer(i)%r = Vtpdt(i)%r - Vt(i)%r
         Aresidual(i)%r = Adelta(i)%r - Adelta_old(i)%r
         Anormalizer(i)%r = Adelta(i)%r
         Wresidual(i)%r = Wtpdt(i)%r - Wtpdt_old(i)%r
         Wnormalizer(i)%r = Wtpdt(i)%r - Wt(i)%r
       END DO

       ! compuute iteration errors and tolerances
       CALL compute_iteration_errors(Rresidual,Vresidual,Aresidual,Wresidual,Rnormalizer,Vnormalizer,Anormalizer,Wnormalizer, &
                                     Rtpdt,Vtpdt,Adelta,Wtpdt,Rerror,Verror,Aerror,Werror,tol1,tol2,tolR,tolV,tolA,tolW)
       
  	   ! print iteration information
  	   WRITE(default_output_unit,20) " iter=", iter, "errR=", Rerror, "errV=", Verror, "errA=", Aerror, "errW=", Werror
       20 FORMAT(TR2,A,I3,TR2,A,ES10.3,TR2,A,ES10.3,TR2,A,ES10.3,TR2,A,ES10.3)

       ! check for convergence
       IF ( (Rerror<=tolR .AND. Verror<=tolV) .AND. (Aerror<=tolA .AND. Werror<=tolW) ) EXIT iterations

       ! check for divergence 
	   IF ( (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) .OR. (Aerror>=1.0E+06 .OR. Werror>=1.0E+06) ) EXIT iterations

     END DO iterations 

     ! stop solution or adapt time ste in case of nonconvergence or divergence
	 IF ( iter>max_niter .OR. (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) .OR. (Aerror>=1.0E+06 .OR. Werror>=1.0E+06) ) THEN
	   WRITE(error_unit,*) " nonconvergence or divergence at step =", s, ", t+dt =", REAL(tpdt,real4)
	   WRITE(error_unit,*) " time step size will be reduced if adaptive time stepping is ON"
	   WRITE(default_output_unit,*) " nonconvergence or divergence at step =", s, ", t+dt =", REAL(tpdt,real4)
	   WRITE(default_output_unit,*) " time step size will be reduced if adaptive time stepping is ON"
	   IF (system%solution_cv%step(s)%adaptive_time_stepping == "off") EXIT steps
	   IF (system%solution_cv%step(s)%adaptive_time_stepping == "on") THEN
         WRITE(error_unit,*) " time adaptivity still not implemented for this solver; program will be terminated"
         STOP
         !CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
	     !tpdt = t + dt
	     !system%solution_cv%step(s)%current_time = tpdt
         !system%solution_cv%step(s)%current_dt = dt
	     !GOTO 10
	   END IF  
	 END IF

     ! update particles positions, velocities, angles, spins and force/moment vectors
     DO i=1,system%no_particles
       system%particle(i)%position = Rtpdt(i)%r 
       system%particle(i)%velocity = Vtpdt(i)%r
       system%particle(i)%angle = (four/(four-DOT_PRODUCT(Adelta(i)%r,At(i)%r)))*(At(i)%r+Adelta(i)%r+half*(Adelta(i)%r.vector.At(i)%r))
       system%particle(i)%spin = Wtpdt(i)%r
       system%particle(i)%nearfield_force = Fnftpdt(i)%r
       system%particle(i)%environment_force = Fenvtpdt(i)%r
     END DO 

     ! print substep results in results file
     IF (tpdt_real4==tprint_real4) THEN 
       IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
       CALL print_system_results(system)
       tprint = tprint + dt_print
       tprint_real4 = REAL(tprint,real4)
     END IF

     ! compute time step size for the next substep (adaptive time stepping)
     IF (system%solution_cv%step(s)%adaptive_time_stepping=="on" .AND. tpdt<system%solution_cv%step(s)%final_time ) THEN
       WRITE(error_unit,*) " time adaptivity still not implemented for this solver; program will be terminated"
       STOP
       !CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
       !system%solution_cv%step(s)%current_dt = dt
     END IF

     ! increment time for next substep
     t = tpdt
     tpdt = tpdt + dt
     tpdt_real4 = REAL(tpdt,real4)

     ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)
     IF (tpdt_real4==tprint_real4) THEN
       tpdt = tprint
     END IF

     ! update time control variable
     system%solution_cv%step(s)%current_time = tpdt

  END DO substeps

END DO steps

! deallocate vectors
DEALLOCATE(Rt,Rtpdt,Rtpdt_old,STAT=dealloc_error_R)
DEALLOCATE(Vt,Vtpdt,Vtpdt_old,STAT=dealloc_error_V)
DEALLOCATE(At,Adelta,Adelta_old,STAT=dealloc_error_A)
DEALLOCATE(Wt,Wtpdt,Wtpdt_old,STAT=dealloc_error_W)
DEALLOCATE(Fnft,Fnftpdt,Fenvt,Fenvtpdt,STAT=dealloc_error_F)
DEALLOCATE(Rresidual,Vresidual,Rnormalizer,Vnormalizer,STAT=dealloc_error_Rres_Vres_and_normalizers)
DEALLOCATE(Aresidual,Wresidual,Anormalizer,Wnormalizer,STAT=dealloc_error_Ares_Wres_and_normalizers)

! check for deallocation error
dealloc_error_R = dealloc_error_R + dealloc_error_V + dealloc_error_A + dealloc_error_W + dealloc_error_F + &
                  dealloc_error_Rres_Vres_and_normalizers + dealloc_error_Ares_Wres_and_normalizers
IF (dealloc_error_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) at the end of the implicit solver"

END SUBROUTINE midpoint_implicit_solver_for_hard_spheres
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE euler_explicit_solver_for_soft_spheres(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind
INTEGER(KIND=int4) :: s, ss, i, j, w, derror_R, derror_V, derror_A, derror_W, derror_F
REAL(KIND=real8)   :: t, tpdt, tprint, tcellupdate, tverletupdate, dt, dt_print, dt_cell_list_update, dt_verlet_list_update, mi, ji, Ref_norm
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4, tcellupdate_real4, tverletupdate_real4, mir
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, ait, aidelta, wit, witpdt, ficont, fifrit, finft, &
                                  fienvt, mifrit, mirolt, figivt, migivt, fiadht, miadht
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, At, Adelta, Wt, Wtpdt, Fnft, Fenvt, Fcont, &
                                                 Ffrict, Fadht, Mfrict, Mrolt, Madht
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHt, PWCHt
REAL(KIND=real8), PARAMETER :: one = 1.0_real8, four = 4.0_real8, half = 0.5_real8

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals to a single particle vector)
! Rt = vector of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = vector of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = vector of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = vector of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnft = vector of near-field force vectors at time t (contains the external force vectors of all particles at time t)
! rit = position vector of particle "i" at time t
! ritpdt = position vector of particle "i" at time t+dt
! vit = velocity vector of particle "i" at time t
! vitpdt = velocity vector of particle "i" at time t+dt
! finft = near-field force vector of particle "i" at time t
! PPCHt = array of particle-particle contact histories at time t
! PWCHt = array of particle-wall contact histories at time t
! etc

! check for for proper solver type
IF (system%general_cv%rotational_dofs=="off") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match rotational dofs switch. Program execution aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match rotational dofs switch. Program execution aborted. "
  STOP
END IF

! allocate arrays
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles))
ALLOCATE(At(1:system%no_particles),Adelta(1:system%no_particles))
ALLOCATE(Wt(1:system%no_particles),Wtpdt(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fenvt(1:system%no_particles))
ALLOCATE(Fcont(1:system%no_particles),Ffrict(1:system%no_particles))
ALLOCATE(Mfrict(1:system%no_particles),Mrolt(1:system%no_particles))
ALLOCATE(Fadht(1:system%no_particles),Madht(1:system%no_particles))
ALLOCATE(PPCHt(1:system%no_particles),PWCHt(1:system%no_particles))

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8
tcellupdate = 0.0_real8
tverletupdate = 0.0_real8

! get particles´ initial conditions
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  system%particle(i)%angle = system%particle(i)%initial_angles
  system%particle(i)%spin = system%particle(i)%initial_spin
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
  At(i)%r = system%particle(i)%angle
  Wt(i)%r = system%particle(i)%spin
  Adelta(i)%r = 0.0_real8
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
END DO

! compute particles´ mass and inertia
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
END DO

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! THINK: maybe the grid, verlet lists and contact lists below shall be initialized within the steps!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! initialize grid and verlet lists
IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
  CALL sort_particles_into_grid_cells(system)
END IF
IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  CALL sort_particles_into_grid_cells(system)
  CALL build_particles_verlet_lists(system)
END IF

! initialize particles´ contact lists
DO i=1,system%no_particles
  CALL initialize_particle_contact_lists(i,system)
END DO  

! initialize rigid walls´ contact lists
DO w=1,system%no_rigid_walls
  CALL initialize_wall_contact_lists(w,system)
END DO  

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_print = system%general_cv%dt_for_results_printing
  dt_cell_list_update = system%grid%dt_for_cell_list_update
  dt_verlet_list_update = system%general_cv%dt_for_verlet_list_update  

  ! initialize time values (tpdt is t+dt, i.e. is the time at the end of the substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  tcellupdate = tcellupdate + dt_cell_list_update
  tverletupdate = tverletupdate + dt_verlet_list_update
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tcellupdate_real4 = REAL(tverletupdate,real4)
  tverletupdate_real4 = REAL(tverletupdate,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)

  ! initialize step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt

  ! initialize substep counter
  ss = 0
    
  ! perform time stepping (one time step is named here a "substep")
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

     ! increment substep counter
     ss = ss + 1
     system%solution_cv%step(s)%current_substep = ss
     
     ! get particles´ arrays at time t
     DO i=1,system%no_particles
       Rt(i)%r = system%particle(i)%position
       Vt(i)%r = system%particle(i)%velocity
       At(i)%r = system%particle(i)%angle
       Wt(i)%r = system%particle(i)%spin
     END DO 
     
     ! clear particle contact-related quantities
     DO i=1,system%no_particles
       system%particle(i)%contact_force = 0.0_real8
       system%particle(i)%friction_force = 0.0_real8
       system%particle(i)%adhesion_force = 0.0_real8
       system%particle(i)%friction_moment = 0.0_real8
       system%particle(i)%rolling_resistance_moment = 0.0_real8
       system%particle(i)%adhesion_moment = 0.0_real8
     END DO

     ! clear wall contact-related forces
     DO w=1,system%no_rigid_walls
       system%wall(w)%contact_force = 0.0_real8
       system%wall(w)%friction_force = 0.0_real8
     END DO
     
     ! compute particles´ contact, friction and adhesion force vectors and arrays at time t
     DO i=1,system%no_particles
       CALL detect_and_resolve_particle_contacts(i,Rt,Vt,Adelta,Wt,sys=system)
       ! notice: at this point, the walls´ contact lists and total contact forces have been fully computed!
       Fcont(i)%r = system%particle(i)%contact_force
       Ffrict(i)%r = system%particle(i)%friction_force
       Fadht(i)%r = system%particle(i)%adhesion_force
       Mfrict(i)%r = system%particle(i)%friction_moment
       Mrolt(i)%r = system%particle(i)%rolling_resistance_moment
       Madht(i)%r = system%particle(i)%adhesion_moment
       CALL update_particle_contact_history_arrays(i,system,PPCHt,PWCHt)
     END DO

     ! compute particles´ nearfield and environment force vectors and arrays at time t
     DO i=1,system%no_particles
       i_kind = system%particle(i)%parameters%kind
       CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
       CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt)
       Fnft(i)%r  = finft
       Fenvt(i)%r = fienvt
     END DO

     ! enforce periodic boundary conditions
     IF (system%general_cv%periodic_bc=="on") THEN
       CALL enforce_periodic_boundary_conditions(Rt,system)
     END IF
     
     ! print substep information on screen
     WRITE(default_output_unit,*)
     WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

    ! loop over particles
    DO i=1,system%no_particles
      mi = system%particle(i)%mass
      ji = system%particle(i)%inertia
      rit = Rt(i)%r
      vit = Vt(i)%r
      ait = At(i)%r
      wit = Wt(i)%r
      finft  = Fnft(i)%r
      fienvt = Fenvt(i)%r
      ficont = Fcont(i)%r
      fifrit = Ffrict(i)%r
      fiadht = Fadht(i)%r
      mifrit = Mfrict(i)%r
      mirolt = Mrolt(i)%r
      miadht = Madht(i)%r
      figivt = 0.0_real8
      migivt = 0.0_real8
      !IF (s==1 .AND. ss==1) figivt = system%particle(i)%initial_given_force
      !IF (s==1 .AND. ss==1) migivt = system%particle(i)%initial_given_moment
      vitpdt = vit + (dt/mi)*(finft+fienvt+figivt) + (dt/mi)*(ficont+fifrit+fiadht)      
      IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
      ritpdt = rit + dt*vitpdt
      IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
      witpdt = wit + (dt/ji)*(mifrit+mirolt+miadht+migivt)
      IF (ANY(system%particle(i)%rotational_dof_codes(4:6)/=0)) CALL enforce_spin_constraints(i,system,witpdt)
      aidelta = dt*witpdt
      IF (ANY(system%particle(i)%rotational_dof_codes(1:3)/=0)) CALL enforce_angle_constraints(i,system,aidelta,witpdt)
      Rtpdt(i)%r = ritpdt
      Vtpdt(i)%r = vitpdt
      Adelta(i)%r = aidelta
      Wtpdt(i)%r = witpdt
    END DO
       
    ! check for numerical instability (not a rigorous check)
	!V_norm = norm(Vtpdt)
	!Ref_norm = norm(Vt)
	!IF (Ref_norm<=1.0E-10) Ref_norm = 1.0_real8
	!Relative_velocity = V_norm/Ref_norm
	!IF (Relative_velocity>=1.0E+08) THEN
	!   WRITE(error_unit,*) " CAUTION: large increase on velocities at substep =", ss, ", t+dt =", REAL(tpdt,real4)
	!   WRITE(error_unit,*) " Be aware of possible instabilities within the explicit solver"
	!   WRITE(default_output_unit,*) " CAUTION: large increase on velocities at substep =", ss, ", t+dt =", REAL(tpdt,real4)
	!   WRITE(default_output_unit,*) " Be aware of possible instabilities within the explicit solver"
	!END IF

    ! update particles positions, velocities, angles and spins
    DO i=1,system%no_particles
      system%particle(i)%position = Rtpdt(i)%r 
      system%particle(i)%velocity = Vtpdt(i)%r
      system%particle(i)%angle = (four/(four-DOT_PRODUCT(Adelta(i)%r,At(i)%r)))*(At(i)%r+Adelta(i)%r+half*(Adelta(i)%r.vector.At(i)%r))
      system%particle(i)%spin = Wtpdt(i)%r
    END DO

    ! print substep results in results file
    IF (tpdt_real4==tprint_real4) THEN 
      IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
      CALL print_system_results(system)
      tprint = tprint + dt_print
      tprint_real4 = REAL(tprint,real4)
    END IF

    ! update cell and verlet lists
    IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
      IF (tpdt_real4==tcellupdate_real4) THEN 
        CALL update_grid_cells_lists(system)
        tcellupdate = tcellupdate + dt_cell_list_update
        tcellupdate_real4 = REAL(tcellupdate,real4)
      END IF
    END IF
    IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
      IF (tpdt_real4==tverletupdate_real4) THEN
        CALL update_grid_cells_lists(system)
        CALL update_particles_verlet_lists(system)
        tverletupdate = tverletupdate + dt_verlet_list_update
        tverletupdate_real4 = REAL(tverletupdate,real4)
      END IF
    END IF

    ! update rigid walls positions, velocities and outside normals (values at t+dt, i.e., for next substep)    
    DO w=1,system%no_rigid_walls
      wkind = system%wall(w)%kind
      CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
      CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
    END DO

    ! increment time for next substep
    t = tpdt
    tpdt = tpdt + dt
    tpdt_real4 = REAL(tpdt,real4)

    ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)
    IF (tpdt_real4==tprint_real4) THEN
      tpdt = tprint
    END IF

    ! update time control variable
    system%solution_cv%step(s)%current_time = tpdt

  END DO substeps

END DO steps

! deallocate arrays
DEALLOCATE(Rt,Rtpdt,STAT=derror_R)
DEALLOCATE(Vt,Vtpdt,STAT=derror_V)
DEALLOCATE(At,Adelta,STAT=derror_A)
DEALLOCATE(Wt,Wtpdt,STAT=derror_W)
DEALLOCATE(Fnft,Fenvt,Fcont,Ffrict,Fadht,Mfrict,Mrolt,Madht,STAT=derror_F)
CALL deallocate_contact_history_arrays(system%no_particles,PPCHt,PWCHt)

! check for deallocation error
derror_R = derror_R + derror_V + derror_A + derror_W + derror_F
IF (derror_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) at the end of the explicit solver"

END SUBROUTINE euler_explicit_solver_for_soft_spheres
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE midpoint_implicit_solver_for_soft_spheres(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind
INTEGER(KIND=int4) :: s, ss, i, j, w, iter, niter, max_niter, desired_niter, derror_R, derror_V, derror_A, derror_W, derror_F, derror_M, &
                      derror_Rres_Vres_and_normalizers, derror_Ares_Wres_and_normalizers, no_contacting_particles, no_contacting_walls
REAL(KIND=real8)   :: t, tpdt, tprint, tcellupdate, tverletupdate, dt, dt_min, dt_max, dt_print, dt_cell_list_update, &
                      dt_verlet_list_update, phi, phi2, mi, ji, tol1, tol2, tolR, tolV, tolA, tolW, Rerror, Verror, Aerror, Werror
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4, tcellupdate_real4, tverletupdate_real4
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, ait, aidelta, wit, witpdt, finft, finftpdt, fienvt, fienvtpdt, &
                                  ficont, ficontpdt, fifrit, fifritpdt, fiadht, fiadhtpdt, mifrit, mifritpdt, mirolt, &
                                  miroltpdt, miadht, miadhtpdt, figivt, migivt
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, At, Adelta, Wt, Wtpdt, Fnft, Fnftpdt, Fenvt, Fenvtpdt, Fcont, Fcontpdt, &
                                                 Ffrict, Ffrictpdt, Mfrict, Mfrictpdt, Mrolt, Mroltpdt, Fadht, Fadhtpdt, Madht, Madhtpdt, &
                                                 Rtpdt_old, Vtpdt_old, Adelta_old, Wtpdt_old, Rresidual, Vresidual, Aresidual, Wresidual, &
                                                 Rnormalizer, Vnormalizer, Anormalizer, Wnormalizer
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHt, PWCHt, PPCHtpdt, PWCHtpdt
REAL(KIND=real8), PARAMETER :: one = 1.0_real8, half=0.5_real8, four=4.0_real8   

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals refer to a single particle vector)
! Rt = array of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = array of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = array of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = array of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnftpdt = array of near-field force vectors at time t+dt (contains the external force vectors of all particles at time t+dt)
! PPCHt = array of particle-particle contact histories at time t
! PWCHt = array of particle-wall contact histories at time t
! rit = position vector of particle i at time t
! ritpdt = position vector of particle i at time t+dt
! vit = velocity vector of particle i at time t
! vitpdt = velocity vector of particle i at time t+dt
! finft = near-field force vector of particle i at time t
! etc

! check for for proper solver type
IF (system%general_cv%rotational_dofs=="off") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match rotational dofs switch. Program execution aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match rotational dofs switch. Program execution aborted. "
  STOP
END IF

! allocate arrays
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles),Rtpdt_old(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles),Vtpdt_old(1:system%no_particles))
ALLOCATE(At(1:system%no_particles),Adelta(1:system%no_particles),Adelta_old(1:system%no_particles))
ALLOCATE(Wt(1:system%no_particles),Wtpdt(1:system%no_particles),Wtpdt_old(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fnftpdt(1:system%no_particles))
ALLOCATE(Fenvt(1:system%no_particles),Fenvtpdt(1:system%no_particles))
ALLOCATE(Fcont(1:system%no_particles),Fcontpdt(1:system%no_particles))
ALLOCATE(Ffrict(1:system%no_particles),Ffrictpdt(1:system%no_particles))
ALLOCATE(Mfrict(1:system%no_particles),Mfrictpdt(1:system%no_particles))
ALLOCATE(Mrolt(1:system%no_particles),Mroltpdt(1:system%no_particles))
ALLOCATE(Fadht(1:system%no_particles),Fadhtpdt(1:system%no_particles))
ALLOCATE(Madht(1:system%no_particles),Madhtpdt(1:system%no_particles))
ALLOCATE(Rresidual(1:system%no_particles),Vresidual(1:system%no_particles))
ALLOCATE(Aresidual(1:system%no_particles),Wresidual(1:system%no_particles))
ALLOCATE(Rnormalizer(1:system%no_particles),Vnormalizer(1:system%no_particles))
ALLOCATE(Anormalizer(1:system%no_particles),Wnormalizer(1:system%no_particles))
ALLOCATE(PPCHt(1:system%no_particles),PPCHtpdt(1:system%no_particles))
ALLOCATE(PWCHt(1:system%no_particles),PWCHtpdt(1:system%no_particles))

! get solution control variables
max_niter = system%solution_cv%max_no_iterations
desired_niter = system%solution_cv%desired_no_iterations
tol1 = system%solution_cv%positions_tolerance
tol2 = system%solution_cv%velocities_tolerance

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8
tcellupdate = 0.0_real8
tverletupdate = 0.0_real8

! get particles´ initial conditions
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  system%particle(i)%angle = system%particle(i)%initial_angles
  system%particle(i)%spin = system%particle(i)%initial_spin
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
  At(i)%r = system%particle(i)%angle
  Wt(i)%r = system%particle(i)%spin
  Adelta(i)%r = 0.0_real8 
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
END DO

! compute particles´ mass and inertia
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
END DO

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! THINK: maybe the grid, verlet lists and contact lists below shall be initialized within the steps!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! initialize grid and verlet lists
IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
  CALL sort_particles_into_grid_cells(system)
END IF
IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  CALL sort_particles_into_grid_cells(system)
  CALL build_particles_verlet_lists(system)
END IF

! initialize particles´ contact lists
DO i=1,system%no_particles
  CALL initialize_particle_contact_lists(i,system)
END DO  

! set phi parameter (generalized mid-point rule)
phi = 0.5_real8
phi2 = 0.25_real8

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_max = system%solution_cv%step(s)%dt_max
  dt_min = system%solution_cv%step(s)%dt_min
  dt_print = system%general_cv%dt_for_results_printing
  dt_cell_list_update = system%grid%dt_for_cell_list_update
  dt_verlet_list_update = system%general_cv%dt_for_verlet_list_update

  ! initialize time values (tpdt is t+dt, i.e., is the time at the end of a substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  tcellupdate = tcellupdate + dt_cell_list_update
  tverletupdate = tverletupdate + dt_verlet_list_update
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tcellupdate_real4 = REAL(tverletupdate,real4)
  tverletupdate_real4 = REAL(tverletupdate,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)

  ! initialize step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt

  ! initialize substep counter
  ss = 0
  
  ! clear contact-related quantities
  DO i=1,system%no_particles
    system%particle(i)%contact_force = 0.0_real8
    system%particle(i)%friction_force = 0.0_real8
    system%particle(i)%adhesion_force = 0.0_real8
    system%particle(i)%friction_moment = 0.0_real8
    system%particle(i)%rolling_resistance_moment = 0.0_real8
    system%particle(i)%adhesion_moment = 0.0_real8
  END DO

  ! initialize particles´ force and moment vectors (forces and moments at the beginning of the step)
  DO i=1,system%no_particles
    i_kind = system%particle(i)%parameters%kind
    CALL detect_and_resolve_particle_contacts(i,Rt,Vt,Adelta,Wt,sys=system)
    CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
    CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt)
    system%particle(i)%nearfield_force = finft
    system%particle(i)%environment_force = fienvt
  END DO

  ! perform time stepping (one time step is named here a "substep")
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

     ! increment substep counter
     ss = ss + 1
     system%solution_cv%step(s)%current_substep = ss

     ! get particles´ arrays at time t
     DO i=1,system%no_particles
       Rt(i)%r = system%particle(i)%position
       Vt(i)%r = system%particle(i)%velocity
       At(i)%r = system%particle(i)%angle
       Wt(i)%r = system%particle(i)%spin
       Fcont(i)%r = system%particle(i)%contact_force
       Ffrict(i)%r = system%particle(i)%friction_force
       Fadht(i)%r = system%particle(i)%adhesion_force
       Fnft(i)%r = system%particle(i)%nearfield_force
       Fenvt(i)%r = system%particle(i)%environment_force
       Mfrict(i)%r = system%particle(i)%friction_moment
       Mrolt(i)%r = system%particle(i)%rolling_resistance_moment
       Madht(i)%r = system%particle(i)%adhesion_moment
       CALL get_particle_contact_history_arrays(i,system,PPCHt,PWCHt)
     END DO

     ! print substep information on screen
  10 WRITE(default_output_unit,*)
     WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

     ! set predictor for positions, velocities, angles and spins at time t+dt
     DO i=1,system%no_particles
       Rtpdt(i)%r = system%particle(i)%position + dt*system%particle(i)%velocity
       Vtpdt(i)%r = system%particle(i)%velocity
       Adelta(i)%r = dt*system%particle(i)%spin !0.0_real8
       Wtpdt(i)%r = system%particle(i)%spin
     END DO

    ! compute rigid walls positions, velocities and outside normals at time t+dt
    DO w=1,system%no_rigid_walls
      wkind = system%wall(w)%kind
      CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
      CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
    END DO

     ! perform iterations
     iterations: DO iter=1,max_niter
       
        ! save number of iterations performed
        niter = iter
       
        ! clear iterative contact-related quantities
        DO i=1,system%no_particles
          system%particle(i)%contact_force = 0.0_real8
          system%particle(i)%friction_force = 0.0_real8
          system%particle(i)%adhesion_force = 0.0_real8
          system%particle(i)%friction_moment = 0.0_real8
          system%particle(i)%rolling_resistance_moment = 0.0_real8
          system%particle(i)%adhesion_moment = 0.0_real8
          CALL reset_particle_contact_lists(i,PPCHt,PWCHt,system)
        END DO

        ! loop over particles
        DO i=1,system%no_particles
          i_kind = system%particle(i)%parameters%kind
          mi = system%particle(i)%mass
          ji = system%particle(i)%inertia
          rit = Rt(i)%r
          vit = Vt(i)%r
          ait = At(i)%r
          wit = Wt(i)%r
          finft = Fnft(i)%r
          fienvt = Fenvt(i)%r
          ficont = Fcont(i)%r
          fifrit = Ffrict(i)%r
          fiadht = Fadht(i)%r
          mifrit = Mfrict(i)%r
          mirolt = Mrolt(i)%r
          miadht = Madht(i)%r
          figivt = 0.0_real8
          migivt = 0.0_real8
          !IF (s==1 .AND. ss==1) figivt = system%particle(i)%initial_given_force
          !IF (s==1 .AND. ss==1) migivt = system%particle(i)%initial_given_moment
          CALL detect_and_resolve_particle_contacts(i,Rtpdt,Vtpdt,Adelta,Wtpdt,sys=system) 
          ficontpdt = system%particle(i)%contact_force
          fifritpdt = system%particle(i)%friction_force
          fiadhtpdt = system%particle(i)%adhesion_force
          mifritpdt = system%particle(i)%friction_moment
          miroltpdt = system%particle(i)%rolling_resistance_moment
          miadhtpdt = system%particle(i)%adhesion_moment
          CALL compute_particle_nearfield_forces(i,i_kind,Rtpdt,Vtpdt,system,finftpdt)
          CALL compute_particle_environment_forces(i,i_kind,Rtpdt,Vtpdt,system,fienvtpdt) 
          vitpdt = vit + (dt/mi)*(phi*(finftpdt+fienvtpdt) + (one-phi)*(finft+fienvt+figivt)) & 
                       + (dt/mi)*(phi*(ficontpdt+fifritpdt+fiadhtpdt) + (one-phi)*(ficont+fifrit+fiadht))
          IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
          ritpdt = rit + (phi*vitpdt + (one-phi)*vit)*dt
          IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
          witpdt = wit + (dt/ji)*(phi*(mifritpdt+miroltpdt+miadhtpdt) + (one-phi)*(mifrit+mirolt+miadht+migivt))
          IF (ANY(system%particle(i)%rotational_dof_codes(4:6)/=0)) CALL enforce_spin_constraints(i,system,witpdt)
          aidelta = (phi*witpdt + (one-phi)*wit)*dt
          IF (ANY(system%particle(i)%rotational_dof_codes(1:3)/=0)) CALL enforce_angle_constraints(i,system,aidelta,witpdt)
          Rtpdt_old(i)%r = Rtpdt(i)%r
          Rtpdt(i)%r = ritpdt
          Vtpdt_old(i)%r = Vtpdt(i)%r
          Vtpdt(i)%r = vitpdt
          Adelta_old(i)%r = Adelta(i)%r
          Adelta(i)%r = aidelta
          Wtpdt_old(i)%r = Wtpdt(i)%r
          Wtpdt(i)%r = witpdt
          Fcontpdt(i)%r = ficontpdt
          Ffrictpdt(i)%r = fifritpdt
          Fadhtpdt(i)%r = fiadhtpdt
          Fnftpdt(i)%r = finftpdt
          Fenvtpdt(i)%r = fienvtpdt
          Mfrictpdt(i)%r = mifritpdt
          Mroltpdt(i)%r = miroltpdt
          Madhtpdt(i)%r = miadhtpdt
          CALL update_particle_contact_history_arrays(i,system,PPCHtpdt,PWCHtpdt)
          Rresidual(i)%r = Rtpdt(i)%r - Rtpdt_old(i)%r
          Rnormalizer(i)%r = Rtpdt(i)%r - Rt(i)%r
          Vresidual(i)%r = Vtpdt(i)%r - Vtpdt_old(i)%r
          Vnormalizer(i)%r = Vtpdt(i)%r - Vt(i)%r
          Aresidual(i)%r = Adelta(i)%r - Adelta_old(i)%r
          Anormalizer(i)%r = Adelta(i)%r
          Wresidual(i)%r = Wtpdt(i)%r - Wtpdt_old(i)%r
          Wnormalizer(i)%r = Wtpdt(i)%r - Wt(i)%r
        END DO

        ! compute iteration errors and tolerances
        CALL compute_iteration_errors(Rresidual,Vresidual,Aresidual,Wresidual,Rnormalizer,Vnormalizer,Anormalizer,Wnormalizer, &
                                      Rtpdt,Vtpdt,Adelta,Wtpdt,Rerror,Verror,Aerror,Werror,tol1,tol2,tolR,tolV,tolA,tolW)

  	    ! print iteration information
  	    WRITE(default_output_unit,20) " iter=", iter, "errR=", Rerror, "errV=", Verror, "errA=", Aerror, "errW=", Werror
        20 FORMAT(TR2,A,I3,TR2,A,ES10.3,TR2,A,ES10.3,TR2,A,ES10.3,TR2,A,ES10.3)

        ! check for convergence
        IF ( (Rerror<=tolR .AND. Verror<=tolV) .AND. (Aerror<=tolA .AND. Werror<=tolW) ) EXIT iterations
        
        ! check for divergence 
        IF ( (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) .OR. (Aerror>=1.0E+06 .OR. Werror>=1.0E+06) ) EXIT iterations

     END DO iterations 

     ! stop solution or adapt time ste in case of nonconvergence or divergence
	 IF ( iter>max_niter .OR. (Rerror>=1.0E+06 .OR. Verror>=1.0E+06) .OR. (Aerror>=1.0E+06 .OR. Werror>=1.0E+06) ) THEN
	   WRITE(error_unit,"(A,I5,A,ES13.6)") "  Failure to converge at substep =", ss, ",   t+dt =", tpdt
	   WRITE(error_unit,*) "  Time step will be reduced if adaptivity=on or else program will be terminated"
	   WRITE(default_output_unit,"(A,I5,A,ES13.6)") " failure to converge at substep =", ss, ", t+dt =", tpdt
	   WRITE(default_output_unit,*) " time step will be reduced if adaptivity=on or else program will be terminated"
       IF (system%solution_cv%step(s)%adaptive_time_stepping=="on") THEN
         CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
	     tpdt = t + dt
         tpdt_real4 = REAL(tpdt,real4)
       
         ! ZZZZZZZZZZZZZZ essa parte talvez não seja necessária porque aqui se está REDUZINDO o dt e não aumentando... ZZZZZZZZZZZZZZZ
         IF (tpdt_real4>tprint_real4) THEN
           dt = tpdt - tprint
           tpdt = tprint
           tpdt_real4 = tprint_real4
         END IF
         IF (system%solution_cv%contact_detection_algorithm=="binning" .AND. tpdt_real4>tcellupdate_real4) THEN
           dt = tpdt - tcellupdate
           tpdt = tcellupdate
           tpdt_real4 = tcellupdate_real4
         END IF
         IF (system%solution_cv%contact_detection_algorithm=="verlet_list" .AND. tpdt_real4>tverletupdate_real4) THEN
           dt = tpdt - tverletupdate
           tpdt = tverletupdate
           tpdt_real4 = tverletupdate_real4
         END IF     
         !ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
         
	     system%solution_cv%step(s)%current_time = tpdt
         system%solution_cv%step(s)%current_dt = dt
         ! loop below is not needed: particles´ forces at t are already stored in Fcont, Ffrict, etc; besides, sys%part%cont_force (etc) is cleared as soon as iterations start!
         !DO i=1,system%no_particles
         !  system%particle(i)%contact_force = Fcont(i)%r
         !  system%particle(i)%friction_force = Ffrict(i)%r
         !  system%particle(i)%adhesion_force = Fadht(i)%r
         !  system%particle(i)%nearfield_force = Fnft(i)%r
         !  system%particle(i)%environment_force = Fenvt(i)%r
         !  system%particle(i)%friction_moment = Mfrict(i)%r
         !  system%particle(i)%rolling_resistance_moment = Mrolt(i)%r
         !  system%particle(i)%adhesion_moment = Madht(i)%r
         !END DO         
         GOTO 10
       END IF
	   EXIT steps
	 END IF

     ! update particles´ positions, velocities, angles, spins, forces, moments and contact lists
     DO i=1,system%no_particles
       system%particle(i)%position = Rtpdt(i)%r 
       system%particle(i)%velocity = Vtpdt(i)%r
       system%particle(i)%angle = (four/(four-DOT_PRODUCT(Adelta(i)%r,At(i)%r)))*(At(i)%r+Adelta(i)%r+half*(Adelta(i)%r.vector.At(i)%r))
       system%particle(i)%spin = Wtpdt(i)%r
       system%particle(i)%nearfield_force = Fnftpdt(i)%r
       system%particle(i)%environment_force = Fenvtpdt(i)%r
       system%particle(i)%contact_force = Fcontpdt(i)%r
       system%particle(i)%friction_force = Ffrictpdt(i)%r
       system%particle(i)%adhesion_force = Fadhtpdt(i)%r
       system%particle(i)%friction_moment = Mfrictpdt(i)%r
       system%particle(i)%rolling_resistance_moment = Mroltpdt(i)%r
       system%particle(i)%adhesion_moment = Madhtpdt(i)%r
       CALL update_particle_contact_lists(i,PPCHtpdt,PWCHtpdt,system)
     END DO

     ! print substep results in results file
     IF (tpdt_real4==tprint_real4) THEN 
       IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
       CALL print_system_results(system)
       tprint = tprint + dt_print
       tprint_real4 = REAL(tprint,real4)
     END IF

     ! update cell and verlet lists
     IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
       IF (tpdt_real4==tcellupdate_real4) THEN 
         CALL update_grid_cells_lists(system)
         tcellupdate = tcellupdate + dt_cell_list_update
         tcellupdate_real4 = REAL(tcellupdate,real4)
       END IF
     END IF
     IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
       IF (tpdt_real4==tverletupdate_real4) THEN
         CALL update_grid_cells_lists(system)
         CALL update_particles_verlet_lists(system)
         tverletupdate = tverletupdate + dt_verlet_list_update
         tverletupdate_real4 = REAL(tverletupdate,real4)
       END IF
     END IF

     ! compute dt for next substep (adaptive time stepping)
     IF (system%solution_cv%step(s)%adaptive_time_stepping=="on" .AND. tpdt<system%solution_cv%step(s)%final_time) THEN
       CALL adapt_time_step_size(dt,dt_min,dt_max,niter,desired_niter)
     END IF

     ! increment time for next substep
     t = tpdt
     tpdt = tpdt + dt
     tpdt_real4 = REAL(tpdt,real4)

     ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)     
     IF (tpdt_real4==tprint_real4) THEN
       tpdt = tprint
     END IF
     
     ! check for consistency in new time value (check if new time exceeds tprint or tcell or tverlet)
     IF (tpdt_real4>tprint_real4) THEN
       dt = tprint - t  !tpdt - tprint
       tpdt = tprint
       tpdt_real4 = tprint_real4
     END IF
     IF (system%solution_cv%contact_detection_algorithm=="binning" .AND. tpdt_real4>tcellupdate_real4) THEN
       dt = tcellupdate - t  !tpdt - tcellupdate
       tpdt = tcellupdate
       tpdt_real4 = tcellupdate_real4
     END IF
     IF (system%solution_cv%contact_detection_algorithm=="verlet_list" .AND. tpdt_real4>tverletupdate_real4) THEN
       dt = tverletupdate - t  !tpdt - tverletupdate
       tpdt = tverletupdate
       tpdt_real4 = tverletupdate_real4
     END IF     

     ! update time control variables
     system%solution_cv%step(s)%current_time = tpdt
     system%solution_cv%step(s)%current_dt = dt
     
  END DO substeps

END DO steps

! deallocate arrays
DEALLOCATE(Rt,Rtpdt,Rtpdt_old,STAT=derror_R)
DEALLOCATE(Vt,Vtpdt,Vtpdt_old,STAT=derror_V)
DEALLOCATE(At,Adelta,STAT=derror_A)
DEALLOCATE(Wt,Wtpdt,STAT=derror_W)
DEALLOCATE(Fcont,Fcontpdt,Ffrict,Ffrictpdt,Fadht,Fadhtpdt,Fnft,Fnftpdt,Fenvt,Fenvtpdt,STAT=derror_F)
DEALLOCATE(Mfrict,Mfrictpdt,Mrolt,Mroltpdt,Madht,Madhtpdt,STAT=derror_M)
DEALLOCATE(Rresidual,Vresidual,Rnormalizer,Vnormalizer,STAT=derror_Rres_Vres_and_normalizers)
DEALLOCATE(Aresidual,Wresidual,Anormalizer,Wnormalizer,STAT=derror_Ares_Wres_and_normalizers)
CALL deallocate_contact_history_arrays(system%no_particles,PPCHt,PWCHt,PPCHtpdt,PWCHtpdt)

! check for deallocation error
derror_R = derror_R + derror_V + derror_A + derror_W + derror_F + derror_M + derror_Rres_Vres_and_normalizers + derror_Ares_Wres_and_normalizers 
IF (derror_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) and arrays at the end of midpoint_implicit solver"


END SUBROUTINE midpoint_implicit_solver_for_soft_spheres
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE newmark_implicit_solver(system)

! dummy arguments
TYPE(system_data) :: system


WRITE(error_unit,*) " Newmark solver not yet implemented!"
WRITE(default_output_unit,*) " Newmark solver not yet implemented!"
STOP

END SUBROUTINE newmark_implicit_solver
!!---------------------------------------------------------------------------------------------------------------------


END MODULE mechanics_solver_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE thermomechanics_solver_class

USE particles_classes
USE rigid_walls_classes
USE external_heating_devices_class
USE periodic_boundary_conditions
USE grid_subroutines
USE verlet_lists_subroutines
USE contact_history_array_class
USE contacts_solver
USE dof_constraints
USE system_properties
USE print_results_for_thermomechanics_solver
USE solvers_utilities

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE thermomechanics_problem_solver(system)

! dummy arguments
TYPE(system_data) :: system


! select time integration solver
SELECT CASE (system%solution_cv%solver_type)

  CASE("euler_explicit_solver_for_soft_spheres")
    CALL euler_explicit_solver_for_soft_spheres(system)
  
  !CASE("midpoint_implicit_solver_for_soft_spheres")
  !  CALL midpoint_implicit_solver_for_soft_spheres(system)

  CASE DEFAULT
    WRITE(error_unit,*) "Solver type not recognized or not yet implemented"
    STOP

END SELECT

END SUBROUTINE thermomechanics_problem_solver
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE euler_explicit_solver_for_soft_spheres(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
CHARACTER(50)      :: i_kind, wkind, dkind
INTEGER(KIND=int4) :: s, ss, i, j, w, k, derror_R, derror_V, derror_A, derror_W, derror_T, derror_F
REAL(KIND=real8)   :: t, tpdt, tprint, tcellupdate, tverletupdate, dt, dt_print, dt_cell_list_update, dt_verlet_list_update, &
                      mi, ji, ci, tit, titpdt, qpienvt, qpicondt, qpiextt, Ref_norm
REAL(KIND=real4)   :: tpdt_real4, tfinal_real4, tprint_real4, tcellupdate_real4, tverletupdate_real4
REAL(KIND=real8), DIMENSION(3) :: rit, ritpdt, vit, vitpdt, ait, aidelta, wit, witpdt, ficont, fifrit, finft, &
                                  fienvt, mifrit, mirolt, figivt, migivt, fiadht, miadht, mienvt
TYPE(array_of_vectors), DIMENSION(:), POINTER :: Rt, Rtpdt, Vt, Vtpdt, At, Adelta, Wt, Wtpdt, Tt, Ttpdt, &
                                                 Fnft, Fenvt, Fcont, Ffrict, Fadht, Mfrict, Mrolt, Madht, Menvt, QPenvt, QPcondt, QPextt 
TYPE(contact_history_array), DIMENSION(:), ALLOCATABLE :: PPCHt, PWCHt
REAL(KIND=real8), PARAMETER :: one = 1.0_real8, four = 4.0_real8, half = 0.5_real8

! Description of variables:
! (vectors beginning with capital letters refer to a whole-system vector; with non-capitals to a single particle vector)
! Rt = vector of position vectors at time t (contains the position vectors of all particles at time t)         
! Rtpdt = vector of position vectors at time t+dt (contains the position vectors of all particles at time t+dt)         
! Vt = vector of velocity vectors at time t (contains the velocity vectors of all particles at time t)
! Vtpdt = vector of velocity vectors at time t+dt (contains the velocity vectors of all particles at time t+dt)         
! Fnft = vector of near-field force vectors at time t (contains the external force vectors of all particles at time t)
! QPenvt = vector of environment heat powers at time t (contains the environment heat powers of all particles at time t)
! rit = position vector of particle "i" at time t
! ritpdt = position vector of particle "i" at time t+dt
! vit = velocity vector of particle "i" at time t
! vitpdt = velocity vector of particle "i" at time t+dt
! finft = near-field force vector of particle "i" at time t
! PPCHt = array of particle-particle contact histories at time t
! PWCHt = array of particle-wall contact histories at time t
! etc

! check for for proper solver type
IF (system%general_cv%rotational_dofs=="off") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match rotational dofs switch. Program execution aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match rotational dofs switch. Program execution aborted. "
  STOP
END IF
IF (system%general_cv%temperature_dofs=="off") THEN
  WRITE(error_unit,*) " ERROR: solver type does not match temperature dofs switch. Program execution aborted. "
  WRITE(default_output_unit,*) " ERROR: solver type does not match temperature dofs switch. Program execution aborted. "
  STOP
END IF

! allocate arrays
ALLOCATE(Rt(1:system%no_particles),Rtpdt(1:system%no_particles))
ALLOCATE(Vt(1:system%no_particles),Vtpdt(1:system%no_particles))
ALLOCATE(At(1:system%no_particles),Adelta(1:system%no_particles))
ALLOCATE(Wt(1:system%no_particles),Wtpdt(1:system%no_particles))
ALLOCATE(Tt(1:system%no_particles),Ttpdt(1:system%no_particles))
ALLOCATE(Fnft(1:system%no_particles),Fenvt(1:system%no_particles))
ALLOCATE(Fcont(1:system%no_particles),Ffrict(1:system%no_particles))
ALLOCATE(Mfrict(1:system%no_particles),Mrolt(1:system%no_particles))
ALLOCATE(Fadht(1:system%no_particles),Madht(1:system%no_particles))
ALLOCATE(Menvt(1:system%no_particles))
ALLOCATE(QPenvt(1:system%no_particles),QPcondt(1:system%no_particles),QPextt(1:system%no_particles))
ALLOCATE(PPCHt(1:system%no_particles),PWCHt(1:system%no_particles))

! clear time values
t = 0.0_real8
tpdt = 0.0_real8
tprint = 0.0_real8
tcellupdate = 0.0_real8
tverletupdate = 0.0_real8

! get particles´ initial conditions and save at corresponding arrays
DO i=1,system%no_particles
  system%particle(i)%position = system%particle(i)%xyz_coordinates
  system%particle(i)%velocity = system%particle(i)%initial_velocity
  system%particle(i)%angle = system%particle(i)%initial_angles
  system%particle(i)%spin = system%particle(i)%initial_spin
  system%particle(i)%temperature = system%particle(i)%initial_temperature
  Rt(i)%r = system%particle(i)%position
  Vt(i)%r = system%particle(i)%velocity
  At(i)%r = system%particle(i)%angle
  Wt(i)%r = system%particle(i)%spin
  Tt(i)%r = system%particle(i)%temperature
  Adelta(i)%r = 0.0_real8
END DO

! get rigid walls´ initial conditions
DO w=1,system%no_rigid_walls
  system%wall(w)%point_position = system%wall(w)%initial_point_position
  system%wall(w)%velocity = system%wall(w)%initial_velocity
  system%wall(w)%outside_normal = normalize(system%wall(w)%initial_outside_normal)
  system%wall(w)%temperature = system%wall(w)%initial_temperature
END DO

! get external heating devices´ initial conditions
DO k=1,system%no_external_heating_devices
  system%external_heating_device(k)%striking_position = system%external_heating_device(k)%initial_striking_position
  system%external_heating_device(k)%velocity = system%external_heating_device(k)%initial_velocity
  system%external_heating_device(k)%intensity = system%external_heating_device(k)%initial_intensity
END DO

! compute particles´ mass, inertia and thermal properties
DO i=1,system%no_particles
  i_kind = system%particle(i)%parameters%kind
  CALL compute_particle_mass_and_inertia(i,i_kind,system)
  CALL compute_particle_thermal_properties(i,i_kind,system)
END DO

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! THINK: maybe the grid, verlet lists and contact lists below shall be initialized within the steps!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! initialize grid and verlet lists
IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
  CALL sort_particles_into_grid_cells(system)
END IF
IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
  CALL sort_particles_into_grid_cells(system)
  CALL build_particles_verlet_lists(system)
END IF

! initialize particles´ contact lists
DO i=1,system%no_particles
  CALL initialize_particle_contact_lists(i,system)
END DO  

! initialize rigid walls´ contact lists
DO w=1,system%no_rigid_walls
  CALL initialize_wall_contact_lists(w,system)
END DO  

! loop over steps (one step is one "load case")
steps: DO s=1,system%solution_cv%no_steps

  ! set time increments
  dt = system%solution_cv%step(s)%initial_dt
  dt_print = system%general_cv%dt_for_results_printing
  dt_cell_list_update = system%grid%dt_for_cell_list_update
  dt_verlet_list_update = system%general_cv%dt_for_verlet_list_update  

  ! initialize time values (tpdt is t+dt, i.e. is the time at the end of the substep)
  t = tpdt
  tpdt = tpdt + dt
  tprint = tprint + dt_print
  tcellupdate = tcellupdate + dt_cell_list_update
  tverletupdate = tverletupdate + dt_verlet_list_update
  
  ! transform time values to single precision (to kill round-off errors)
  tpdt_real4 = REAL(tpdt,real4)
  tprint_real4 = REAL(tprint,real4)
  tcellupdate_real4 = REAL(tverletupdate,real4)
  tverletupdate_real4 = REAL(tverletupdate,real4)
  tfinal_real4 = REAL(system%solution_cv%step(s)%final_time,real4)

  ! initialize step control variables
  system%solution_cv%current_step = s
  system%solution_cv%step(s)%current_substep = 0
  system%solution_cv%step(s)%current_time = tpdt
  system%solution_cv%step(s)%current_dt = dt

  ! initialize substep counter
  ss = 0
    
  ! perform time stepping (one time step is named a substep here)
  substeps: DO WHILE (tpdt_real4<=tfinal_real4)

     ! increment substep counter
     ss = ss + 1
     system%solution_cv%step(s)%current_substep = ss
     
     ! get particles´ arrays at time t
     DO i=1,system%no_particles
       Rt(i)%r = system%particle(i)%position
       Vt(i)%r = system%particle(i)%velocity
       At(i)%r = system%particle(i)%angle
       Wt(i)%r = system%particle(i)%spin
       Tt(i)%r = system%particle(i)%temperature
     END DO 
     
     ! clear particles´ contact-related quantities (forces, moments and heat power)
     DO i=1,system%no_particles
       system%particle(i)%contact_force = 0.0_real8
       system%particle(i)%friction_force = 0.0_real8
       system%particle(i)%adhesion_force = 0.0_real8
       system%particle(i)%friction_moment = 0.0_real8
       system%particle(i)%rolling_resistance_moment = 0.0_real8
       system%particle(i)%adhesion_moment = 0.0_real8
       system%particle(i)%conduction_heat_power = 0.0_real8
     END DO

     ! clear walls´ contact-related quantities
     DO w=1,system%no_rigid_walls
       system%wall(w)%contact_force = 0.0_real8
       system%wall(w)%friction_force = 0.0_real8
       !system%wall(w)%conduction_heat_power = 0.0_real8
     END DO
     
     ! compute particles´ contact-related quantites (forces, moments and heat power) at time t and save at corresponding arrays
     DO i=1,system%no_particles
       CALL detect_and_resolve_particle_contacts(i,Rt,Vt,Adelta,Wt,Tt,system)
       ! notice: at this point, the walls´ contact lists and total contact forces have been fully computed!
       Fcont(i)%r = system%particle(i)%contact_force
       Ffrict(i)%r = system%particle(i)%friction_force
       Fadht(i)%r = system%particle(i)%adhesion_force
       Mfrict(i)%r = system%particle(i)%friction_moment
       Mrolt(i)%r = system%particle(i)%rolling_resistance_moment
       Madht(i)%r = system%particle(i)%adhesion_moment
       QPcondt(i)%r = system%particle(i)%conduction_heat_power
       CALL update_particle_contact_history_arrays(i,system,PPCHt,PWCHt)
     END DO

     ! compute particles´ nearfield and environment forces, moments and heat powers at time t and save at corresponding arrays
     DO i=1,system%no_particles
       i_kind = system%particle(i)%parameters%kind
       CALL compute_particle_nearfield_forces(i,i_kind,Rt,Vt,system,finft)
       CALL compute_particle_environment_forces(i,i_kind,Rt,Vt,system,fienvt,mienvt,Wt)
       CALL compute_particle_environment_heat_power(i,i_kind,Rt,Vt,Tt,system,qpienvt)
       CALL compute_particle_external_device_heat_power(i,i_kind,Rt,Vt,system,qpiextt)
       Fnft(i)%r  = finft
       Fenvt(i)%r = fienvt
       Menvt(i)%r = mienvt
       QPenvt(i)%r(1) = qpienvt
       QPextt(i)%r(1) = qpiextt
     END DO

     ! enforce periodic boundary conditions
     IF (system%general_cv%periodic_bc=="on") THEN
       CALL enforce_periodic_boundary_conditions(Rt,system)
     END IF
     
     ! print substep information on screen
     WRITE(default_output_unit,*)
     WRITE(default_output_unit,"(A,I2,TR6,A,I9,TR6,A,ES13.6)") " step =", s, " substep =", ss, " t+dt =", tpdt

    ! loop over particles
    DO i=1,system%no_particles
      mi = system%particle(i)%mass
      ji = system%particle(i)%inertia
      ci = system%particle(i)%specific_heat
      rit = Rt(i)%r
      vit = Vt(i)%r
      ait = At(i)%r
      wit = Wt(i)%r
      tit = Tt(i)%r(1)
      finft  = Fnft(i)%r
      fienvt = Fenvt(i)%r
      ficont = Fcont(i)%r
      fifrit = Ffrict(i)%r
      fiadht = Fadht(i)%r
      mifrit = Mfrict(i)%r
      mirolt = Mrolt(i)%r
      miadht = Madht(i)%r
      mienvt = Menvt(i)%r
      qpienvt = QPenvt(i)%r(1)
      qpicondt = QPcondt(i)%r(1)
      qpiextt = QPextt(i)%r(1)
      figivt = 0.0_real8
      migivt = 0.0_real8
      !IF (s==1 .AND. ss==1) figivt = system%particle(i)%initial_given_force
      !IF (s==1 .AND. ss==1) migivt = system%particle(i)%initial_given_moment
      vitpdt = vit + (dt/mi)*(finft+fienvt+figivt) + (dt/mi)*(ficont+fifrit+fiadht)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !IF (SQRT(DOT_PRODUCT(vitpdt,vitpdt))>=0.1_real8) THEN
      !  vitpdt = 0.0_real8
      !END IF  
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF (ANY(system%particle(i)%translational_dof_codes(4:6)/=0)) CALL enforce_velocity_constraints(i,system,vitpdt)
      ritpdt = rit + dt*vitpdt
      IF (ANY(system%particle(i)%translational_dof_codes(1:3)/=0)) CALL enforce_position_constraints(i,system,ritpdt,vitpdt)
      witpdt = wit + (dt/ji)*(mifrit+mirolt+miadht+mienvt+migivt)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !IF (SQRT(DOT_PRODUCT(witpdt,witpdt))>=0.1_real8) THEN
      !  witpdt = 0.0_real8
      !END IF  
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF (ANY(system%particle(i)%rotational_dof_codes(4:6)/=0)) CALL enforce_spin_constraints(i,system,witpdt)
      aidelta = dt*witpdt
      IF (ANY(system%particle(i)%rotational_dof_codes(1:3)/=0)) CALL enforce_angle_constraints(i,system,aidelta,witpdt)
      titpdt = tit + (dt/(mi*ci))*(qpienvt+qpicondt+qpiextt)
      IF (system%particle(i)%temperature_dof_codes/=0) CALL enforce_temperature_constraints(i,system,titpdt)
      Rtpdt(i)%r = ritpdt
      Vtpdt(i)%r = vitpdt
      Adelta(i)%r = aidelta
      Wtpdt(i)%r = witpdt
      Ttpdt(i)%r(1) = titpdt
    END DO
       
    ! check for numerical instability (not a rigorous check)
	!V_norm = norm(Vtpdt)
	!Ref_norm = norm(Vt)
	!IF (Ref_norm<=1.0E-10) Ref_norm = 1.0_real8
	!Relative_velocity = V_norm/Ref_norm
	!IF (Relative_velocity>=1.0E+08) THEN
	!   WRITE(error_unit,*) " CAUTION: large increase on velocities at substep =", ss, ", t+dt =", REAL(tpdt,real4)
	!   WRITE(error_unit,*) " Be aware of possible instabilities within the explicit solver"
	!   WRITE(default_output_unit,*) " CAUTION: large increase on velocities at substep =", ss, ", t+dt =", REAL(tpdt,real4)
	!   WRITE(default_output_unit,*) " Be aware of possible instabilities within the explicit solver"
	!END IF

    ! update particles´ positions, velocities, angles, spins and temperatures
    DO i=1,system%no_particles
      system%particle(i)%position = Rtpdt(i)%r 
      system%particle(i)%velocity = Vtpdt(i)%r
      system%particle(i)%angle = (four/(four-DOT_PRODUCT(Adelta(i)%r,At(i)%r)))*(At(i)%r+Adelta(i)%r+half*(Adelta(i)%r.vector.At(i)%r))
      system%particle(i)%spin = Wtpdt(i)%r
      system%particle(i)%temperature = Ttpdt(i)%r(1)
    END DO

    ! print substep results in results file
    IF (tpdt_real4==tprint_real4) THEN 
      IF (system%general_cv%compute_system_properties=="yes") CALL compute_system_properties(system)
      CALL print_system_thermomechanics_results(system)
      tprint = tprint + dt_print
      tprint_real4 = REAL(tprint,real4)
    END IF

    ! update cell and verlet lists
    IF (system%solution_cv%contact_detection_algorithm=="binning") THEN
      IF (tpdt_real4==tcellupdate_real4) THEN 
        CALL update_grid_cells_lists(system)
        tcellupdate = tcellupdate + dt_cell_list_update
        tcellupdate_real4 = REAL(tcellupdate,real4)
      END IF
    END IF
    IF (system%solution_cv%contact_detection_algorithm=="verlet_list") THEN
      IF (tpdt_real4==tverletupdate_real4) THEN
        CALL update_grid_cells_lists(system)
        CALL update_particles_verlet_lists(system)
        tverletupdate = tverletupdate + dt_verlet_list_update
        tverletupdate_real4 = REAL(tverletupdate,real4)
      END IF
    END IF

    ! update rigid walls´ positions, velocities, outside normals and temperatures (values at t+dt, i.e., for the next substep)    
    DO w=1,system%no_rigid_walls
      wkind = system%wall(w)%kind
      CALL compute_rigid_wall_position_and_velocity(w,wkind,tpdt,dt,system)
      CALL compute_rigid_wall_outside_normal(w,wkind,tpdt,dt,system)
      CALL compute_rigid_wall_temperature_and_heating_rate(w,wkind,tpdt,dt,system)
    END DO

    ! update external heating devices´ positions, velocities and intensities (values at t+dt, i.e., for the next substep)    
    DO k=1,system%no_external_heating_devices
      dkind = system%external_heating_device(k)%kind
      CALL compute_external_heating_device_position_and_velocity(k,dkind,tpdt,dt,system)
      CALL compute_external_heating_device_intensity(k,dkind,tpdt,dt,system)
    END DO

    ! increment time for next substep
    t = tpdt
    tpdt = tpdt + dt
    tpdt_real4 = REAL(tpdt,real4)

    ! reduce round-off errors on accumulated time (especially useful when using time adaptivity)
    IF (tpdt_real4==tprint_real4) THEN
      tpdt = tprint
    END IF

    ! update time control variable
    system%solution_cv%step(s)%current_time = tpdt

  END DO substeps

END DO steps

! deallocate arrays
DEALLOCATE(Rt,Rtpdt,STAT=derror_R)
DEALLOCATE(Vt,Vtpdt,STAT=derror_V)
DEALLOCATE(At,Adelta,STAT=derror_A)
DEALLOCATE(Wt,Wtpdt,STAT=derror_W)
DEALLOCATE(Tt,Ttpdt,STAT=derror_T)
DEALLOCATE(Fnft,Fenvt,Fcont,Ffrict,Fadht,Mfrict,Mrolt,Madht,Menvt,QPenvt,QPcondt,QPextt,STAT=derror_F)
CALL deallocate_contact_history_arrays(system%no_particles,PPCHt,PWCHt)

! check for deallocation error
derror_R = derror_R + derror_V + derror_A + derror_W + derror_T + derror_F
IF (derror_R/=0) WRITE(error_unit,*) "deallocation failure of vector(s) at the end of the explicit solver"

END SUBROUTINE euler_explicit_solver_for_soft_spheres
!!---------------------------------------------------------------------------------------------------------------------


END MODULE thermomechanics_solver_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE solver_class

USE system_data_types
USE mechanics_solver_class
USE thermomechanics_solver_class
! USE thermics_solver_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE solver(system)

! dummy arguments
TYPE(system_data) :: system

! local variables
INTEGER(KIND=int4) :: i


problem_type: SELECT CASE (system%solution_cv%problem_type)

  CASE("mechanics")
    CALL mechanics_problem_solver(system)

  CASE("thermomechanics")
    CALL thermomechanics_problem_solver(system)

!  !CASE("thermics")
!    !CALL thermics_problem_solver(system)

  CASE DEFAULT
    WRITE(error_unit,*) "Problem type not recognized or not yet implemented"
    STOP

END SELECT problem_type

END SUBROUTINE solver
!!---------------------------------------------------------------------------------------------------------------------


END MODULE solver_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE input_data_class

USE particle_data_types
USE system_data_types
USE pefmat

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE capture_keyword(keyword)

! dummy arguments
CHARACTER(LEN=*), INTENT(IN) :: keyword

! local variables
CHARACTER(LEN=LEN(keyword)) :: key

READ(input_file_unit,*) key
IF(ADJUSTL(key)/=keyword)THEN
  WRITE(error_unit,*) "keyword ", TRIM(keyword), " not found"
  STOP
END IF

END SUBROUTINE capture_keyword
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE read_system_attribute(name,ni,nr,i,r)

CHARACTER(LEN=*) :: name
INTEGER(KIND=int4) :: ni, nr
INTEGER(KIND=int4), DIMENSION(ni) :: i
REAL(KIND=real4),   DIMENSION(nr) :: r

INTEGER(KIND=int4) :: row


IF(ni==0.AND.nr>0)THEN
  READ(input_file_unit,*)row, name, r
ELSE IF(ni>0.AND.nr>0)THEN
  READ(input_file_unit,*)row, name, i, r
ELSE IF(ni>0.AND.nr==0)THEN
  READ(input_file_unit,*)row, name, i
ELSE IF(ni==0.AND.nr==0)THEN
  READ(input_file_unit,*)row, name
END IF

END SUBROUTINE read_system_attribute
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE read_input_system_data(system)

! dummy arguments
TYPE(system_data), INTENT(OUT) :: system

! local variables
INTEGER(KIND=int4) :: i, row, r, ls
REAL(KIND=real4) :: kij_attr, kij_rep, expij_attr, expij_rep, dcutoff, deq, dinst


! read system name
CALL capture_keyword("$system_name")
READ(input_file_unit,*) system%name

! read number of particles
CALL capture_keyword("$no_particles")
READ(input_file_unit,*) system%no_particles
ALLOCATE(system%particle(1:system%no_particles))

! rotational_dofs switch (on/off)
CALL capture_keyword("$rotational_dofs")
READ(input_file_unit,*) system%general_cv%rotational_dofs

! temperature_dofs switch (on/off)
CALL capture_keyword("$temperature_dofs")
READ(input_file_unit,*) system%general_cv%temperature_dofs
WRITE(error_unit,*) "CAUTION: make sure input data is compatible with the temperature dofs switch"

! read particles´ shape
CALL capture_keyword("$particles_shape")
READ(input_file_unit,*) system%solution_cv%particles_shape

! read particles data (type, radius, charge and attributes)
CALL capture_keyword("$particles_attributes_and_data")
DO i=1,system%no_particles
  READ(input_file_unit,*) row, system%particle(i)%parameters%kind, system%particle(i)%radius, system%particle(i)%charge, &
                          system%particle(i)%material_set_number, system%particle(i)%pressure_surface_number
END DO                                                                                                         

! read number of material sets
CALL capture_keyword("$no_material_sets")
READ(input_file_unit,*) system%no_material_sets
ALLOCATE(system%material_set(1:system%no_material_sets))

! read material sets
CALL capture_keyword("$material_sets_data")
DO i=1,system%no_material_sets
  READ(input_file_unit,*) row, system%material_set(i)%mass_density, system%material_set(i)%charge_density, &
                          system%material_set(i)%elasticity_modulus, system%material_set(i)%poisson_coeff, &
                          system%material_set(i)%contact_damping_ratio, system%material_set(i)%friction_damping_ratio, &
                          system%material_set(i)%coefficient_of_restitution, system%material_set(i)%static_friction_coeff, &
                          system%material_set(i)%dynamic_friction_coeff, system%material_set(i)%rolling_resistance_coeff, &
                          system%material_set(i)%rolling_resistance_damping_ratio 
END DO

IF (system%general_cv%temperature_dofs=="on") THEN
  WRITE(error_unit,*) "CAUTION: make sure all thermal properties are provided for the material(s) and wall(s)"
  CALL capture_keyword("$material_sets_thermal_properties")
  DO i=1,system%no_material_sets
    READ(input_file_unit,*) row, system%material_set(i)%specific_heat, system%material_set(i)%thermal_conductivity,  &
                            system%material_set(i)%drag_heating_efficiency, system%material_set(i)%radiative_efficiency, &
                            system%material_set(i)%absorptance, system%material_set(i)%degrading_temperature
    IF (system%material_set(i)%specific_heat==0.0) THEN
      WRITE(error_unit,*) "ERROR: specific heat is zero in material set", i, "; program will be terminated"
      STOP
    END IF  
  END DO
END IF

! read particle coordinates and initial velocities
CALL capture_keyword("$particle_coordinates_and_initial_velocities")
DO i=1,system%no_particles
  IF (system%general_cv%rotational_dofs=="on") THEN
    READ(input_file_unit,*) row, system%particle(i)%xyz_coordinates, &
                            system%particle(i)%initial_velocity, system%particle(i)%initial_spin
    system%particle(i)%initial_angles = 0.0_real8
  ELSE
    READ(input_file_unit,*) row, system%particle(i)%xyz_coordinates, system%particle(i)%initial_velocity
  END IF
END DO

! read particle initial temperatures
IF (system%general_cv%temperature_dofs=="on") THEN
  CALL capture_keyword("$particle_initial_temperatures")
  DO i=1,system%no_particles
    READ(input_file_unit,*) row, system%particle(i)%initial_temperature
  END DO
END IF

! initialize particles´ dof codes
DO i=1,system%no_particles
  system%particle(i)%translational_dof_codes = 0
  IF (system%general_cv%rotational_dofs=="on") system%particle(i)%rotational_dof_codes = 0
END DO

! read number of constrained particles
CALL capture_keyword("$no_constrained_particles")
READ(input_file_unit,*) system%no_constrained_particles

! read dof codes of constrained particles
IF (system%no_constrained_particles>0) THEN
  CALL capture_keyword("$constrained_particles_dof_codes")
  DO i=1,system%no_constrained_particles
    IF (system%general_cv%rotational_dofs=="on") THEN
      IF (system%general_cv%temperature_dofs=="on") THEN
        READ(input_file_unit,*) row, system%particle(row)%translational_dof_codes, system%particle(row)%rotational_dof_codes, &
                                     system%particle(row)%temperature_dof_codes
      ELSE
        READ(input_file_unit,*) row, system%particle(row)%translational_dof_codes, system%particle(row)%rotational_dof_codes
      END IF  
    ELSE
      IF (system%general_cv%temperature_dofs=="on") THEN
        READ(input_file_unit,*) row, system%particle(row)%translational_dof_codes, system%particle(row)%temperature_dof_codes
      ELSE
        READ(input_file_unit,*) row, system%particle(row)%translational_dof_codes
      END IF    
    END IF
  END DO                                                                                                         
END IF

! read number of constrained particles with harmonic constraints
IF (system%no_constrained_particles>0) THEN
  CALL capture_keyword("$no_constrained_particles_with_harmonic_constraints")
  READ(input_file_unit,*) system%no_constrained_particles_with_harmonic_constraints
END IF

! read harmonic constraints data
IF (system%no_constrained_particles_with_harmonic_constraints>0) THEN
  CALL capture_keyword("$harmonic_constraints_data")
  DO i=1,system%no_constrained_particles_with_harmonic_constraints
    READ(input_file_unit,*) row, system%particle(row)%harmonic_constraints_data
  END DO
END IF

! read number of particles with initial given forces and moments (particles with given forces and moments at t=0)
CALL capture_keyword("$no_particles_with_initial_given_forces_and_moments")
READ(input_file_unit,*) system%no_particles_with_initial_given_forces_and_moments

! read initial given forces and moments (given forces and moments at t=0)
IF (system%no_particles_with_initial_given_forces_and_moments>0) THEN
  CALL capture_keyword("$particle_initial_given_forces_and_moments")
  DO i=1,system%no_particles_with_initial_given_forces_and_moments
    READ(input_file_unit,*) row, system%particle(row)%initial_given_force, system%particle(row)%initial_given_moment 
  END DO                                                                                                         
ELSE
  DO i=1,system%no_particles
    system%particle(i)%initial_given_force = 0.0_real8
    system%particle(i)%initial_given_moment = 0.0_real8
  END DO
END IF

! read number of rigid walls
CALL capture_keyword("$no_rigid_walls")
READ(input_file_unit,*) system%no_rigid_walls

! read rigid walls properties
IF (system%no_rigid_walls>0) THEN
  ALLOCATE(system%wall(1:system%no_rigid_walls))
  CALL capture_keyword("$rigid_walls_properties")
  DO i=1,system%no_rigid_walls
    READ(input_file_unit,*) row, system%wall(i)%kind, system%wall(i)%charge, system%wall(i)%contact_damping_ratio, &
                            system%wall(i)%friction_damping_ratio, system%wall(i)%coefficient_of_restitution,  &
                            system%wall(i)%static_friction_coeff, system%wall(i)%dynamic_friction_coeff,  &
                            system%wall(i)%rolling_resistance_coeff, system%wall(i)%rolling_resistance_damping_ratio
  END DO                                                                                                         
END IF

! read rigid walls geometry data
IF (system%no_rigid_walls>0) THEN
  CALL capture_keyword("$rigid_walls_geometry_and_velocity_data")
    DO i=1,system%no_rigid_walls
      IF (system%wall(i)%kind=="flat_rigid_wall") THEN
        READ(input_file_unit,*) row, system%wall(i)%initial_point_position, system%wall(i)%initial_outside_normal, &
                                system%wall(i)%initial_velocity
      ELSE IF (system%wall(i)%kind=="harmonic_flat_rigid_wall") THEN
        READ(input_file_unit,*) row, system%wall(i)%initial_point_position, system%wall(i)%initial_outside_normal, &
                                system%wall(i)%harmonic_amplitudes, system%wall(i)%harmonic_frequencies
      ELSE IF (system%wall(i)%kind=="forced_flat_rigid_wall") THEN
        WRITE (error_unit,*) "CAUTION: please don´t forget to provide initial acceleration and driving force for forced walls"
        WRITE (error_unit,*) "REMEMBER: initial acceleration and driving force must be parallel to forced wall´s outside normal!"
        READ(input_file_unit,*) row, system%wall(i)%initial_point_position, system%wall(i)%initial_outside_normal, &
                                system%wall(i)%initial_velocity, system%wall(i)%initial_acceleration, system%wall(i)%initial_driving_force 
      ELSE IF (system%wall(i)%kind=="cylindrical_rigid_wall") THEN
        READ(input_file_unit,*) row, system%wall(i)%initial_point_position, system%wall(i)%radius, &
                                system%wall(i)%initial_velocity
      ELSE  
        WRITE(error_unit,*) "rigid wall kind not recognized or not yet implemented"
        STOP
      END IF
    END DO
END IF

! read rigid walls thermal data
IF (system%no_rigid_walls>0 .AND. system%general_cv%temperature_dofs=="on") THEN
  CALL capture_keyword("$rigid_walls_thermal_properties")
  DO i=1,system%no_rigid_walls
    READ(input_file_unit,*) row, system%wall(i)%thermal_kind, system%wall(i)%initial_temperature, system%wall(i)%initial_heating_rate
  END DO
END IF

! read external fields vectors
CALL capture_keyword("$external_fields_vectors")
READ(input_file_unit,*) system%external_fields%gravity_accel_vector, &
                        system%external_fields%electric_field_vector, system%external_fields%magnetic_field_vector

! read electric and magnetic fields domain limits
IF (ANY(system%external_fields%electric_field_vector(1:3)/=0.0) .OR. ANY(system%external_fields%magnetic_field_vector(1:3)/=0.0)) THEN
  CALL capture_keyword("$electric_and_magnetic_fields_domain_limits")
  READ(input_file_unit,*) system%external_fields%emfields_xi, system%external_fields%emfields_xf, system%external_fields%emfields_yi, &
                          system%external_fields%emfields_yf, system%external_fields%emfields_zi, system%external_fields%emfields_zf
ELSE
  WRITE(error_unit,*) "CAUTION: please recall that if electric and magnetic fields are zero their domain limits must not be defined in input file"
END IF

! read environment fluid properties
CALL capture_keyword("$environment_fluid_properties")
READ(input_file_unit,*) system%environment_fluid%density, system%environment_fluid%viscosity, system%environment_fluid%velocity
IF (system%environment_fluid%viscosity==0.0) THEN
  WRITE(error_unit,*) "CAUTION: fluid viscosity is zero. Drag forces and Re number will be ignored"
END IF

! read environment fluid thermal properties
IF (system%general_cv%temperature_dofs=="on") THEN
  CALL capture_keyword("$environment_fluid_thermal_properties")
  READ(input_file_unit,*) system%environment_fluid%specific_heat, system%environment_fluid%thermal_conductivity, system%environment_fluid%temperature
  IF (system%environment_fluid%thermal_conductivity==0.0) THEN
    WRITE(error_unit,*) "CAUTION: fluid thermal conductiviy is zero. Convective heat transfer and Pr number will be ignored"
  END IF
END IF

! read nearfield forces switch (on/off)
CALL capture_keyword("$nearfield_forces_switch")
READ(input_file_unit,*) system%general_cv%nearfields_switch

! read nearfields data (coefficients, exponents and cutoff distance)
IF (system%general_cv%nearfields_switch=="on") THEN
  CALL capture_keyword("$number_of_nearfields_sets")
  READ(input_file_unit,*) system%no_nearfields_sets
  IF (system%no_nearfields_sets>0) THEN
    ALLOCATE (system%nearfields_set(1:system%no_nearfields_sets))
    CALL capture_keyword("$nearfields_sets_data")
    DO i=1,system%no_nearfields_sets
      READ(input_file_unit,*) row, system%nearfields_set(i)%kij_attr, system%nearfields_set(i)%expij_attr, &
                              system%nearfields_set(i)%kij_rep, system%nearfields_set(i)%expij_rep, &
                              system%nearfields_set(i)%cutoff_distance
      kij_attr = system%nearfields_set(i)%kij_attr
      expij_attr = system%nearfields_set(i)%expij_attr
      kij_rep = system%nearfields_set(i)%kij_rep
      expij_rep = system%nearfields_set(i)%expij_rep
      dcutoff = system%nearfields_set(i)%cutoff_distance
      IF (kij_attr/=0.0_real8 .AND. (ABS(expij_rep)-ABS(expij_attr))/=0.0_real8) THEN 
        deq = (kij_rep/kij_attr)**(1.0_real8/(ABS(expij_rep)-ABS(expij_attr)))
        IF (dcutoff<deq) WRITE(error_unit,*) "CAUTION: nearfields cutoff distance is less than equilibrium distance in nearfields set", i
      ELSE IF (kij_attr==0.0_real8 .OR. (ABS(expij_rep)-ABS(expij_attr))==0.0_real8) THEN 
        WRITE(error_unit,*) "ERROR: nearfields coefficients lead to division by zero in the computation of deq in nearfields set", i, "; program will be terminated"
        STOP
      END IF
      IF ((kij_attr/=0.0_real8 .AND. expij_attr/=0.0_real8) .AND. (ABS(expij_rep)-ABS(expij_attr))/=0.0_real8) THEN
        dinst = ((kij_rep*ABS(expij_rep))/(kij_attr*ABS(expij_attr)))**(1.0_real8/(ABS(expij_rep)-ABS(expij_attr)))
        IF (dcutoff>dinst) WRITE(error_unit,*) "CAUTION: nearfields cutoff distance exceeds instability distance in nearfields set", i
      ELSE IF ((kij_attr==0.0_real8 .OR. expij_attr==0.0_real8) .OR. (ABS(expij_rep)-ABS(expij_attr))==0.0_real8) THEN
        WRITE(error_unit,*) "ERROR: nearfields coefficients lead to division by zero in the computation of dinst in nearfields set", i, "; program will be terminated"
        STOP
      END IF
    END DO
  ELSE
    WRITE(error_unit,*) "Error: number of nearfields sets is not valid"
    STOP
  END IF
END IF

! read adhesion forces switch (on/off)
CALL capture_keyword("$adhesion_forces_switch")
READ(input_file_unit,*) system%general_cv%adhesion_switch
WRITE(error_unit,*) "CAUTION: make sure the contact model is compatible with the adhesion forces switch"

! read adhesion forces data (model and parameters) (coefficients, exponents and min/max deformations)
IF (system%general_cv%adhesion_switch=="on") THEN
  CALL capture_keyword("$adhesion_model_type_and_number_of_parameters")
  READ(input_file_unit,*) system%adhesion_model%model_type, system%adhesion_model%no_parameters
  IF (system%adhesion_model%no_parameters>0) THEN
    ALLOCATE (system%adhesion_model%model_parameters(1:system%adhesion_model%no_parameters))
    CALL capture_keyword("$adhesion_model_parameters")
    READ(input_file_unit,*) system%adhesion_model%model_parameters
  END IF
END IF

! read global damping coefficient
CALL capture_keyword("$global_damping_coefficient")
READ(input_file_unit,*) system%global_damping_coefficient

! read number of pressurized surfaces (allowed only for certain kinds of spheres)
CALL capture_keyword("$no_pressure_surfaces")
READ(input_file_unit,*) system%no_pressure_surfaces

! read pressurized surface initial values
IF (system%no_pressure_surfaces>0) THEN
  ALLOCATE(system%pressure_surface(1:system%no_pressure_surfaces))
  CALL capture_keyword("$pressure_surface_data")
  DO i=1,system%no_pressure_surfaces
    READ(input_file_unit,*) row, system%pressure_surface(i)%initial_pressure, &
                            system%pressure_surface(i)%initial_influence_area
  END DO                                                                                                         
END IF

! read number of springs
CALL capture_keyword("$no_springs")
READ(input_file_unit,*) system%no_springs

! read spring properties sets
IF (system%no_springs>0) THEN
  ALLOCATE(system%spring(1:system%no_springs))
  CALL capture_keyword("$no_spring_properties_sets")
  READ(input_file_unit,*) system%no_spring_properties_sets
  ALLOCATE(system%spring_properties_set(1:system%no_spring_properties_sets))
  CALL capture_keyword("$spring_properties_sets")
  DO i=1,system%no_spring_properties_sets
    READ(input_file_unit,*) row, system%spring_properties_set(i)%stiffness, system%spring_properties_set(i)%elasticity_modulus, &
                            system%spring_properties_set(i)%critical_stress, system%spring_properties_set(i)%dashpot_constant 
  END DO
END IF

! read spring connectivities and attributes
IF (system%no_springs>0) THEN
  CALL capture_keyword("$spring_connectivities_and_attributes")
  DO i=1,system%no_springs
    READ(input_file_unit,*) row, system%spring(i)%connectivities, system%spring(i)%properties_set_number, &
                            system%spring(i)%initial_length, system%spring(i)%connecting_wall_number, &
                            system%spring(i)%pressure_surface_number
  END DO                                                                                                         
END IF

! read external heating devices´ data
IF (system%general_cv%temperature_dofs=="on") THEN
  CALL capture_keyword("$no_external_heating_devices")
  READ(input_file_unit,*) system%no_external_heating_devices
  IF (system%no_external_heating_devices>0) THEN
    ALLOCATE(system%external_heating_device(1:system%no_external_heating_devices))
    CALL capture_keyword("$external_heating_devices_properties")
    DO i=1,system%no_external_heating_devices
      READ(input_file_unit,*) row, system%external_heating_device(i)%kind, system%external_heating_device(i)%initial_intensity, &
                              system%external_heating_device(i)%initial_intensity_increase_rate, system%external_heating_device(i)%cross_sectional_area, &
                              system%external_heating_device(i)%attenuation_coeff, system%external_heating_device(i)%max_penetration 
    END DO    
    CALL capture_keyword("$external_heating_devices_geometry_and_velocity_data")
    DO i=1,system%no_external_heating_devices
       READ(input_file_unit,*) row, system%external_heating_device(i)%initial_striking_position, system%external_heating_device(i)%axial_direction, &
                               system%external_heating_device(i)%initial_velocity
    END DO    
  END IF
END IF

! periodic boundary conditions
CALL capture_keyword("$periodic_boundary_conditions")
READ(input_file_unit,*) system%general_cv%periodic_bc

! read domain limits for periodic boundary conditions
IF (system%general_cv%periodic_bc=="on") THEN
  CALL capture_keyword("$periodic_bc_domain_boundaries")
  READ(input_file_unit,*) system%general_cv%perbc_xmin, system%general_cv%perbc_xmax, system%general_cv%perbc_ymin, &
                          system%general_cv%perbc_ymax, system%general_cv%perbc_zmin, system%general_cv%perbc_zmax
END IF

! read grid data
CALL capture_keyword("$grid_data")
READ(input_file_unit,*) system%grid%xbeg, system%grid%xend, system%grid%ybeg, system%grid%yend, system%grid%zbeg, system%grid%zend, &
                        system%grid%ndivx, system%grid%ndivy, system%grid%ndivz, system%grid%dt_for_cell_list_update

! read verlet lists data
CALL capture_keyword("$verlet_lists_data")
READ(input_file_unit,*) system%general_cv%verlet_distance, system%general_cv%dt_for_verlet_list_update

! read solution control variables 1
CALL capture_keyword("$solution_control_variables_1")
READ(input_file_unit,*) system%solution_cv%problem_type, system%solution_cv%solver_type, system%solution_cv%contact_detection_algorithm, &
                        system%solution_cv%contact_model, system%solution_cv%rolling_resistance_model

! read solution control variables 2
CALL capture_keyword("$solution_control_variables_2")
READ(input_file_unit,*) system%solution_cv%positions_tolerance, system%solution_cv%velocities_tolerance, &
                        system%solution_cv%desired_no_iterations, system%solution_cv%max_no_iterations

! read number of steps
CALL capture_keyword("$no_steps")
READ(input_file_unit,*) system%solution_cv%no_steps
ALLOCATE(system%solution_cv%step(1:system%solution_cv%no_steps))

! read step control variables
CALL capture_keyword("$step_control_variables")
DO i=1,system%solution_cv%no_steps
  READ(input_file_unit,*) row, system%solution_cv%step(i)%initial_dt, system%solution_cv%step(i)%dt_min, &
                          system%solution_cv%step(i)%dt_max, system%solution_cv%step(i)%final_time, &
                          system%solution_cv%step(i)%adaptive_time_stepping, system%solution_cv%step(i)%collisions_duration_parameter
END DO

! read output control variables
CALL capture_keyword("$output_control_variables")
WRITE(error_unit,*) "CAUTION: please be careful as to specify all output control variables (e.g. do not miss the print_extra_results_file field)"
READ(input_file_unit,*) system%general_cv%results_file_format, system%general_cv%dt_for_results_printing, &
                        system%general_cv%compute_system_properties, system%general_cv%print_rotational_dofs, &
                        system%general_cv%print_extra_results_file

! end reading data
CALL capture_keyword("$end")

END SUBROUTINE read_input_system_data
!!---------------------------------------------------------------------------------------------------------------------


END MODULE input_data_class

!!=============================================================================================================================
!!=============================================================================================================================


!MODULE data_base_input_class

!!USE data_base_types
!USE input_data_class

!IMPLICIT NONE

!CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


!SUBROUTINE input_data_base(db)

!! dummy arguments
!TYPE(data_base), INTENT(OUT) :: db

!! local variables
!INTEGER(KIND=int4) :: i

!CALL capture_keyword("$title")
!READ(input_unit,*) db%title

!CALL capture_keyword("$no_systems")
!READ(input_unit,*) db%no_systems

!! allocate systems
!ALLOCATE(db%system(1:db%no_systems))

!! loop over systems
!DO i=1,db%no_systems
!  CALL read_input_system_data(db%system(i))
!END DO

!CALL capture_keyword("$end")


!END SUBROUTINE input_data_base
!!---------------------------------------------------------------------------------------------------------------------


!END MODULE data_base_input_class

!!=============================================================================================================================
!!=============================================================================================================================


MODULE open_and_close_session

USE input_data_class
USE solver_class

IMPLICIT NONE

CONTAINS
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE open_session

CHARACTER(LEN=8)   :: name
INTEGER(KIND=int4) :: ios

! print program header on screen
CALL header(default_output_unit)

! ask for input file
WRITE(default_output_unit,*)
WRITE(default_output_unit,*)" Enter name of input file (8 characters, no extension)"
READ(default_input_unit,*) name

! set extensions to files
input_file       = TRIM(name)//".inp"
error_file       = TRIM(name)//".err"
psy_results_file = TRIM(name)//".rst"
GiD_results_file = TRIM(name)//".post.res"
GiD_mesh_file    = TRIM(name)//".post.msh"
extra_results_file = TRIM(name)//".xrs"

! open error file (and print header on it)
OPEN(UNIT=error_unit,FILE=error_file,STATUS="UNKNOWN",IOSTAT=ios)
IF(ios/=0)THEN
  WRITE(*,*) "error during opening of error file"
  STOP
END IF
CALL header(error_unit)

! open input file
OPEN(UNIT=input_file_unit,FILE=input_file,STATUS="OLD",IOSTAT=ios)
IF(ios/=0)THEN
  WRITE(error_unit,*) "error during opening of input file"
  STOP
END IF

! open psy results file (and print header on it)
OPEN(UNIT=psy_results_unit,FILE=psy_results_file,STATUS="UNKNOWN",IOSTAT=ios)
IF(ios/=0)THEN
  WRITE(error_unit,*) "error during opening of results file"
  STOP
END IF
CALL header(psy_results_unit)

! open GiD results file
OPEN(UNIT=GiD_results_unit,FILE=GiD_results_file,STATUS="UNKNOWN",IOSTAT=ios)
IF(ios/=0)THEN
  WRITE(error_unit,*) "error during opening of GiD result file"
  STOP
END IF

! open GiD mesh file
OPEN(UNIT=GiD_mesh_unit,FILE=GiD_mesh_file,STATUS="UNKNOWN",IOSTAT=ios)
IF(ios/=0)THEN
  WRITE(error_unit,*) "error during opening of GiD mesh file"
  STOP
END IF

! open extra results file (and print header on it)
OPEN(UNIT=extra_results_unit,FILE=extra_results_file,STATUS="UNKNOWN",IOSTAT=ios)
IF(ios/=0)THEN
  WRITE(error_unit,*) "error during opening of extra results file"
  STOP
END IF
CALL header(extra_results_unit)

END SUBROUTINE open_session
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE close_session

! close input (.inp) file
CLOSE (input_file_unit)

! print footer on psy results file (.rst) and close it
CALL footer(psy_results_unit)
CLOSE (psy_results_unit)

! print footer on extra results file (.xrs) and close it
CALL footer(extra_results_unit)
CLOSE (extra_results_unit)

! close GiD output files (.msh and .post.res)
CLOSE (GiD_mesh_unit)
CLOSE (GiD_results_unit)

! print footer on error file (.err) and close it
CALL footer(error_unit)
CLOSE (error_unit)

! print footer on screen
CALL footer(default_output_unit)

END SUBROUTINE close_session
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE header(u)

INTEGER(KIND=int4), INTENT(IN) :: u

WRITE(u,*) " ##############################################################"
WRITE(u,*) " #                            PSY                             #"
WRITE(u,*) " #             Particle Systems Analysis Program              #"
WRITE(u,*) " ##############################################################"
WRITE(u,*)
WRITE(u,*) "                                               Version 1.7     "
WRITE(u,*) "   (c) 2019            Developed by EDUARDO M. B. CAMPELLO     "
WRITE(u,*)

END SUBROUTINE header
!!---------------------------------------------------------------------------------------------------------------------


SUBROUTINE footer(u)

INTEGER(KIND=int4), INTENT(IN) :: u

WRITE(u,*)
WRITE(u,*) " ##############################################################"
WRITE(u,*) " #                           END                              #"
WRITE(u,*) " ##############################################################" 

END SUBROUTINE footer
!!---------------------------------------------------------------------------------------------------------------------


END MODULE open_and_close_session

!!=============================================================================================================================
!!=============================================================================================================================


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PROGRAM psy_main

USE open_and_close_session
USE basic_types
USE io_files
USE clock
USE deallocate_arrays

IMPLICIT NONE

TYPE(system_data) :: system
INTEGER(int4)     :: t
CHARACTER(3)      :: unit

! open files and units
CALL open_session

! read input data
CALL read_input_system_data(system)

! solve and print results
CALL start_clock(t)
CALL solver(system)
CALL stop_clock(t,unit)

! print computer time for solution
WRITE(default_output_unit,*)
WRITE(default_output_unit,"(A61,I8,A)") "Time Required for this analysis:", t, unit
WRITE(default_output_unit,*) 

! deallocate system arrays
CALL deallocate_system_arrays(system)

! close files and units
CALL close_session

END PROGRAM psy_main

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++