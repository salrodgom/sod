module vector_module
 implicit none
 public
 type  :: vector
  sequence
  real :: x
  real :: y
  real :: z
 end type vector
 contains
!
 pure type(vector) function Array2Vector(a)
  implicit none
  real,intent(in) :: a(1:3)
  array2vector%x = a(1)
  array2vector%y = a(2)
  array2vector%z = a(3)
 end function Array2Vector
!
 pure function Vector2Array(v) result(a)
  implicit none
  type(vector),intent(in) :: v
  real :: a(1:3)
  a(1) = v%x
  a(2) = v%y
  a(3) = v%z
 end function Vector2Array
!
 pure type(vector) function PointCoordinatesAtDistanceFromOriginAndDirection(o,u,d)
  implicit none
  type(vector),intent(in) :: o,u
  real,intent(in)         :: d
  PointCoordinatesAtDistanceFromOriginAndDirection%x = o%x + u%x*d
  PointCoordinatesAtDistanceFromOriginAndDirection%y = o%y + u%y*d
  PointCoordinatesAtDistanceFromOriginAndDirection%z = o%z + u%z*d
 end function PointCoordinatesAtDistanceFromOriginAndDirection
!
 pure type (vector) function vector_unitary(v1)
  implicit none
  type(vector), intent(in) :: v1
  real                     :: norma
  norma=absvec(v1)
  vector_unitary%x = v1%x/norma
  vector_unitary%y = v1%y/norma
  vector_unitary%z = v1%z/norma
 end function vector_unitary
!
 pure type (vector) function vector_scale(v1, r )
  implicit none
  type(vector), intent(in) :: v1
  real,intent(in)          :: r
  vector_scale%x = v1%x*r
  vector_scale%y = v1%y*r
  vector_scale%z = v1%z*r
 end function vector_scale
!
 pure type (vector) function vector_add(v1, v2)
  implicit none
  type (vector), intent(in) :: v1, v2
  vector_add%x = v1%x + v2%x
  vector_add%y = v1%y + v2%y
  vector_add%z = v1%z + v2%z
 end function vector_add
!
 pure type (vector) function vector_sub(v1, v2)
  implicit none
  type (vector), intent(in) :: v1, v2
  vector_sub%x = v1%x - v2%x
  vector_sub%y = v1%y - v2%y
  vector_sub%z = v1%z - v2%z
 end function vector_sub
!
 pure type(vector) function cross(a, b)
  implicit none
  type(vector),intent(in) :: a, b
  cross%x = a%y * b%z - a%z * b%y
  cross%y = a%z * b%x - a%x * b%z
  cross%z = a%x * b%y - a%y * b%x
 end function cross
!
 pure real function dot(a, b)
  implicit none
  type(vector), intent(in) :: a, b
  dot = a%x * b%x + a%y*b%y + a%z*b%z
 end function dot
!
 pure real function absvec(a)
  type (vector), intent (in) :: a
  absvec = sqrt(a%x*a%x + a%y*a%y + a%z*a%z)
 end function absvec
!
 pure real function Dist2Vectors(a, b)
  type(vector), intent(in) :: a,b
  dist2vectors = absvec( vector_sub(b,a) )
 end function  Dist2Vectors
!
 pure type (vector) function bisectriz_vector(v1,v2)
  implicit none
  type (vector), intent(in) :: v1,v2
  type (vector)             :: unitary
  real                      :: r
  bisectriz_vector = vector_add(v2,vector_scale( vector_unitary(vector_sub(v1,v2)),0.5*Dist2Vectors(v1,v2)))
 end function bisectriz_vector
!
 pure real function angle2vectors(a, b)
  ! vectors are the  r_ij vector, (r_ij, r_jk)
  type(vector), intent(in) :: a, b
  angle2vectors = acos( dot(a,b)/(absvec(a)*absvec(b)) )
 end function  angle2vectors
!
 pure real function angle3vectors_jik(r_i, r_j, r_k)
  type(vector), intent(in) :: r_i, r_j, r_k
  type(vector)             :: r_ij, r_ik !, r_jk
  real                     :: x
  r_ij = vector_sub(r_i,r_j)
  !r_jk = vector_sub(r_j,r_k)
  r_ik = vector_sub(r_i,r_k)
  x = dot(r_ij,r_ik)/(absvec(r_ij)*absvec(r_ik))
  if (x > 1.0) then
   angle3vectors_jik = acos(1.0)
  elseif ( x < -1.0 ) then
   angle3vectors_jik = acos(-1.0)
  else
   angle3vectors_jik = acos(x)
  end if
  return
 end function  angle3vectors_jik
!
 pure real function dihedral_angle4vectors_ijkl(r_i, r_j, r_k, r_l)
  type(vector), intent(in) :: r_i, r_j, r_k, r_l
  type(vector)             :: r_ij, r_jk, r_kl
  real                     :: x
  r_ij = vector_sub(r_i,r_j)
  r_jk = vector_sub(r_j,r_k)
  r_kl = vector_sub(r_k,r_l)
  x = dot(cross(r_ij,r_jk),cross(r_jk,r_kl))/(absvec(cross(r_ij,r_jk))*absvec(cross(r_jk,r_kl)))
  if (x > 1.0) then
   dihedral_angle4vectors_ijkl = acos(1.0)
  elseif ( x < -1.0 ) then
   dihedral_angle4vectors_ijkl = acos(-1.0)
  else
   dihedral_angle4vectors_ijkl = acos(x)
  end if
  return
 end function  dihedral_angle4vectors_ijkl
!
 pure real function dihedral_angle4vectors_ijkm(r_i, r_j, r_k, r_m)
  type(vector), intent(in) :: r_i, r_j, r_k, r_m
  type(vector)             :: b_0,b_1,b_2
  real                     :: x
  b_0 = vector_sub(r_j,r_i)
  b_1 = vector_sub(r_k,r_i)
  b_2 = vector_sub(r_m,r_i)
  x = dot(cross(b_0,b_1),cross(b_0,b_2))/(absvec(cross(b_0,b_1))*absvec(cross(b_0,b_2)))
  if (x > 1.0) then
   dihedral_angle4vectors_ijkm = acos(1.0)
  elseif ( x < -1.0 ) then
   dihedral_angle4vectors_ijkm = acos(-1.0)
  else
   dihedral_angle4vectors_ijkm = acos(x)
  end if
  return
 end function  dihedral_angle4vectors_ijkm
end module vector_module
!--------------------------------------------------------------------------------------------------------------------------------
module geometric
 implicit none
 public
 contains
 !
 pure function Crystal2BoxCoordinates(rv, r_c) result (r_b)
  implicit none
  real,intent(in)  :: r_c(1:3), rv(1:3,1:3)
  real             :: r_b(1:3)
  integer          :: j
  forall ( j=1:3 )
   r_b(j) = rv(j,1)*r_c(1)  + rv(j,2)*r_c(2)  + rv(j,3)*r_c(3)
  end forall
 end function Crystal2BoxCoordinates
!
 pure function Box2CrystalCoordinates(vr,r_b) result (r_c)
  implicit none
  real,intent(in)  :: r_b(1:3), vr(1:3,1:3)
  real             :: r_c(1:3)
  integer          :: j
  forall ( j=1:3 )
   r_c(j) = mod(vr(j,1)*r_b(1)  + vr(j,2)*r_b(2)  + vr(j,3)*r_b(3) + 100.0,1.0)
  end forall
 end function Box2CrystalCoordinates
!
 subroutine TwoCoordinatesInSameSpace(rv,r1,r2,v1,v2,s)
  implicit none
  real,intent(in)  :: r1(1:3),r2(1:3),rv(1:3,1:3)!
  real,intent(out) :: v1(1:3),v2(1:3)
  real             :: s
  v2 = r2
  call make_distances(.true.,r1,v2,rv,v1,s)
 end subroutine TwoCoordinatesInSameSpace
!
 subroutine ThreeCoordinatesInSameSpace(rv,r1,r2,r3,v1,v2,v3)
  implicit none
  real,intent(in)  :: r1(1:3),r2(1:3),r3(1:3),rv(1:3,1:3)!
  real,intent(out) :: v1(1:3),v2(1:3),v3(1:3)
  real             :: s
  v2 = r2
  call make_distances(.true.,r1,v2,rv,v1,s)
  call make_distances(.true.,r3,v2,rv,v3,s)
 end subroutine ThreeCoordinatesInSameSpace
!
 subroutine FourCoordinatesInSameSpace(rv,r1,r2,r3,r4,v1,v2,v3,v4)
  implicit none
  real,intent(in)  :: r1(1:3),r2(1:3),r3(1:3),r4(1:3),rv(1:3,1:3)!
  real,intent(out) :: v1(1:3),v2(1:3),v3(1:3),v4(1:3)
  real             :: s
  v1 = r1 ! First "pivote" V1 -> V2 -> V3 -> V4
  call make_distances(.true.,r2,v1,rv,v2,s)
  call make_distances(.true.,r3,v2,rv,v3,s)
  call make_distances(.true.,r4,v3,rv,v4,s)
 end subroutine FourCoordinatesInSameSpace
!
 subroutine Cell(rv,vr,cell_0)
 implicit none
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = acos(-1.0)
 real :: alp, bet
 real :: cosa, cosb, cosg
 real :: gam, sing
 real :: degtorad
 degtorad=pi/180.0
 if(cell_0(4) == 90.0) then
  cosa = 0.0
 else
  alp=cell_0(4)*degtorad
  cosa=cos(alp)
 endif
 if(cell_0(5) == 90.0) then
  cosb = 0.0
 else
  bet = cell_0(5)*degtorad
  cosb = cos(bet)
 endif
 if(cell_0(6) == 90.0) then
  sing = 1.0
  cosg = 0.0
 else
  gam = cell_0(6)*degtorad
  sing = sin(gam)
  cosg = cos(gam)
 endif
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
 return
 end subroutine cell
!
 subroutine uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  real, parameter    :: pi = acos(-1.0)
  real, parameter    :: radtodeg = 180.0/pi
  do i = 1,3
   temp(i) = 0.0
   do j = 1,3
    temp(i) = temp(i) + rv(j,i)*rv(j,i)
   end do
   temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
   temp(3+i) = 0.0
  enddo
  do j = 1,3
   temp(4) = temp(4) + rv(j,2)*rv(j,3)
   temp(5) = temp(5) + rv(j,1)*rv(j,3)
   temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  do i=4,6
   if (abs(cell_0(i) - 90.0 ) < 0.00001) cell_0(i) = 90.0
   if (abs(cell_0(i) - 120.0) < 0.00001) cell_0(i) = 120.0
  end do
  return
 end subroutine uncell
!
 subroutine Inverse(a,c,n)
  implicit none
  integer    :: n
  real       :: a(n,n), c(n,n)
  real       :: l(n,n), u(n,n), b(n), d(n), x(n)
  real       :: coeff
  integer    :: i,j,k
  l=0.0 ; u=0.0 ; b=0.0
  do k=1,n-1
   do i=k+1,n
    coeff=a(i,k)/a(k,k)
    l(i,k) = coeff
    do j=k+1,n
     a(i,j) = a(i,j) - coeff*a(k,j)
    end do
   end do
  end do
  do i=1,n
   l(i,i) = 1.0
  end do
  do j=1,n
   do i=1,j
     u(i,j) = a(i,j)
   end do
  end do
  do k=1,n
   b(k)=1.0
   d(1) = b(1)
   do i=2,n
     d(i)=b(i)
     do j=1,i-1
       d(i) = d(i) - l(i,j)*d(j)
     end do
   end do
   x(n)=d(n)/u(n,n)
   do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
     x(i)=x(i)-u(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
   end do
   do i=1,n
    c(i,k) = x(i)
   end do
   b(k)=0.0
  end do
  return
 end subroutine Inverse
!
subroutine make_distances(flag,r2,r1,rv,r3,dist)
 use vector_module, only: dist2vectors, array2vector
 implicit none
 real,    intent(in)  :: r1(3),r2(3),rv(3,3)    ! coordenadas y matriz de cambio
 real,    intent(out) :: dist,r3(1:3)                     ! matriz de distancias n x n
 real                 :: d_image(1:27),image(3,27)        ! array de distancias
 real                 :: distance,rcm(3),phi
 integer              :: k,l,m,n,o,i,j                    ! variables mudas
 real                 :: atom(3),ouratom(3)               ! coordenadas preparadas
 logical              :: flag                             ! out the coordinate of the atom
! {{ calculamos la matriz de distancias }}
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom = r1
         atom(1) = r2(1) + real(l)
         atom(2) = r2(2) + real(m)
         atom(3) = r2(3) + real(n)
         d_image(k) = dist2vectors( array2vector(Crystal2BoxCoordinates(rv,ouratom)),&
                                    array2vector(Crystal2BoxCoordinates(rv,atom   )) )
         !d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)       ! r3 sera la imagen de atom a la menor distancia
          image(i,k) = atom(i) ! entre todas las imagenes.
         end forall
     enddo
   enddo
  enddo
  dist=minval(d_image)
  if(flag)then
   phi=9999999.999999 ! initial infinite distance
   k=1                ! selected image
   do l=1,27
    if(d_image(l)<=phi)then
      phi=d_image(l) ! seleccionamos la imagen con la menor distancia
      k=l            !
    endif
   enddo
   forall ( l=1:3)
     r3(l)=image(l,k)
   end forall
  else
   r3(1:3)=0.0
  end if
end subroutine make_distances
!
subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
 implicit none
 integer,intent(in) :: n
 real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
 real,intent(out)   :: dist_matrix(n,n)
 integer            :: i,j,k
 real               :: r1(3),r2(3),s,v1(3)
 DO i=1,n
    dist_matrix(i,i)=0.0
    DO j=i+1,n
       forall ( k=1:3 )
        r1(k)=x(k,i)
        r2(k)=x(k,j)
       end forall
       call make_distances(.false.,r1,r2,rv,v1,s)
       !    make_distances(.true.,r1,v2,rv,v1,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
 return
end subroutine make_dist_matrix
!
pure real function volume(rv)
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
!
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  return
end function
!
 pure real function distance(atom,ouratom,rv)
  implicit none
  integer :: j
  real,intent(in) :: atom(3), ouratom(3), rv(3,3)
  real            :: dist(3),o_atom(3),o_ouratom(3)
  forall ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  end forall
  distance = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 end function
end module geometric
!------------------------------------------------------------------------------------------------------------------------------

program cif2lammps
 use iso_fortran_env
 use vector_module
 use geometric
 implicit none
! locals
 integer             :: i,j,k,l,h,hh,m,n,ierr
 integer             :: ii,jj,kk,ll,nn
 real                :: r = 1.0e12
!parameters
 real,parameter      :: k_B = 8.617332478e-5
 real,parameter      :: r_min_criteria_connectivity=0.15,shell_mass=0.1
 INTEGER, PARAMETER  :: muchisimo=100000
! variables
 integer             :: num_args
 integer             :: n_atoms = 0
 real                :: r0=r_min_criteria_connectivity
 real                :: ratom(3),rouratom(3),ratomInSameSpace(3)
 real                :: rv(3,3),vr(3,3),cell_0(1:6) = 0.0
 real                :: rr1(3),rr2(3),vv1(3),vv2(3)
 real                :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
 real                :: xy,xz,yz
 character(len=8)    :: spams
 character(len=100)  :: line
 character(len=20)   :: spam,simulation_type="single"
 character(len=80)   :: string_stop_head= "_atom_site_fract_z" !"_atom_site_charge"
 character(len=100)  :: CIFFilename=" "
 character(len=100)  :: filename=" "
 CHARACTER (LEN=18)  :: string(muchisimo)
 integer,parameter   :: n_atom_types_max=100
 integer             :: n_atom_types = 0
 integer,allocatable :: check_topol_id(:)
 real,allocatable    :: check_topol_distance(:)
 character(len=4)    :: atom_types(n_atom_types_max)
 ! type
 integer,parameter   :: max_n_componets=2
 type  :: node
  integer           :: element
  integer           :: type_
  character(len=2)  :: label_element
  character(len=4)  :: label
  character(len=4)  :: new_label="Xxxx"
  character(len=6)  :: label_from_CIFFile="Xxxxxx"
  character(len=50) :: topol_label="Xxxx"
  character(len=10) :: hybridization="Unknown"
  integer          :: degree
  real             :: charge
  real             :: radius
  real             :: mass
  integer          :: n_components
  character(5)     :: clabel(max_n_componets)
  real             :: xyzs(1:3,max_n_componets)
  real             :: xyzc(1:3,max_n_componets)
 end type
 integer, allocatable, dimension(:)  :: core_shell_body
 integer                             :: n_core_shell_bodies
 ! allocatable
 type(node),allocatable,dimension(:) :: atom
 type(node),allocatable,dimension(:) :: guest_atom
 real,allocatable    :: DistanceMatrix(:,:)
 logical,allocatable :: ConnectedAtoms(:,:)
 integer :: bond_types_max=100, bend_types_max=100, tors_types_max=100, impr_types_max=100
! bonds
 character(len=8),dimension(:), allocatable   :: bond_type_string
 integer,allocatable                          :: bond_type_histogram(:)
 character(len=8)                             :: bond_string
 integer                                      :: n_bond_types = 0
 character(len=80),dimension(:), allocatable  :: bond
 !character(len=8),dimension(:), allocatable   :: bond_type
 integer                                      :: n_bonds = 0
 integer                                      :: n_bonds_max = 100000
! bends
 character(len=12),dimension(:), allocatable  :: bend_type_string
 integer,allocatable                          :: bend_type_histogram(:)
 character(len=12)                            :: bend_string
 integer                                      :: n_bend_types = 0
 character(len=80),dimension(:), allocatable  :: bend
 !character(len=20),dimension(:), allocatable  :: bend_type
 integer                                      :: n_bends = 0
 integer                                      :: n_bends_max = 100000
! torsion (dihedrals)
 character(len=16),dimension(:), allocatable  :: tors_type_string
 integer,allocatable                          :: tors_type_histogram(:)
 character(len=16)                            :: tors_string
 integer                                      :: n_tors_types = 0
 character(len=80),dimension(:), allocatable  :: tors
 !character(len=25),dimension(:), allocatable  :: tors_type
 integer                                      :: n_torss = 0
 integer                                      :: n_torss_max = 100000
! torsion (impropers)
 character(len=16),dimension(:), allocatable  :: impr_type_string
 integer,allocatable                          :: impr_type_histogram(:)
 character(len=16)                            :: impr_string
 integer                                      :: n_impr_types = 0
 character(len=80),dimension(:), allocatable  :: impr
 !character(len=25),dimension(:), allocatable  :: impr_type
 integer                                      :: n_imprs = 0
 integer                                      :: n_imprs_max = 100000
! arguments in line
 character(len=100),dimension(:), allocatable :: args
 logical                                      :: charges_flag         = .true.
 logical                                      :: modify_topology_flag = .false.
 logical                                      :: fill_with_Ar_flag    = .false.
 logical                                      :: adsorption_fast_atom_saturation_INPUT = .false.
 logical                                      :: input_from_RASPA =.false.
 logical                                      :: input_from_MS = .false.
 logical                                      :: flag_naming
 logical                                      :: flag_lammps = .false.
 logical                                      :: ase_flag = .false.
 logical                                      :: add_oxygen_atoms = .false.
 ! Element table and Matrix Topology:
 integer,allocatable  :: TopologyMatrix(:,:)
 integer,parameter    :: max_number_of_elements = 13
 integer              :: element_table(1:max_number_of_elements)
!
 !call random_seed()
!  Si O  H  C  sh He Ne Ar Kr Xe Ge Na Al
!  1  2  3  4  5  6  7  8  9  10 11 12 13
!  14 8  1  6  0  2  10 18 36 54 32 11 13
!Si                     O                      H                      C
 element_table(1) = 14 ; element_table(2) = 8 ; element_table(3) = 1 ; element_table(4) = 6
!sh                     He                     Ne                     Ar
 element_table(5) = 0 ; element_table(6) =2 ; element_table(7) =10   ; element_table(8) =18
!Kr                     Xe                     Ge                     Na
 element_table(9) = 36 ; element_table(10)=54 ; element_table(11)= 32; element_table(12)=11
!Al
 element_table(13) = 13
 !
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-f','--filled','--filled-with-guest','-Ar')
    adsorption_fast_atom_saturation_INPUT=.true.
   case ('-O','--add-oxygen-atoms')
    add_oxygen_atoms = .true.
   case ('-h','--help')
    call print_help()
    stop
   case ('-c','--cif','--in','-i')
    CIFFilename=args(i+1)
    filename=CIFFilename(1:Clen_trim(CIFFilename)-4)
    write(6,'(a)') filename
   case ('-wq','--no-charges')
    charges_flag=.false.
    if(adsorption_fast_atom_saturation_INPUT)then
     string_stop_head="_atom_site_charge"
    else
     string_stop_head= "_atom_site_fract_z"
    end if
   case ('-S','--search-and-modify-topology')
    modify_topology_flag = .true.
   case ('-F','--fill-with-Ar')
    fill_with_Ar_flag = .true.
   case ("-l","--lammps")
    flag_lammps = .true.
   case ('-A','--from-ase')
    ase_flag = .true.
    string_stop_head = "  _atom_site_occupancy"
   case ('-R','--from-RASPA')
    input_from_RASPA = .true.
    string_stop_head = "_atom_site_charge"
   case ('-MS')
    input_from_MS = .true.
    string_stop_head = "_atom_site_occupancy"
  end select
 end do
 open(100,file=CIFfilename,iostat=ierr)
 write(6,'(a)') CIFfilename
 if(ierr/=0) stop 'Error opening CIF file'
 read_cif: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_cif
  if(line(1:14)=="_cell_length_a")then
   read(line,*)spam,cell_0(1)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_b")then
   read(line,*)spam,cell_0(2)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_c")then
   read(line,*)spam,cell_0(3)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_alpha")then
   read(line,*)spam,cell_0(4)
   cycle read_cif
  end if
  if(line(1:16)=="_cell_angle_beta")then
   read(line,*)spam,cell_0(5)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_gamma")then
   read(line,*)spam,cell_0(6)
   cycle read_cif
  end if
  if(line(1:)==string_stop_head) exit read_cif
 end do read_cif
 call cell(rv,vr,cell_0)
 n_atoms=0
 if(add_oxygen_atoms) ii=0
 read_natoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0.or.line(1:5)=="loop_".or.len_trim(line)==0) exit read_natoms
  if (add_oxygen_atoms) then
   n_atoms=n_atoms+5 ! 1 Si, 2 O and 2 shO
  else
   n_atoms=n_atoms+1
  end if
 end do read_natoms
!
 allocate(atom(n_atoms))
 atom(1:n_atoms)%new_label="Xxxx"
 allocate( ConnectedAtoms(n_atoms,n_atoms))
 allocate( DistanceMatrix(n_atoms,n_atoms))
 allocate( TopologyMatrix(n_atoms,max_number_of_elements) )
 TopologyMatrix(1:n_atoms,1:max_number_of_elements) = 0
 allocate( core_shell_body(1:n_atoms))
!
 rewind(100)
 write(6,*)'Atoms:',n_atoms
 do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  if(line(1:)==string_stop_head) exit
 end do
 forall (i=1:n_atom_types_max)
  forall (j=1:4)
   atom_types(i)(j:j)=" "
  end forall
 end forall
 i=0
 read_atoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0)                     exit read_atoms ! the end of the CIF File
  if(i>=1.and.line(1:5)=="loop_") exit read_atoms ! final loop in CIF files from Material Studio
  if(i>=1.and.len_trim(line)==0)  exit read_atoms ! finals lines in CIF files from RASPA
  i=i+1
  atom(i)%n_components=1
  if(input_from_RASPA)then
!loop_
!_atom_site_label
!_atom_site_type_symbol
!_atom_site_fract_x
!_atom_site_fract_y
!_atom_site_fract_z
!_atom_site_charge
   read(line,*,iostat=ierr)atom(i)%label_from_CIFFile, atom(i)%label, (atom(i)%xyzs(j,1),j=1,3), atom(i)%charge
  elseif(ase_flag)then
!loop_
! _atom_site_type_symbol
! _atom_site_label
! _atom_site_symmetry_multiplicity
! _atom_site_fract_x
! _atom_site_fract_y
! _atom_site_fract_z
! _atom_site_occupancy
! Si  Si1       1.0  0.12326  0.05099  0.19816  1.0000
! Si  Si2       1.0  0.37323  0.30098  0.44815  1.0000
! Si  Si3       1.0  0.12316  0.19719  0.44818  1.0000
! Si  Si4       1.0  0.37317  0.44721  0.19816  1.0000
! Si  Si5       1.0  0.12326  0.44721  0.30177  1.0000
! Si  Si6       1.0  0.37327  0.19719  0.05179  1.0000
! Si  Si7       1.0  0.12319  0.30098  0.05179  1.0000
! Si  Si8       1.0  0.37317  0.05099  0.30180  1.0000
! Si  Si9       1.0  0.19646  0.12418  0.05189  1.0000
! Si  Si10      1.0  0.44647  0.37420  0.30190  1.0000
! Si  Si11      1.0  0.05010  0.19730  0.12499  1.0000
   read(line,*,iostat=ierr)atom(i)%label,atom(i)%label_from_CIFFile,xz,(atom(i)%xyzs(j,1),j=1,3),xy
  elseif(input_from_MS)then
!loop_
!_atom_site_label
!_atom_site_type_symbol
!_atom_site_fract_x
!_atom_site_fract_y
!_atom_site_fract_z
!_atom_site_U_iso_or_equiv
!_atom_site_adp_type
!_atom_site_occupancy
   read(line,*,iostat=ierr)atom(i)%label_from_CIFFile, atom(i)%label, (atom(i)%xyzs(j,1),j=1,3), xy, spams, xz
  elseif(charges_flag) then
   read(line,*,iostat=ierr)atom(i)%label,(atom(i)%xyzs(j,1),j=1,3)
  else if(adsorption_fast_atom_saturation_INPUT) then
   read(line,*,iostat=ierr)atom(i)%label_from_CIFFile, atom(i)%label, (atom(i)%xyzs(j,1),j=1,3), atom(i)%charge
  else
   read(line,*,iostat=ierr)atom(i)%label,(atom(i)%xyzs(j,1),j=1,3)
   atom(i)%charge=0.0
  end if
  if(ierr/=0) exit read_atoms
  call CheckAtom(atom(i)%label,atom(i)%mass,atom(i)%radius,atom(i)%element,atom(i)%label_element)
  atom(i)%xyzc(1:3,1) = Crystal2BoxCoordinates(rv, atom(i)%xyzs)  !; atom(i)%xyzc(1:3,1) = rr1(1:3)
  if(add_oxygen_atoms) ii=ii+1
 end do read_atoms
 close(100)
 DistanceMatrix=0.0
 ConnectedAtoms=.false.
 atom(:)%degree= 0
! add oxygen-atoms (and shells
 if (add_oxygen_atoms) then
  write(6,'(a)')'Add Oxygen-atom subroutine:'
  write(6,*)int(n_atoms/5.0), 'Si-atoms detected.'
  loop_Si1: do i=1,int(n_atoms/5.0)
   k=0
   loop_Si2: do j=i+1,int(n_atoms/5.0)
    forall (h=1:3)
     rouratom(h)=atom(i)%xyzs(h,1)
     ratom(h)=atom(j)%xyzs(h,1)
    end forall
    call TwoCoordinatesInSameSpace(rv,ratom,rouratom,vv1,vv2, DistanceMatrix(i,j) )
    DistanceMatrix(j,i) = DistanceMatrix(i,j)
    !if (i/=j .and. DistanceMatrix(i,j) <= 3.8 .and. ConnectedAtoms(i,j).eqv..false. ) then
    if ( i/=j .and. DistanceMatrix(i,j) <= 3.8 ) then
     ConnectedAtoms(j,i) = .true.; ConnectedAtoms(i,j) = .true.
     atom(i)%degree = atom(i)%degree + 1
     atom(j)%degree = atom(j)%degree + 1
     !k = k + 1   ! Si1-degree
     ll = ll + 1 ! Discovered O-atom
     write(6,*)i,j,DistanceMatrix(i,j),atom(i)%degree, atom(j)%degree
     !
!    Core:
     kk=int(n_atoms/5.0) + ll ! local index
     atom(kk)%label="O2  "
     call CheckAtom(atom(kk)%label,atom(kk)%mass,atom(kk)%radius,atom(kk)%element,atom(kk)%label_element)
     atom(kk)%xyzc(1:3,1) = &
      Vector2Array(bisectriz_vector(Array2Vector(Crystal2BoxCoordinates(rv,vv1)),Array2Vector(Crystal2BoxCoordinates(rv,vv2))))
     atom(kk)%xyzs(1:3,1) = Box2CrystalCoordinates(vr, atom(kk)%xyzc(1:3,1) )
     ! 
     !STOP
!    Shell:
     nn=int(n_atoms/5.0) + ll + int(n_atoms/5.0)*2
     atom(nn)%label="shO2"
     call CheckAtom(atom(nn)%label,atom(nn)%mass,atom(nn)%radius,atom(nn)%element,atom(nn)%label_element)
     call random_number(rr1)
     forall (h=1:3)
      atom(nn)%xyzc(h,1) = rr1(h)*0.001 + atom(kk)%xyzc(h,1)
     end forall
     atom(nn)%xyzs(1:3,1) = Box2CrystalCoordinates(vr, atom(nn)%xyzc(1:3,1) )
     !
     if ( atom(i)%degree == 4 ) exit loop_Si2
    end if
   end do loop_Si2
   !STOP
  end do loop_Si1
  ! Detect the closet 2-Si-atoms:
  ! -----------------------------
 else
! topology:
  do i=1,n_atoms
   k=0
   do j=1,n_atoms
    forall (h=1:3)
     rouratom(h)=atom(i)%xyzs(h,1)
     ratom(h)=atom(j)%xyzs(h,1)
    end forall
    call make_distances(.false.,ratom,rouratom,rv,ratomInSameSpace,r)
    DistanceMatrix(i,j)=r
    if(i/=j.and.r<=atom(i)%radius+atom(j)%radius+r_min_criteria_connectivity)then
     k=k+1
     ConnectedAtoms(i,j)=.true.
    end if
   end do
   atom(i)%degree=k
  end do
 end if
! recheck:
! No bonds by constrains:
!    Si O  H  C  sh He Ne Ar Kr Xe Rn Na Al
   ! 1  2  3  4  5  6  7  8  9  10 11 12 13
   ! 14 8  1  6  0  2  10 18 36 54 86 11 13
 do i=1,n_atoms
  do j=1,n_atoms
   if(i/=j.and.ConnectedAtoms(i,j))then
    if(((atom(i)%element==14.and.atom(j)%element==1).or. &
        (atom(i)%element==1.and.atom(j)%element==14)).or.&
        ((atom(i)%element==11.and.atom(j)%element==8).or. &
            (atom(i)%element==8.and.atom(j)%element==11)).or.&
            ((atom(i)%element==11.and.atom(j)%element==0).or. &
                (atom(i)%element==0.and.atom(j)%element==11)).or.&
                ((atom(i)%element==11.and.atom(j)%element==14).or. &
                    (atom(i)%element==14.and.atom(j)%element==11)).or.&
       ((atom(i)%element==14.and.atom(j)%element==6).or. &
        (atom(i)%element==6.and.atom(j)%element==14)).or. &
       ((atom(i)%element==13.and.atom(j)%element==1).or. &
        (atom(i)%element==1.and.atom(j)%element==13)).or.&
       ((atom(i)%element==13.and.atom(j)%element==6).or. &
        (atom(i)%element==6.and.atom(j)%element==13)).or. &
       ((atom(i)%element==13.and.atom(j)%element==14).or. &
        (atom(i)%element==14.and.atom(j)%element==13)).or. &
       ((atom(i)%element==0.and.atom(j)%element==1).or. &
        (atom(i)%element==1.and.atom(j)%element==0)).or. &
       ((atom(i)%element==6.and.atom(j)%element==0).or. &
        (atom(i)%element==0.and.atom(j)%element==6)).or.&
        (atom(i)%element==14.and.atom(j)%element==14).or. &
        (atom(i)%element==8.and.atom(j)%element==8).or. &
        (atom(i)%element==1.and.atom(j)%element==1).or. &
        (atom(i)%element==0.and.atom(j)%element==0).or. &
        (atom(i)%element==11.and.atom(j)%element==11).or. &
        (atom(i)%element==13.and.atom(j)%element==13) &
       )  then
     ConnectedAtoms(i,j)=.false.
     ConnectedAtoms(j,i)=.false.
     atom(i)%degree=atom(i)%degree-1
     atom(j)%degree=atom(j)%degree-1
    end if
   end if
  end do
 end do
! call check_bonds_with_degree(0,3)  ! sh - O
! rename atoms with for certain forcefield
 n_core_shell_bodies = 0
 if ( modify_topology_flag ) then
  !
  write(6,'(a4,1x,a3,1x,13(i2,1x))')'Atom','Z',(element_table(i),i=1,max_number_of_elements)
  absolute_scan: do i=1,n_atoms
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j) )then
     do k=1,max_number_of_elements
      if(  atom(j)%element == element_table(k) )  TopologyMatrix(i,k) = TopologyMatrix(i,k) + 1
     end do
    end if
   end do
   write(atom(i)%topol_label,'(a,a,13(i1,a))')atom(i)%label_element,'@',&
   (TopologyMatrix(i,k),'_', k=1,max_number_of_elements )
   atom(i)%topol_label=adjustl(trim(atom(i)%topol_label) )
   write(6,*) atom(i)%topol_label, (TopologyMatrix(i,k), k=1, max_number_of_elements)
  select case(atom(i)%topol_label)
   ! code:
   ! element@Si_O_H_C_sh_He_Ne_Ar_Kr_Xe_Rn_Na_Al_
   ! Silicon:
   case("Si@0_4_0_0_4_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Si2 "
    atom(i)%charge        = +4.0
    atom(i)%hybridization = "Si_SiO2"
   case("Si@0_4_0_0_3_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Si1 "
    atom(i)%charge        = +4.0
    atom(i)%hybridization = "Si_SiOH"
   case("Si@0_4_0_0_2_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Si1 "
    atom(i)%charge        = +4.0
    atom(i)%hybridization = "Si_SiOH"
   case("Si@0_4_0_0_1_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Si1 "
    atom(i)%charge        = +4.0
    atom(i)%hybridization = "Si_SiOH"
   case("Si@0_4_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Si1 "
    atom(i)%charge        = +4.0
    atom(i)%hybridization = "Si_SiOH"
   ! Aluminium:
   case("Al@0_4_0_0_4_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "Al1 "
    atom(i)%charge        = +3.0
    atom(i)%hybridization = "Al_AlO2"
   ! Cations:
   case("Na@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'Na  '
    atom(i)%charge        = +1.0
    atom(i)%hybridization = "Na"
   ! Oxygens:
   case("O@2_0_0_0_1_0_0_0_0_0_0_0_0_","O@1_0_0_0_1_0_0_0_0_0_0_0_1_")
    atom(i)%new_label     = "O2  "
    atom(i)%charge        = +0.86902
    atom(i)%hybridization = "O_SiO2_core"
    atom(i)%mass          = atom(i)%mass - shell_mass
    n_core_shell_bodies = n_core_shell_bodies + 1
    core_shell_body(i)  = n_core_shell_bodies
   case("O@1_0_1_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "O1  "
    atom(i)%charge        = -1.42600
    atom(i)%hybridization = "O_SiOH"
   case("O@0_0_2_0_1_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "O8  "
    atom(i)%charge        = +1.250000
    atom(i)%hybridization = "O_H2O_core"
    atom(i)%mass          = atom(i)%mass - shell_mass
    n_core_shell_bodies = n_core_shell_bodies + 1
    core_shell_body(i)  = n_core_shell_bodies
   ! Shell:
   case("sh@2_1_0_0_0_0_0_0_0_0_0_0_0_","sh@1_1_0_0_0_0_0_0_0_0_0_0_1_")
    atom(i)%new_label     = "shO2"
    atom(i)%charge        = -2.86902
    atom(i)%hybridization = "O_SiO2_shell"
    atom(i)%mass          = shell_mass
   case("sh@0_1_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "shO8"
    atom(i)%charge        = -2.050000
    atom(i)%hybridization = "O_H2O_shell"
    atom(i)%mass          = shell_mass
   ! Hydrogen:
   case("H@0_1_0_0_0_0_0_0_0_0_0_0_0_","H@0_0_0_1_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "H  "
    atom(i)%charge        = +0.42600
    atom(i)%hybridization = "H_SiOH"
   ! Carbons:  ! UA
   ! ---------------------------------------------------------------------------
   !case("C@0_0_0_0_0_0_0_0_0_0_0_0_0_")
   ! atom(i)%new_label     = "C4  "
   ! atom(i)%charge        =  0.0
   ! atom(i)%hybridization = "C_CH4_sp3"
   !case("C@0_0_0_1_0_0_0_0_0_0_0_0_0_")
   ! atom(i)%new_label     = "C3  "
   ! atom(i)%charge        =  0.0
   ! atom(i)%hybridization = "C_CH3_sp3"
   ! !puede ser:              C_CH2_sp2
   !case("C@0_0_0_2_0_0_0_0_0_0_0_0_0_")
   ! atom(i)%new_label     = "C2  "
   ! atom(i)%charge        = 0.0
   ! atom(i)%hybridization = "C_CH2_sp3"
   ! ! puede ser:             C_CH_sp2
   ! Carbons: All Atoms

   case("C@0_0_4_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C4  "
    atom(i)%charge        = -0.200
    atom(i)%hybridization = "C_CH4_sp3"
   case("C@0_0_3_1_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C3  "
    atom(i)%charge        = +0.006
    atom(i)%hybridization = "C_CH3_sp3"
   case("C@0_0_2_1_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C2  "
    atom(i)%charge        = -0.284
    atom(i)%hybridization = "C_CH2_sp2"
   ! Carbons: United Atoms:
   case("C@0_0_0_1_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = "C2  "
    atom(i)%charge        = 0.0
    atom(i)%hybridization = "C_CH2_sp2"
   case("He@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'He  '
    atom(i)%charge        = +0.0
    atom(i)%hybridization = "He"
   case("Ne@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'Ne  '
    atom(i)%charge        = +0.0
    atom(i)%hybridization = "Ne"
   case("Ar@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'Ar  '
    atom(i)%charge        = +0.0
    atom(i)%hybridization = "Ar"
   case("Kr@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'Kr  '
    atom(i)%charge        = +0.0
    atom(i)%hybridization = "Kr"
   case("Xe@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'Xe  '
    atom(i)%charge        = +0.0
    atom(i)%hybridization = "Xe"
   case("Rn@0_0_0_0_0_0_0_0_0_0_0_0_0_")
    atom(i)%new_label     = 'Rn  '
    atom(i)%charge        = +0.0
    atom(i)%hybridization = "Rn"
   case default
    write(6,*)  "Atom not found"
    write(6,'(a)') atom(i)%topol_label
    do j=1,n_atoms
     if( i/=j .and. DistanceMatrix(i,j)<= 2.5 ) then
      write(6,*) 'Distance:',   DistanceMatrix(i,j)
      write(6,*) 'Connected:',  ConnectedAtoms(i,j)
      write(6,*) 'Name:', atom(j)%element, atom(j)%label_element
      write(6,'(a)')'======================='
     end if
    end do
    STOP
   end select
  end do absolute_scan
! Rename Hydrogens:
  do i=1,n_atoms
   if( atom(i)%element == 1 ) then
    do j=1,n_atoms
     if(i/=j.and.ConnectedAtoms(i,j))then
      if(atom(j)%element == 6) then
      select case( atom(j)%topol_label )
       case("C@0_0_4_0_0_0_0_0_0_0_0_0_0_")
        atom(i)%new_label     = 'H4  '
        atom(i)%charge        = +0.05
        atom(i)%hybridization = "H_CH4_sp3"
       case("C@0_0_3_1_0_0_0_0_0_0_0_0_0_")
        atom(i)%new_label     = 'H3  '
        atom(i)%charge        = -0.002
        atom(i)%hybridization = "H_CH3_sp3"
       case("C@0_0_2_1_0_0_0_0_0_0_0_0_0_")
        atom(i)%new_label     = 'H2  '
        atom(i)%charge        = +0.142
        atom(i)%hybridization = "H_CH2_sp2"
       case default
        write(6,*)  "Hydrogen-atom not found"
        STOP
      end select
      else if ( atom(j)%element == 8 ) then
       select case( atom(j)%topol_label )
       case("O@0_0_2_0_1_0_0_0_0_0_0_0_0_")
        atom(i)%new_label     = 'H8  '
        atom(i)%charge        = +0.400000
        atom(i)%hybridization = "H_H2O"
       case default
        write(6,*)  "Hydrogen-atom not found"
        write(6,*)  atom(j)%topol_label
        STOP
       end select
      end if
     end if
    end do
   end if
  end do
! Reconnect: H8 with shO8:
 do i=1,n_atoms
  if(atom(i)%new_label == "shO8" ) then
   do j=1,n_atoms
    if(i/=j.and.atom(j)%new_label == "H8  ".and.&
     DistanceMatrix(i,j)<=1.2)then
     ConnectedAtoms(i,j)=.true.
     ConnectedAtoms(j,i)=.true.
     atom(i)%degree=atom(i)%degree+1
     atom(j)%degree=atom(j)%degree+1
     write(6,'(a)')'Reconnecting ShO8 - H8 [...]'
    end if
   end do
  end if
 end do
! Remove the links: O2 - Si1,Si2
 do i=1,n_atoms
  if( atom(i)%new_label == "O2  ".or.atom(i)%new_label == "O8  ") then
  do j=1,n_atoms
   if(i/=j.and.ConnectedAtoms(i,j))then
    if(atom(j)%element==14.or.atom(j)%element == 1.or.atom(j)%element==13)  then
     ConnectedAtoms(i,j)=.false.
     ConnectedAtoms(j,i)=.false.
     atom(i)%degree=atom(i)%degree-1
     atom(j)%degree=atom(j)%degree-1
    end if
   end if
  end do
  end if
 end do
 do i=1,n_atoms
  if(atom(i)%new_label=="Xxxx")then
   STOP
  else
   atom(i)%label=atom(i)%new_label
  end if
 end do
  ! }}
 end if
! Core-Shell bodies:
 do i=1,n_atoms
  if(atom(i)%new_label == "O2  " .or. atom(i)%new_label == "O8  " ) then
    do j=1,n_atoms
     if(i/=j.and.(atom(j)%new_label == "shO2".or.atom(j)%new_label == "shO8")&
        .and.ConnectedAtoms(i,j))then
      core_shell_body(j) = core_shell_body(i)
     end if 
    end do
  end if
 end do
! atom types:
 do i = 1,n_atoms
  check_atom_type: if(i==1)then
   n_atom_types=1
   atom_types(n_atom_types)=atom(i)%label
   atom(i)%type_=n_atom_types
  else
   n_atom_types=n_atom_types+1
   do h=1,n_atom_types-1
    if(atom(i)%label==atom_types(h))then
     n_atom_types=n_atom_types-1
     atom(i)%type_=h
     exit check_atom_type
    end if
   end do
   atom_types(n_atom_types)=atom(i)%label
   atom(i)%type_=n_atom_types
  end if check_atom_type
  write(6,'(a4,1x,3(f10.8,1x),a2,1x,3(f10.6,1x),1x,a20)') &
   atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),&
   atom(i)%label_element,atom(i)%charge,atom(i)%radius,atom(i)%mass, atom(i)%topol_label
 end do
! degree of each node:
 write(6,*)'=========='
 write(6,*)'atom_types:',n_atom_types
 do i=1,n_atom_types
  write(6,*)atom_types(i)
 end do
 open(987,file="atom_types_for_dump.txt")
 write(987,'(100(a4,1x))') ( atom_types(i), i=1,n_atom_types)
 close(987)
! lammps:
 if ( flag_lammps ) then
! bends:
 allocate(bend_type_string(bend_types_max))
 allocate(bend_type_histogram(bend_types_max))
 allocate(bend(n_bends_max))
 bend_type_histogram=0
 bend_string(1:12)="            "
 do i=1,bend_types_max
  do j=1,12
   write(bend_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_bends_max)
  forall (j=1:80)
   bend(i)(j:j)=" "
  end forall
 end forall
 do_i_bend: do i=1,n_atoms ! central atom
  do_j_bend: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_bend: do k=1,n_atoms
     if(ConnectedAtoms(i,k).and.k/=j.and.k/=i)then
      write(bend_string(1:12),'(3a4)') atom(j)%label(1:4),atom(i)%label(1:4),atom(k)%label(1:4)
      if(n_bend_types==0)then
       n_bends=n_bends+1
       n_bend_types=n_bend_types+1
       bend_type_string(1)=bend_string
       bend_type_histogram(1)=1
       write(bend(n_bends)(1:28),'(4i7)')n_bend_types,j,i,k
       write(bend(n_bends)(29:31),'(a3)')' # '
       write(bend(n_bends)(32:),'(a)') bend_string
       cycle do_k_bend
      else
       n_bends=n_bends+1
       n_bend_types=n_bend_types+1 ! try
       check_bend_type: do h=n_bend_types-1,1,-1
        if(bend_string(1:12)==bend_type_string(h)(1:12).or.&
           bend_string(1:12)==bend_type_string(h)(9:12)//bend_type_string(h)(5:8)//&
           bend_type_string(h)(1:4))then
         n_bend_types=n_bend_types-1
         write(bend(n_bends)(1:28),'(4i7)')h,j,i,k
         write(bend(n_bends)(29:31),'(a3)')' # '
         write(bend(n_bends)(32:),'(a)') bend_string
         check_bend: do l=1,n_bends-1
          read(bend(l)(8:28),'(3i7)')jj,ii,kk
          if((jj==j.and.ii==i.and.kk==k).or.&
             (kk==j.and.ii==i.and.jj==k))then
           n_bends=n_bends-1
           cycle do_k_bend
          end if
         end do check_bend
         bend_type_histogram(h)=bend_type_histogram(h)+1
         cycle do_k_bend
        end if
       end do check_bend_type
       bend_type_histogram(n_bend_types)=1
       bend_type_string(n_bend_types)=bend_string
       write(bend(n_bends)(1:28),'(4i7)')n_bend_types,j,i,k
       write(bend(n_bends)(29:31),'(a3)')' # '
       write(bend(n_bends)(32:),'(a)') bend_string
      end if
     end if
    end do do_k_bend
   end if
  end do do_j_bend
 end do do_i_bend
! ...
! Remove the links: O2 - Si1,Si2
 ! Si1 - shO2
 ! Si2 - shO2
 ! Si1 - O1
 do i=1,n_atoms
  select case( atom(i)%new_label )
  case("Si1 ","Si2  ","Al1 ")
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j).and.(atom(j)%element==0.or.atom(j)%element==8))then
     ConnectedAtoms(i,j)=.false.
     ConnectedAtoms(j,i)=.false.
     atom(i)%degree=atom(i)%degree-1
     atom(j)%degree=atom(j)%degree-1
    end if
   end do
  case("shO2","O1  ")
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j).and.(atom(j)%element==14.or.atom(j)%element==13))then
     ConnectedAtoms(i,j)=.false.
     ConnectedAtoms(j,i)=.false.
     atom(i)%degree=atom(i)%degree-1
     atom(j)%degree=atom(j)%degree-1
    end if
   end do
  case("O8")
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j).and.atom(j)%element==1)then
      ConnectedAtoms(i,j)=.false.
      ConnectedAtoms(j,i)=.false.
      atom(i)%degree=atom(i)%degree-1
      atom(j)%degree=atom(j)%degree-1
    end if
   end do
  end select
 end do
 write(6,*)'=========='
 write(6,'(a)') 'Connectivity array for each atom:'
 write(6,'(1000(i1))')(atom(i)%degree,i=1,n_atoms)
 write(6,*)'=========='
 allocate(bond_type_string(bond_types_max))
 allocate(bond_type_histogram(bond_types_max))
 allocate(bond(n_bonds_max))
 bond_type_histogram=0
 bond_string(1:8)="        "
 do i=1,bond_types_max
  do j=1,8
   write(bond_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_bonds_max)
  forall (j=1:80)
   bond(i)(j:j)=" "
  end forall
 end forall
 do_i: do i=1,n_atoms
  do_j: do j=i+1,n_atoms
   if(ConnectedAtoms(i,j))then
    write(bond_string(1:8),'(2a4)') atom(i)%label(1:4),atom(j)%label(1:4)
    forall (l=1:100)
     line(l:l)=" "
    end forall
    if(n_bond_types==0)then
     n_bonds=n_bonds+1
     n_bond_types=n_bond_types+1
     bond_type_string(1)=bond_string
     bond_type_histogram(1)=1
     write(bond(n_bonds)(1:21),'(3i7)')n_bond_types,i,j
     write(bond(n_bonds)(22:24),'(a3)')' # '
     write(bond(n_bonds)(25:),'(a)') bond_string
     cycle do_j
    else
     n_bonds=n_bonds+1
     n_bond_types=n_bond_types+1 ! try
     check_bond: do h=n_bond_types-1,1,-1
      if(bond_string(1:8)==bond_type_string(h)(1:8).or.&
         bond_string(1:8)==bond_type_string(h)(5:8)//bond_type_string(h)(1:4))then
       n_bond_types=n_bond_types-1
       bond_type_histogram(h)=bond_type_histogram(h)+1
       write(bond(n_bonds)(1:21),'(3i7)')h,i,j
       write(bond(n_bonds)(22:24),'(a3)')' # '
       write(bond(n_bonds)(25:),'(a)') bond_string
       cycle do_j
      end if
     end do check_bond
     bond_type_histogram(n_bond_types)=1
     bond_type_string(n_bond_types)=bond_string
     write(bond(n_bonds)(1:21),'(3i7)')n_bond_types,i,j
     write(bond(n_bonds)(22:24),'(a3)')' # '
     write(bond(n_bonds)(25:),'(a)') bond_string
    end if
   end if
  end do do_j
 end do do_i
! Torsions:
! dihedrals
 allocate(tors_type_string(tors_types_max))
 allocate(tors_type_histogram(tors_types_max))
 allocate(tors(n_torss_max))
 tors_type_histogram=0
 tors_string(1:16)="                "
 do i=1,tors_types_max
  do j=1,16
   write(tors_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_torss_max)
  forall (j=1:80)
   tors(i)(j:j)=" "
  end forall
 end forall
 do_i_tors: do i=1,n_atoms
  do_j_tors: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_tors: do k=1,n_atoms
     if(ConnectedAtoms(j,k).and.k/=j.and.k/=i)then
      do_l_tors: do l=1,n_atoms
       if(ConnectedAtoms(k,l).and.l/=k.and.l/=j.and.l/=i)then
        write(tors_string(1:16),'(4a4)')atom(i)%label(1:4),atom(j)%label(1:4),&
                                        atom(k)%label(1:4),atom(l)%label(1:4)
        if(n_tors_types==0)then
         n_torss=n_torss+1
         n_tors_types=n_tors_types+1
         tors_type_string(1)=tors_string
         tors_type_histogram(1)=1
         write(tors(n_torss)(1:35),'(5i7)')n_tors_types,i,j,k,l
         write(tors(n_torss)(36:38),'(a3)')' # '
         write(tors(n_torss)(39:),'(a)') tors_string
         cycle do_l_tors
        else
         n_torss=n_torss+1
         n_tors_types=n_tors_types+1 ! try
         check_tors_type: do h=n_tors_types-1,1,-1
          if(tors_string(1:16)==tors_type_string(h)(1:16).or.&
             tors_string(1:16)==tors_type_string(h)(13:16)//tors_type_string(h)(9:12)//&
                                tors_type_string(h)(5:8)//tors_type_string(h)(1:4))then
           n_tors_types=n_tors_types-1
           write(tors(n_torss)(1:35),'(5i7)')h,i,j,k,l
           write(tors(n_torss)(36:38),'(a3)')' # '
           write(tors(n_torss)(39:),'(a)') tors_string
           check_tors: do m=1,n_torss-1
            read(tors(m)(8:35),'(4i7)')ii,jj,kk,ll
            if((ii==i.and.jj==j.and.kk==k.and.ll==l).or.&
               (ll==i.and.kk==j.and.jj==k.and.ii==l))then
             n_torss=n_torss-1
             cycle do_l_tors
            end if
           end do check_tors
           tors_type_histogram(h)=tors_type_histogram(h)+1
           cycle do_l_tors
          end if
         end do check_tors_type
         tors_type_histogram(n_tors_types)=1
         tors_type_string(n_tors_types)=tors_string
         write(tors(n_torss)(1:35),'(5i7)') n_tors_types,i,j,k,l
         write(tors(n_torss)(36:38),'(a3)')' # '
         write(tors(n_torss)(39:),'(a)') tors_string
        end if
       end if
      end do do_l_tors
     end if
    end do do_k_tors
   end if
  end do do_j_tors
end do do_i_tors
! impropers
! =========
 allocate(impr_type_string(impr_types_max))
 allocate(impr_type_histogram(impr_types_max))
 allocate(impr(n_imprs_max))
 impr_type_histogram=0
 impr_string(1:16)="                "
 do i=1,impr_types_max
  do j=1,16
   write(impr_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_imprs_max)
  forall (j=1:80)
   impr(i)(j:j)=" "
  end forall
 end forall
!   l
!  .|.
!   i - j
!   |
!   k
 do_i_impr: do i=1,n_atoms
  do_j_impr: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_impr: do k=1,n_atoms
     if(ConnectedAtoms(i,k).and.k/=j.and.k/=i)then
      do_l_impr: do l=1,n_atoms
       if(ConnectedAtoms(i,l).and.l/=k.and.l/=j.and.l/=i)then
        write(impr_string(1:16),'(4a4)')atom(i)%label(1:4),atom(j)%label(1:4),&
                                        atom(k)%label(1:4),atom(l)%label(1:4)
        if(n_impr_types==0)then
         n_imprs=n_imprs+1
         n_impr_types=n_impr_types+1
         impr_type_string(1)=impr_string
         impr_type_histogram(1)=1
         write(impr(n_imprs)(1:35),'(5i7)')n_impr_types,i,j,k,l
         write(impr(n_imprs)(36:38),'(a3)')' # '
         write(impr(n_imprs)(39:),'(a)') impr_string
         cycle do_l_impr
        else
         n_imprs=n_imprs+1
         n_impr_types=n_impr_types+1 ! try
         check_impr_type: do h=n_impr_types-1,1,-1
          if(impr_string(1:16)==impr_type_string(h)(1:16).or.&
             impr_string(1:16)==impr_type_string(h)(1:4)//impr_type_string(h)(9:12)//&
                                impr_type_string(h)(5:8)//impr_type_string(h)(13:16).or.&
             impr_string(1:16)==impr_type_string(h)(1:4)//impr_type_string(h)(5:8)//&
                                impr_type_string(h)(13:16)//impr_type_string(h)(9:12).or.&
             impr_string(1:16)==impr_type_string(h)(1:4)//impr_type_string(h)(13:16)//&
                                impr_type_string(h)(9:12)//impr_type_string(h)(5:8) )then
           n_impr_types=n_impr_types-1
           write(impr(n_imprs)(1:35),'(5i7)')h,i,j,k,l
           write(impr(n_imprs)(36:38),'(a3)')' # '
           write(impr(n_imprs)(39:),'(a)') impr_string
           check_impr: do m=1,n_imprs-1
            read(impr(m)(8:35),'(4i7)')ii,jj,kk,ll
            if((ii==i.and.jj==j.and.kk==k.and.ll==l).or.&
               (ii==i.and.kk==j.and.jj==k.and.ll==l).or.&
               (ii==i.and.jj==j.and.ll==k.and.kk==l).or.&
               (ii==i.and.ll==j.and.jj==k.and.kk==l).or.&
               (ii==i.and.kk==j.and.ll==k.and.jj==l).or.&
               (ii==i.and.ll==j.and.kk==k.and.jj==l))then
             n_imprs=n_imprs-1
             cycle do_l_impr
            end if
           end do check_impr
           impr_type_histogram(h)=impr_type_histogram(h)+1
           cycle do_l_impr
          end if
         end do check_impr_type
         impr_type_histogram(n_impr_types)=1
         impr_type_string(n_impr_types)=impr_string
         write(impr(n_imprs)(1:35),'(5i7)') n_impr_types,i,j,k,l
         write(impr(n_imprs)(36:38),'(a3)')' # '
         write(impr(n_imprs)(39:),'(a)') impr_string
        end if
       end if
      end do do_l_impr
     end if
    end do do_k_impr
   end if
  end do do_j_impr
 end do do_i_impr
! ...
 write(6,*)'Bond types:',n_bond_types
 write(6,*)'bonds:',n_bonds
 do i=1,n_bond_types
  write(6,*) bond_type_string(i),bond_type_histogram(i)
 end do
 write(6,*)'Bend types:',n_bend_types
 write(6,*)'bends:',n_bends
 do i=1,n_bend_types
  write(6,*) bend_type_string(i),bend_type_histogram(i)
 end do
 write(6,*)'Dihedral types:',n_tors_types
 write(6,*)'dihedrals:',n_torss
 do i=1,n_tors_types
  write(6,*) tors_type_string(i),tors_type_histogram(i)
 end do
 write(6,*)'Improper types:',n_impr_types
 write(6,*)'impropers:',n_imprs
 do i=1,n_impr_types
  write(6,*) impr_type_string(i),impr_type_histogram(i)
 end do
! output:
 call cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
      xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
 call output_lammps()
 deallocate(bend_type_string)
 deallocate(bend_type_histogram)
 deallocate(bend)
 deallocate(bond_type_string)
 deallocate(bond_type_histogram)
 end if
 call output_gulp()
 !call output_pdb()
 call output_CIF()
 deallocate(TopologyMatrix)
 deallocate(DistanceMatrix)
 deallocate(ConnectedAtoms)
 deallocate(atom)
 stop 'Done'
 contains
!
 subroutine remove_bond_in_atom( III )
  implicit none
  integer,intent(in) :: III
  integer            :: i,j,k
  write(6,*)'=============================='
  write(6,*)'Strange ',atom(III)%element,'-atom:',III,atom(III)%degree, atom(III)%label
  allocate( check_topol_id(1:atom(III)%degree))
  allocate( check_topol_distance(1:atom(III)%degree) )
  j=0
  do i=1,n_atoms
   if(III/=i.and.ConnectedAtoms(III,i))then
    j=j+1
    check_topol_id(j)=i
    check_topol_distance(j)=DistanceMatrix(III,i)
    write(6,*)atom(i)%label,DistanceMatrix(III,i),j,atom(i)%degree
   end if
  end do
  ghj: do i=1,atom(III)%degree
   j=check_topol_id(i)
   if( DistanceMatrix(III,j) == maxval( check_topol_distance ) ) then
    ConnectedAtoms(III,j)=.false.
    ConnectedAtoms(j,III)=.false.
    atom(III)%degree=atom(III)%degree-1
    atom(j)%degree=atom(j)%degree-1
    write(6,*)'Solution:'
    write(6,*)'remove bond:',III,j,atom(III)%label,atom(j)%label, DistanceMatrix(III,j), maxval( check_topol_distance )
    exit ghj
   end if
  end do ghj
  deallocate( check_topol_id )
  deallocate( check_topol_distance )
  return
 end subroutine remove_bond_in_atom
!
 subroutine check_bonds_with_degree(ZZZ,KKK)
 implicit none
 integer,intent(in) :: ZZZ,KKK
 do i=1,n_atoms
  do while ( atom(i)%element == ZZZ .and. atom(i)%degree > KKK )
   write(6,*)'=============================='
   write(6,*)'Strange ',ZZZ,'-atom:',i,atom(i)%degree
   h=0
   allocate( check_topol_id(1:atom(i)%degree) )
   allocate( check_topol_distance(1:atom(i)%degree) )
   do j=1,n_atoms
    if(i/=j.and.ConnectedAtoms(i,j))then
     h=h+1
     check_topol_id(h)=j
     check_topol_distance(h)=DistanceMatrix(i,j)
     write(6,*)atom(j)%label,DistanceMatrix(i,j),h,r,atom(j)%degree
    end if
   end do
   do h=1,atom(i)%degree
    j=check_topol_id(h)
    if( DistanceMatrix(i,j) == maxval( check_topol_distance ) ) then
     m=j
     ConnectedAtoms(i,j)=.false.
     ConnectedAtoms(j,i)=.false.
     atom(i)%degree=atom(i)%degree-1
     atom(j)%degree=atom(j)%degree-1
    end if
   end do
   write(6,*)'Solution:'
   write(6,*)'remove bond:', i,m,atom(i)%label,atom(m)%label
   deallocate( check_topol_id )
   deallocate( check_topol_distance )
  end do
 end do
 return
 end subroutine check_bonds_with_degree
 subroutine search_forcefield(string,interaction,label)
  implicit none
  character(len=80),intent(out):: string
  real                         :: p1,p2,p3,x,y
  integer                      :: p0,p4
  character(len=10)            :: type_fff
  character(len=4),intent(in)  :: interaction
  character(len=4),intent(in)  :: label(1:4)
  character(len=4)             :: ourlabel(1:4)
  character(len=80)            :: units
  character(len=80)            :: passing
  character(len=80)            :: variables
  integer                      :: u=123,ierr=0
  forall (ii=1:80)
   string(ii:ii)=" "
   units(ii:ii)=" "
   passing(ii:ii)=" "
   variables(ii:ii)=" "
  end forall
  open(u,file='forcefield.lib')
  fff: select case (interaction)
   case('bond')
    initbond: do
     read(u,'(a)') line
     if(line(1:11)=='Bond Coeffs') then
      read(line(12:),'(a)') units
      units=adjustl(trim(units))
      exit initbond
     end if
    end do initbond
    do
     read(u,'(a)') line
     if(line(1:12)=='Angle Coeffs'.or.line(1:3)=="End") then
      write(string,*) 0.0,0.0
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(2(a4,1x),a)')ourlabel(1),ourlabel(2),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("harmonic")
         read(variables,*) p1,p2
         select case(units)
          case('kcal/mol','kcal/mol/A/A')
           p1=kcalmol2ev(p1) ! kcal/mol/A/A
           p2=p2             ! A
          case('K','K/A/A','RASPA')
           p1=0.5*K2eV(p1)  ! lammps:     K*(r-r_o)**2
           p2=p2            ! raspa:  0.5*K*(r-r_o)**2
         end select
         write(string,*) p1,p2,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown bond potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        p1=0.0
        p2=1.0
        write(string,*) p1,p2,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('bend')
    initbend: do
     read(u,'(a)') line
     if(line(1:12)=='Angle Coeffs') then
      read(line(13:),'(a)') units
      units=adjustl(trim(units))
      !write(6,*) units
      exit initbend
     end if
    end do initbend
    do
     read(u,'(a)') line
     if(line(1:15)=='Dihedral Coeffs'.or.line(1:3)=="End") then
      write(string,*) 0.0,0.0
      exit fff
     end if
     if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(3(a4,1x),a)')ourlabel(1),ourlabel(2),ourlabel(3),passing
      if((adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("harmonic")
         read(variables,*) p1,p2
         select case(units)
          case('kcal/mol','kcal/mol/rad/rad')
           p1=kcalmol2ev(p1) ! kcal/mol/A/A
           p2=p2             ! A
          case('K','K/rad/rad','RASPA')
           p1=0.5*K2eV(p1)
           p2=p2
         end select
         write(string,*) p1,p2,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown bend potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        p1=0.0
        p2=90.0
        write(string,*) p1,p2,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('tors')
    inittors: do
     read(u,'(a)')line
     if(line(1:15)=='Dihedral Coeffs') then
      read(line(16:),'(a)') units
      units=adjustl(trim(units))
      exit inittors
     end if
    end do inittors
    do
     read(u,'(a)') line
     if(line(1:15)=='Improper Coeffs'.or.line(1:3)=="End") then
      type_fff="none"
      write(string,*)1,0.0,1,0.0
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(4(a4,1x),a)')(ourlabel(ii),ii=1,4),passing
      if((adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("fourier")
         read(variables,*) p0,p1,p4,p2
         select case(units)
          case('kcal/mol')
           p0=p0
           p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
           p4=p4
           p2=p2             ! deg
          case('K','K/rad/rad','RASPA')
           p0=p0
           p1=K2eV(p1)
           p2=p2
           p4=p4
         end select
         write(string,*) p0,p1,p4,p2,' # ',type_fff
        case default
         write(6,'(a,1x,a)')'Unknown dihedral potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        write(string,*)1,0.0,1,0.0
       end if
       exit fff
      end if
     end if
    end do
   case('impr')
    initimpr: do
     read(u,'(a)') line
     if(line(1:15)=='Improper Coeffs') then
      read(line(16:),'(a)') units
      units=adjustl(trim(units))
      exit initimpr
     end if
    end do initimpr
    do
     read(u,'(a)',iostat=ierr) line
     if (ierr/=0) exit fff
     if(line(1:11)=='Pair Coeffs'.or.line(1:3)=="End") then
      type_fff="none"
      write(string,*) 0.0,1,1
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(4(a4,1x),a)')(ourlabel(ii),ii=1,4),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i j k l
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&  ! 1 2 3 4
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i j l k
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&  ! 1 2 4 3
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i k j l
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&  ! 1 3 2 4
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i k l j
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&  ! 1 3 4 2
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i l j k
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(4))).and.&  ! 1 4 2 3
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&  ! i l k j
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(4))).and.&  ! 1 4 3 2
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(2)))) )then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("cossq")        ! 0.5*p1*cos(x-p2)^2
         read(variables,*) p1,p2
         select case(units)
          case('kcal/mol')
           p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
           p2=p2             ! [-]
         end select
         write(string,*) p1,p2,' # ',type_fff
        case("cvff")         ! p1*(1+p0*cos(p4*x))
         read(variables,*) p1,p0,p4
         select case(units)
          case('kcal/mol')
           p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
           p0=p0             ! [-]
           p4=p4             ! [-]
          case('K','K/rad/rad','RASPA')
           p1=K2eV(p1)
           p0=p0
           p4=p4
         end select
         write(string,*) p1,p0,p4
        case default
         write(6,'(a,1x,a)')'Unknown improper potential:', adjustl(trim(type_fff))
         stop
        end select
       else
        write(string,*) 0.0,1,1
       end if
       exit fff
      end if
     end if
    end do
   case('pair')
    initpair: do
     read(u,'(a)') line
     if(line(1:11)=='Pair Coeffs') then
      read(line(12:),'(a)') units
      units=adjustl(trim(units))
      exit initpair
     end if
    end do initpair
    do
     read(u,'(a)',iostat=ierr) line
     if(ierr/=0) exit fff
     if(line(1:3)=="End") then
      type_fff = "none"
      !write(string,*)'lj/cut/coul/long',0.0,1.0
      write(string,*)'coul/long'
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(2(a4,1x),a)')ourlabel(1),ourlabel(2),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        select case(adjustl(trim(type_fff)))
        case("buck/long","buck/cut","buck","buck/cs")
         read(variables,*) p1,p2,p3
         select case(units)
          case('RASPA')
           p1=K2eV(p1)   ! A
           p2=1.0/p2     ! rho
           p3=K2eV(p3)   ! C
         end select
         write(string,*)'buck/coul/long/cs',p1,p2,p3
        case("lj/long","lj/cut","lennard")
         read(variables,*) p1,p2
         select case(units)
          case('RASPA')
           p1=K2eV(p1)  ! epsilon
           p2=p2        ! sigma
         end select
         write(string,*)'lj/cut/coul/long',p1,p2
        case("6_12_pot")
         read(variables,*) p1,p2
         select case(units)
          case('RASPA')
           p1=K2eV(p1)  ! A
           p2=K2eV(p2)  ! B
           x = p2*p2/(4.0*p1)
           y = (p1/p2)**(1.0/6.0)
           p1 = x
           p2 = y
         end select
         write(string,*)'lj/cut/coul/long',p1,p2
        end select
       else
        type_fff="lj/cut"
        write(string,*)'lj/cut/coul/long',0.0,1.0
       end if
       exit fff
      end if
     end if
    end do
  end select fff
  close(u)
  return
 end subroutine search_forcefield
 subroutine output_CIF()
 implicit none
 character(len=120) :: CIFFilenameNew
 integer           :: u=1000
 integer           :: i
 character(len=12)  :: extension="_topol.cif"
 CIFFilenameNew=filename(1:Clen_trim(filename))//extension
 CIFFilenameNew=adjustl(CIFFilenameNew)
 open(u,file=CIFFilenameNew)
 write(u,'(a)')'data_subtitutions'
 write(u,'(a)')'_audit_creation_method  iGOR'
 write(u,'(a)')"_audit_author_name 'Sponge Bob'"
 write(u,'(a,5x,f14.7)')'_cell_length_a',cell_0(1)
 write(u,'(a,5x,f14.7)')'_cell_length_b',cell_0(2)
 write(u,'(a,5x,f14.7)')'_cell_length_c',cell_0(3)
 write(u,'(a,5x,f14.7)')'_cell_angle_alpha',cell_0(4)
 write(u,'(a,5x,f14.7)')'_cell_angle_beta',cell_0(5)
 write(u,'(a,5x,f14.7)')'_cell_angle_gamma',cell_0(6)
 write(u,'(a,5x,f14.7)')'_cell_volume',volume(rv)
 write(u,'(a)')"_symmetry_space_group_name_hall 'p 1'"
 write(u,'(a)')"_symmetry_space_group_name_h-m 'p 1'"
 write(u,'(a)')'_symmetry_int_tables_number 1'
 write(u,'(a)')"_symmetry_equiv_pos_as_xyz 'x,y,z'"
 write(u,'(a)')'loop_'
 write(u,'(a)')'_atom_site_label'
 write(u,'(a)')'_atom_site_fract_x'
 write(u,'(a)')'_atom_site_fract_y'
 write(u,'(a)')'_atom_site_fract_z'
 write(u,'(a)')'_atom_site_charge'
  do i=1,n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7,1x)')&
   atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),atom(i)%charge
  end do
 close(u)
 end subroutine output_CIF
 subroutine output_pdb()
  implicit none
  integer           :: u=333,i
  character(len=4)  :: extension=".pdb"
  character(len=100) :: PDBFilename
  PDBFilename=filename(1:Clen_trim(filename))//extension
  PDBFilename=adjustl(PDBfilename)
  open(u,file=PDBFilename)
  write(u,'(a5,i5)')'MODEL',k
  write(u,'(a6,9(f14.7,1x))')'REMARK',&
   xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,xz,yz
  write(u,'(a6,1x,a,1x,i6,1x,a,i6)')'REMARK','# atoms',n_atoms,'# timestep',int(0.0)
  write(u,'(a6,1x,a,1x,f14.7)')'REMARK','Volume:',volume(rv)
  write(u,'(a6,3f9.3,3f7.2)')'CRYST1',(cell_0(j),j=1,6)
  do i=1,n_atoms
   write(u,'(a6,i5,1x,a2,3x,a4,1x,i4,4x,3f8.3,2f6.2,10x,a2)') &
   'ATOM  ',i,atom(i)%label,'MOL ',0,(atom(i)%xyzc(j,1),j=1,3),0.0,0.0,atom(i)%label
  end do
  write(u,'(a6)')'ENDMDL'
  close(u)
 end subroutine output_pdb
!
 subroutine output_gulp()
  implicit none
  character(len=100) :: GULPFilename
  integer           :: u=444
  integer           :: i
  real              :: mmm,rrr
  integer           :: zzz
  character(len=2)  :: zlz
  character(len=6)  :: atomtype_library
  character(len=4)  :: extension=".gin"
  character(len=10) :: newlabel
  GULPFilename=filename(1:Clen_trim(filename))//extension
  GULPFilename=adjustl(GULPfilename)
  !adjustl(
  open(u,file=GULPFilename)
  write(u,'(a)')'opti conp prop'
  write(u,'(A)')'cell'
  write(u,'(6(f9.5,1x))') (cell_0(j) , j=1,6)
  write(u,'(A)')'fractional'
  do i=1,n_atoms
   select case(atom(i)%label)
    case("Si  ",'Si0 ':'Si99')
     newlabel="Si  core  "
    case("O   ",'O0  ':'O999')
     newlabel="O2  core  "
    case("shO ",'shO0':'shO9')
     newlabel = "O2  shel  "
    case default
     write(6,'(a)') 'Atom nof found in the library'
     write(6,'(a)') atom(i)%label
   end select
   write(u,'(a10,1x,3(f14.7,1x),1x,f14.7)')newlabel,(atom(i)%xyzs(j,1),j=1,3),atom(i)%charge
  end do
   write(u,'(a,1x,i3)')'species',n_atom_types
   write(u,'(a)')'Si core  Si    core'
   write(u,'(a)')'O2 core  O_O2- core'
   write(u,'(a)')'O2 shel  O_O2- shel'
   write(u,'(a)')'end'
  !do i=1,n_atom_types
  !call checkatomtype_gulp(atom_types(i), atomtype_library)
  !select case(atom_types(i))
  ! case("shO2")
  !  write(u,'(a4,1x,a4,1x,a10)')'O2 
  ! case("H0  ":"H999")
  !  write(u,'(a4,1x,a4,1x,a6)') "H_  ", "core", atomtype_library
  ! case default
  !  write(u,'(a4,1x,a4,1x,a6)') atom_types(i), "core", atomtype_library
  !end select
  !end do
  write(u,'(a)') "library catlow"
  write(u,'(a)') "switch_minimiser rfo gnorm 0.05"
  write(u,'(a)') "stepmx opt 0.1"
  write(u,'(a)') 'dump every 1 out.grs'
  write(u,'(a)') 'output cif out.cif'
  close(u)
 end subroutine output_gulp
 !
 subroutine output_lammps()
  implicit none
  integer            :: day,hour,i4_huge,milli,minute,month,second,year
  integer            :: ii,jj
  character(len=10)  :: time
  character(len=8)   :: date
  integer            :: u=222
  real               :: mmm,rrr
  integer            :: zzz
  integer,parameter  :: seed = 86456
  character(len=2)   :: zlz
  character(len=120) :: DataFilename
  character(len=5)   :: extension=".data"
  character(len=4)   :: label(4)
  character(len=80)  :: str
  call srand(seed)
  DataFilename=filename(1:Clen_trim(filename))//extension
  DataFilename=adjustl(trim(DataFilename))
  forall (ii=1:4)
   forall (jj=1:4)
    label(ii)(jj:jj)=" "
   end forall
  end forall
  open(u,file=DataFilename)
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  write(u,'(a,1x,i2,1x,i2,1x,i4,1x,i2,a1,i2,a1,i2)')'Created on',day,month,year,hour,':',minute,':',second
  write(u,*)' '
  write(u,*)n_atoms,' atoms'
  if (n_bonds>0) write(u,*)n_bonds,' bonds'
  if (n_bends>0) write(u,*)n_bends,' angles'
  if (n_torss>0) write(u,*)n_torss,' dihedrals'
  if (n_imprs>0) write(u,*)n_imprs,' impropers'
  write(u,*)' '
  write(u,*)n_atom_types,' atom types'
  if (n_bonds>0) write(u,*)n_bond_types,' bond types'
  if (n_bends>0) write(u,*)n_bend_types,' angle types'
  if (n_torss>0) write(u,*)n_tors_types,' dihedral types'
  if (n_imprs>0) write(u,*)n_impr_types,' improper types'
  write(u,*)' '
  write(u,*)xlo_bound,xhi_bound,' xlo xhi'
  write(u,*)ylo_bound,yhi_bound,' ylo yhi'
  write(u,*)zlo_bound,zhi_bound,' zlo zhi'
  write(u,*)xy,xz,yz,' xy xz yz'
  write(u,*)' '
  write(u,'(a)')'Masses'
  write(u,*)' '
  do i=1,n_atom_types
   call CheckAtom(atom_types(i),mmm,rrr,zzz,zlz)
   ! O_core -> O_atom - O_shell
   if(atom_types(i)=="O2  ") mmm=mmm-shell_mass
   write(u,'(i4,1x,f14.7,1x,a,a)')i,mmm,' # ',atom_types(i)
  end do
  if( n_bonds>0) then
  write(u,*)' '
  write(u,'(a)')'Bond Coeffs'
  write(u,*)' '
  do i=1,n_bond_types
   read(bond_type_string(i),'(2a4)')(label(ii),ii=1,2)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'bond',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',bond_type_string(i)
  end do
  end if
  if( n_bends>0 ) then
  write(u,*)' '
  write(u,'(a)')'Angle Coeffs'
  write(u,*)' '
  do i=1,n_bend_types
   read(bend_type_string(i),'(3a4)')(label(ii),ii=1,3)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'bend',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',bend_type_string(i)
  end do
  end if
  if ( n_torss>0 ) then
  write(u,*)' '
  write(u,'(a)')'Dihedral Coeffs'
  write(u,*)' '
  do i=1,n_tors_types
   read(tors_type_string(i),'(4a4)')(label(ii),ii=1,4)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'tors',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',tors_type_string(i)
  end do
  end if
  if( n_imprs>0 ) then
  write(u,*)' '
  write(u,'(a)')'Improper Coeffs'
  write(u,*)' '
  do i=1,n_impr_types
   read(impr_type_string(i),'(4a4)')(label(ii),ii=1,4)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'impr',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',impr_type_string(i)
  end do
  end if
  write(u,*)' '
  write(u,'(a)')'PairIJ Coeffs'
  write(u,*)' '
  do i=1,n_atom_types
   do j=i,n_atom_types
    label(1)=atom_types(i)
    label(2)=atom_types(j)
    forall (ii=1:80)
     str(ii:ii)=" "
    end forall
    call search_forcefield(str,'pair',label)
    write(u,*)i,j,str,' # ', (label(ii),ii=1,2)
   end do
  end do
  write(u,*)' '
  write(u,'(a)')'Atoms'
  write(u,*)' '
  do i=1,n_atoms
   select case(atom(i)%hybridization)
   case("O_SiO2_shell")
    write(u,'(i5,1x,i3,1x,i3,1x,f10.7,1x,3(f14.5,3x),a1,1x,a,1x,a)')&
   i,1,atom(i)%type_,atom(i)%charge,(atom(i)%xyzc(j,1)-rand()*0.001,j=1,3),'#',atom(i)%hybridization,atom(i)%topol_label
   case default
    write(u,'(i5,1x,i3,1x,i3,1x,f10.7,1x,3(f14.7,1x),a1,1x,a,1x,a)')&
   i,1,atom(i)%type_,atom(i)%charge,(atom(i)%xyzc(j,1),j=1,3),'#',atom(i)%hybridization,atom(i)%topol_label
   end select 
  end do
  write(u,*)' '
  if(n_bonds>0)then
   write(u,'(a)')'Bonds'
   write(u,*)' '
   do i=1,n_bonds
    write(u,'(i8,a)')i,bond(i)
   end do
   write(u,*)' '
  end if
  if(n_bends>0)then
   write(u,'(a)')'Angles'
   write(u,*)' '
   do i=1,n_bends
    write(u,'(i8,a)')i,bend(i)
   end do
   write(u,*)' '
  end if
  if( n_torss> 0)then
   write(u,'(a)')'Dihedrals'
   write(u,*)' '
   do i=1,n_torss
    write(u,'(i8,a)')i,tors(i)
   end do
   write(u,*)' '
  end if
  if(n_imprs>0)then
   write(u,'(a)')'Impropers'
   write(u,*)' '
   do i=1,n_imprs
    write(u,'(i8,a)')i,impr(i)
   end do
   write(u,*)' '
  end if
  !if(n_core_shell_bodies > 0) then
  !  write(u,'(a)')'CS-Info'
  !  write(u,*)' '
  !  ! Core-Shell bodies:
  !  do i=1,n_atoms
  !   if(atom(i)%new_label == "O2  ".or. atom(i)%new_label == "O8  ".or.&
  !      atom(i)%new_label == "shO2".or. atom(i)%new_label == "shO8") then
  !    write(u,'(i8,1x,i8)')i, core_shell_body(i)
  !   else
  !    write(u,'(i8,1x,a)')i, 'NULL'
  !   end if
  !  end do 
  !end if
  !
  close(u)
 end subroutine output_lammps
 subroutine cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
                              xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
  implicit none
  real,intent(in)  :: cell_0(6)
  real,intent(out) :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
  real,intent(out) :: xy,xz,yz
  real,parameter :: pi=acos(-1.0)
  real,parameter :: radtodeg = 180.0/pi
  real,parameter :: degtorad = pi/180.0
  real :: xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz,cosa,cosb,cosg
  real :: alp,bet,gam,sing
  IF(cell_0(4) == 90.0) THEN
    cosa = 0.0
  ELSE
    ALP=cell_0(4)*degtorad
    COSA=cos(ALP)
  ENDIF
  IF(cell_0(5) == 90.0) THEN
    cosb = 0.0
  ELSE
    bet = cell_0(5)*degtorad
    cosb = cos(bet)
  ENDIF
  IF(cell_0(6) == 90.0) then
    sing = 1.0
    cosg = 0.0
  ELSE
    gam = cell_0(6)*degtorad
    sing = sin(gam)
    cosg = cos(gam)
  ENDIF
  lx=cell_0(1)
  xy=cell_0(2)*cosg
  xz=cell_0(3)*cosb
  ly=sqrt(cell_0(2)*cell_0(2)-xy*xy)
  yz=(cell_0(2)*cell_0(3)*cosa-xy*xz)/ly
  lz=sqrt(cell_0(3)*cell_0(3)-xz*xz-yz*yz)
  xlo_bound=0.0
  ylo_bound=0.0
  zlo_bound=0.0
  xhi_bound=lx
  yhi_bound=ly
  zhi_bound=lz
  !xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
  !xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
  !ylo_bound = ylo + MIN(0.0,yz)
  !yhi_bound = yhi + MAX(0.0,yz)
  !zlo_bound = zlo
  !zhi_bound = zhi
  return
 end subroutine cellnormal2lammps

 subroutine checkatomtype_gulp(label,newlabel)
  implicit none
  character(len=4),intent(in)  :: label
  character(len=10),intent(out) :: newlabel
  select case(label)
   case("H   ","H1  ":"H999")
    newlabel = "H    "
   case("O1  ")
    newlabel = "O1   "
   case("O2  ")
    newlabel = "O_O2- core"
   case("shO2")
    newlabel = "O_O2- shel"
   case("Si  ","Si1 ":"Si99")
    newlabel = "Si    core"
   case("Al  ","Al1 ":"Al99")
    newlabel = "Al   "
   case("Na  ","Na1 ":"Na99")
    newlabel = "Na   "
   case("C1  ","C4  ")
    newlabel = "C_R2  "
   case("C2  ","C6  ")
    newlabel = "C_R   "
   case("C5  ")
    newlabel = "C_R3  "
   case("C3  ","C8  ")
    newlabel = "C_3   "
   case("C7  ")
    newlabel = "C_32  "
   case default
    write(6,'(a,1x,a)')"Atom unknowed:",label
    newlabel = label
    STOP
  end select
  return
 end subroutine checkatomtype_gulp

 subroutine checkatom(Label,m,s,Z,Zlabel)
  implicit none
  character(len=4),intent(in)  :: Label
  real,intent(out)             :: m,s
  integer,intent(out)          :: Z
  character(len=2),intent(out) :: ZLabel
  select case(Label)
   ! 1  2 3 4 5 6 7  8  9  10 11 12
   ! 14 8 1 6 0 2 10 18 36 54 86 11
   case('Na  ','Na0 ':'Na99')
    Z=11
    m=22.98976928
    s=0.0
    Zlabel='Na'
   case('C   ','C0  ':'C999')
    Z=6
    m=12.0107
    s=1.0
    ZLabel=' C'
   case('H   ','H0  ':'H999')
    Z=1
    m=1.00794
    s=0.620
    Zlabel=' H'
   case('O   ','O0  ':'O999','OH  ')
    Z=8
    m=15.9994
    s=0.70
    Zlabel=' O'
   case('Si  ','Si0 ':'Si99')
    Z=14
    m=28.0855
    s=1.14
    zlabel='Si'
   case('Al  ','Al1 ':'Al99')
    Z=13
    m=26.9800
    s=1.14
    zlabel='Al'
   case('He  ')
    Z=2
    m=4.0026
    s=0.0001
    Zlabel='He'
   case('Ne  ')
    Z=10
    m=20.179
    s=0.0001
    Zlabel='Ne'
   case('Ar  ')
    Z=18
    m=39.948
    s=0.0001
    Zlabel='Xe'
   case('Kr  ')
    Z=36
    m=83.80
    s=0.0001
    Zlabel='Xe'
   case('Xe  ')
    Z=54
    m=131.293
    s=0.0001
    Zlabel='Xe'
   case('Ge  ')
    Z=32
    m=72.64
    s=1.14
    Zlabel='Xe'
   case('shO ','shO1':'shO9','ShO ')
    Z=0
    m=shell_mass
    s=1.2
    Zlabel='sh'
   case default
    write(6,'(a1,a4,a1)')"'",label,"'"
    STOP 'Atom unknowed in checkatom subroutine'
  end select
 end subroutine checkatom

 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
! SUBROUTINE cell(rv,vr,cell_0)
! implicit none
! integer :: i,j
! real, intent(in)  :: cell_0(6)
! real, intent(out) :: rv(3,3),vr(3,3)
! real, parameter   :: pi = ACOS(-1.0)
! real :: alp,bet
! real :: cosa,cosb,cosg
! real :: gam,sing
! real :: DEGTORAD
! DEGTORAD=pi/180.0
! IF(cell_0(4) == 90.0) THEN
!   cosa = 0.0
! ELSE
!   ALP=cell_0(4)*degtorad
!   COSA=cos(ALP)
! ENDIF
! IF(cell_0(5) == 90.0) THEN
!   cosb = 0.0
! ELSE
!   bet = cell_0(5)*degtorad
!   cosb = cos(bet)
! ENDIF
! IF(cell_0(6) == 90.0) then
!   sing = 1.0
!   cosg = 0.0
! ELSE
!   gam = cell_0(6)*degtorad
!   sing = sin(gam)
!   cosg = cos(gam)
! ENDIF
! rv(1,1) = cell_0(1)
! rv(1,2) = cell_0(2)*cosg
! rv(1,3) = cell_0(3)*cosb
! rv(2,1) = 0.0
! rv(2,2) = cell_0(2)*sing
! rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
! rv(3,1) = 0.0
! rv(3,2) = 0.0
! rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
! call inverse(rv,vr,3)
! WRITE(6,'(a)') 'Cell:'
! WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
! WRITE(6,'(a)')'Linear Transformation Operator:'
! DO i=1,3
!  WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
! ENDDO
! WRITE(6,'(a)')'----------------------------------------'
! WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
! DO i=1,3
!  WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
! ENDDO
! WRITE(6,'(a)')'----------------------------------------'
! RETURN
! END SUBROUTINE cell
!!
! SUBROUTINE uncell(rv,cell_0)
!  implicit none
!  real,intent(out)   :: cell_0(6)
!  real,intent(in)    :: rv(3,3)
!  integer            :: i,j
!  real               :: temp(6)
!  REAL               :: radtodeg
!  REAL, PARAMETER    :: pi=ACOS(-1.0)
!  radtodeg=180.0/PI
!  do i = 1,3
!    temp(i) = 0.0
!    do j = 1,3
!      temp(i) = temp(i) + rv(j,i)**2
!    enddo
!    temp(i) = sqrt(temp(i))
!  enddo
!  cell_0(1) = abs(temp(1))
!  cell_0(2) = abs(temp(2))
!  cell_0(3) = abs(temp(3))
!  do i = 1,3
!    temp(3+i) = 0.0
!  enddo
!  do j = 1,3
!    temp(4) = temp(4) + rv(j,2)*rv(j,3)
!    temp(5) = temp(5) + rv(j,1)*rv(j,3)
!    temp(6) = temp(6) + rv(j,1)*rv(j,2)
!  enddo
!  temp(4) = temp(4)/(temp(2)*temp(3))
!  temp(5) = temp(5)/(temp(1)*temp(3))
!  temp(6) = temp(6)/(temp(1)*temp(2))
!  cell_0(4) = radtodeg*acos(temp(4))
!  cell_0(5) = radtodeg*acos(temp(5))
!  cell_0(6) = radtodeg*acos(temp(6))
!  DO i=4,6
!     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
!     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
!  ENDDO
!  RETURN
! END SUBROUTINE uncell
!!
! SUBROUTINE inverse(a,c,n)
! implicit none
! integer n
! real a(n,n), c(n,n)
! real L(n,n), U(n,n), b(n), d(n), x(n)
! real coeff
! integer i, j, k
! L=0.0
! U=0.0
! b=0.0
! do k=1, n-1
!   do i=k+1,n
!      coeff=a(i,k)/a(k,k)
!      L(i,k) = coeff
!      do j=k+1,n
!         a(i,j) = a(i,j)-coeff*a(k,j)
!      end do
!   end do
! end do
! do i=1,n
!  L(i,i) = 1.0
! end do
! do j=1,n
!  do i=1,j
!    U(i,j) = a(i,j)
!  end do
! end do
! do k=1,n
!  b(k)=1.0
!  d(1) = b(1)
!  do i=2,n
!    d(i)=b(i)
!    do j=1,i-1
!      d(i) = d(i) - L(i,j)*d(j)
!    end do
!  end do
!  x(n)=d(n)/U(n,n)
!  do i = n-1,1,-1
!    x(i) = d(i)
!    do j=n,i+1,-1
!      x(i)=x(i)-U(i,j)*x(j)
!    end do
!    x(i) = x(i)/u(i,i)
!  end do
!  do i=1,n
!    c(i,k) = x(i)
!  end do
!  b(k)=0.0
! end do
! RETURN
! END SUBROUTINE inverse
!
!subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
! implicit none
! integer,intent(in) :: n
! real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
! real,intent(out)   :: dist_matrix(n,n)
! integer            :: i,j,k
! real               :: r1(3),r2(3),s
! DO i=1,n
!    dist_matrix(i,i)=0.0
!    DO j=i+1,n
!       forall ( k=1:3 )
!        r1(k)=x(k,i)
!        r2(k)=x(k,j)
!       end forall
!       call make_distances(cell_0,r1,r2,rv,s)
!       dist_matrix(i,j)=s
!       dist_matrix(j,i)=dist_matrix(i,j)
!    END DO
! END DO
! return
!end subroutine make_dist_matrix
!
! SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
! IMPLICIT NONE
! REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
! REAL,    intent(out) :: dist
! REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
! INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
! REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
!  k=0
!  do l=-1,1
!   do m=-1,1
!      do n=-1,1
!         k = k + 1
!         ouratom(1) = r1(1)
!         ouratom(2) = r1(2)
!         ouratom(3) = r1(3)
!         atom(1) = r2(1) + l
!         atom(2) = r2(2) + m
!         atom(3) = r2(3) + n
!         d_image(k) = distance(atom,ouratom,rv)
!         forall ( i=1:3)
!           image(i,k) = atom(i)
!         end forall
!     enddo
!   enddo
!  enddo
!  dist=MINVAL(d_image)
!  RETURN
! END SUBROUTINE
!!
!real function volume(rv)
!  implicit none
!  real, intent(in)  :: rv(3,3)
!  real       :: r1x
!  real       :: r1y
!  real       :: r1z
!  real       :: r2x
!  real       :: r2y
!  real       :: r2z
!  real       :: r3x
!  real       :: r3y
!  real       :: r3z
!  real       :: vol
!!
!  r1x = rv(1,1)
!  r1y = rv(2,1)
!  r1z = rv(3,1)
!  r2x = rv(1,2)
!  r2y = rv(2,2)
!  r2z = rv(3,2)
!  r3x = rv(1,3)
!  r3y = rv(2,3)
!  r3z = rv(3,3)
!  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
!  volume = abs(vol)
!  RETURN
!end function
!
 !REAL FUNCTION DISTANCE(atom,ouratom,rv)
 ! IMPLICIT NONE
 ! INTEGER :: j
 ! REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
 ! REAL :: rv(3,3)
 ! FORALL ( j=1:3 )
 !  o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
 !  o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
 !  dist(j) = o_ouratom(j) - o_atom(j)
 ! END FORALL
 ! DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 !END FUNCTION
!
 real function K2eV(x)
  implicit none
  real    :: x
  real,parameter :: kB_1=11604.52211052
  K2eV=x/kB_1
  return
 end function
!
 real function kcalmol2ev(x)
  implicit none
  real            :: x
  kcalmol2ev = 0.04336*x
  return
 end function  kcalmol2ev
 real function kjmol2ev(x)
  implicit none
  real,intent(in) :: x
  kjmol2ev=x*(1036.427/100000.0)
  return
 end  function kjmol2ev
 real function kjmolnmnm2evAA(x)
  real, intent(in) :: x
  kjmolnmnm2evaa = kjmol2ev(x)/100.0
  return
 end function  kjmolnmnm2evaa
 subroutine print_help()
    print '(a)', '  -h,   --help,   print usage information and exit'
    print '(a)', '  -c,   --cif,    CIF File input'
    print '(a)', '  -S,             Scan mode'
    print '(a)', '  -wq,            Do not read charges from CIF File'
 end subroutine print_help
end program cif2lammps
