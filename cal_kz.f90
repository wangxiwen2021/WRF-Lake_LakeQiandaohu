!!!!!!!!!!!!!!!!!!!!
!compile: gfortran -o cal_kz cal_kz.f90 ${NETCDF_LIB}/libnetcdff.a ${NETCDF_LIB}/libnetcdf.a -I${NETCDF_INC}
!!!!!!!!!!!!!!!!!!!!

program main
  use netcdf
  implicit none
  real, parameter :: denh2o = 1.000e3
  real, parameter :: cpliq = 4.188e3
  !real, parameter :: dtime = 24*3600. 
  real, parameter :: dtime = 30*24*3600. 
  real, parameter :: dz = 0.5
  real, parameter :: beta = 0.4                 ! fraction solar rad absorbed at surface
  real, parameter :: za = 0.6                   !  base of surface absorption layer (m)
  integer, parameter :: defval = -999           ! missing value
  integer, parameter :: NO = 7

  real, allocatable :: tlak(:,:)  ! [time, depth]
  real, allocatable :: dTdt(:,:), dTdz(:,:), sumdTdt(:,:), dphidumdz(:,:)       ! [time, depth]
  real, allocatable :: phidum(:,:)      ! [time, depth]
  real, allocatable :: ke(:,:)      ! [time, depth]
  real, allocatable :: rhow(:,:)        ! [time, depth]
  real, allocatable :: zlak(:)          ! [depth]
  real, allocatable :: sabg(:)          ! [time]
  integer, allocatable :: junk(:,:)     ! [time, 4]
  real :: zin, zout, rsfin, rsfout
  real :: eta(12)       !monthly

  integer :: ncid, varid
  integer :: Timeid, Timelen
  integer :: Depthsid, Depthslen
  integer :: i, k       ! i for time, k for depth
  integer :: start(2), counts(2)
  character(len=*), parameter :: Timename = 'Time_d' 
  character(len=*), parameter :: Depthsname = 'Depth' 
  !character(len=*), parameter :: filename = './Data/QDH-Daba-2016.nc'
  character(len=*), parameter :: filename = './Data/QDH-Daba-2016-mon.nc'
  character(len=*), parameter :: varname = 't_lake_utc04'

  integer :: ncid_o, Timeid_o, Depthsid_o, varid_o1, varid_o2, varid_o3, varid_o4, varid_o5, varid_o6
  !character(len=*), parameter :: filename_o = './Data/kme-2016.nc'
  character(len=*), parameter :: filename_o = './Data/kme-2016-mon-with-phidum.nc'
  character(len=*), parameter :: varname_o = 'kme'

!----read TLAK from netcdf-----
  call check( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
  call check( nf90_inq_dimid(ncid, Timename, Timeid) )
  call check( nf90_inquire_dimension(ncid, Timeid, len=Timelen) )
  call check( nf90_inq_dimid(ncid, Depthsname, Depthsid) )
  call check( nf90_inquire_dimension(ncid, Depthsid, len=Depthslen) )

  print*, Timelen
  print*, Depthslen
  allocate( tlak(Depthslen, Timelen) )  !A variable with shape [time, level,latitude,longitude] should be declared as var(longitude, latitude, level,time) in Fortran.
  allocate( dTdt(Depthslen, Timelen) )
  allocate( dTdz(Depthslen, Timelen) )
  allocate( sumdTdt(Depthslen, Timelen) )
  allocate( phidum(Depthslen, Timelen) )
  allocate( dphidumdz(Depthslen, Timelen) )
  allocate( rhow(Depthslen, Timelen) )
  allocate( ke(Depthslen, Timelen) )
  allocate( zlak(Depthslen) )
  allocate( sabg(Timelen) )
  allocate( junk(Timelen, 4))
  counts = (/Depthslen, Timelen/)
  start = (/1, 1/)

  call check( NF90_INQ_VARID(ncid, trim(varname), varid) )
  call check( NF90_GET_VAR(ncid, varid, tlak, start=start, count=counts) )
  call check( NF90_CLOSE(ncid) )
  !print*, tlak(4, 5)  !(4, 3) in ncl
  !print*, tlak(4, 7)  !(6, 3) in ncl
  !print*, tlak(6, 9)  !(8, 5) in ncl

!----read SABG from txt-----
  !open(unit=13, file="./Data/era5-utc04.txt", status="old")
  open(unit=13, file="./Data/era5-utc04-mon.txt", status="old")
  do i=1, size(sabg)
    !read(13,'(I6, 3I7, f7.2)') junk(i,1), junk(i,2), junk(i,3), junk(i,4), sabg(i)      !year, month, day, hour 
    read(13, '(f7.2)') sabg(i)
  end do

!----initialize other values----
  eta = (/0.266986311, 0.260304327, 0.28458672, 0.43999031, 0.601756239, 0.443918799, 0.451999672 &
        , 0.421387264, 0.374253468, 0.324082232, 0.243322242, 0.268365931/)
  zlak = (/(k, k=5, 640, 5)/)*0.1
  rhow = 1000.*( 1.0 - 1.9549e-05*(abs(tlak-277.)) )**1.68
  dTdt = defval; dTdz = defval; sumdTdt = defval 

!----start calculation----
  do k = 1, Depthslen
    dTdt(k,:) = der2_forward(y=tlak(k,:), dx=dtime)
  end do
  do i = 1, Timelen
    dTdz(:,i) = der2_center(y=tlak(:,i), dx=dz)
  end do
  do i = 2, Timelen-1
    do k = 1, Depthslen-1
      sumdTdt(k, i) = integrate(x=zlak(k:Depthslen), y=dTdt(k:Depthslen, i))
    end do
  end do

  do k=1, Depthslen 
    zin = zlak(k) - 0.5*dz
    zout = zlak(k) + 0.5*dz
    do i=1, Timelen
     rsfin = exp( -eta(i)*max( zin-za, 0.) )
     rsfout = exp( -eta(i)*max( zout-za, 0.) )
     phidum(k, i) = (rsfin-rsfout) * sabg(i) * (1.-beta)   
    end do
  end do

  do i = 1, Timelen
    dphidumdz(:,i) = der2_center(y=phidum(:,i), dx=dz)          
  end do

  print*, 'dTdt: ', dTdt(:,NO) 
  print*, ''
  print*, 'sumdTdt: ', sumdTdt(:,NO) 
  print*, ''
  print*, 'tlak: ', tlak(:,NO)
  print*, ''
  print*, 'dTdz: ', dTdz(:,NO)
  print*, ''

  do i=1, Timelen
    do k=1, Depthslen
      if (dTdz(k,i).eq.0. .or. dTdz(k,i).eq.defval .or. sumdTdt(k,i).eq.defval .or. phidum(k,i).eq.defval) then
        ke(k,i) = defval
      else
        ke(k,i) = ( sumdTdt(k,i) - phidum(k,i)/(rhow(k,i)*cpliq) ) / (-dTdz(k,i))
        !ke(k,i) = sumdTdt(k,i) / (-dTdz(k,i))
      end if
    end do
  end do
  print*, 'ke: ', ke(:,NO)
  print*, ''
  
  call check( nf90_create(trim(filename_o), NF90_CLOBBER, ncid_o) )
  call check( nf90_def_dim(ncid_o, Timename, Timelen, Timeid_o))
  call check( nf90_def_dim(ncid_o, Depthsname, Depthslen, Depthsid_o))
  call check( nf90_def_var(ncid_o, "kme", NF90_FLOAT, (/Depthsid_o, Timeid_o/), varid_o1))
  call check( nf90_def_var(ncid_o, "sumdTdt", NF90_FLOAT, (/Depthsid_o, Timeid_o/), varid_o2))
  call check( nf90_def_var(ncid_o, "dTdz", NF90_FLOAT, (/Depthsid_o, Timeid_o/), varid_o3))
  call check( nf90_def_var(ncid_o, "phidum", NF90_FLOAT, (/Depthsid_o, Timeid_o/), varid_o4))
  call check( nf90_def_var(ncid_o, "dphidumdz", NF90_FLOAT, (/Depthsid_o, Timeid_o/), varid_o5))
  call check( nf90_def_var(ncid_o, "dTdt", NF90_FLOAT, (/Depthsid_o, Timeid_o/), varid_o6))
  call check( nf90_enddef(ncid_o) )
  call check( nf90_put_var(ncid_o, varid_o1, ke))
  call check( nf90_put_var(ncid_o, varid_o2, sumdTdt))
  call check( nf90_put_var(ncid_o, varid_o3, dTdz))
  call check( nf90_put_var(ncid_o, varid_o4, phidum))
  call check( nf90_put_var(ncid_o, varid_o5, dphidumdz))
  call check( nf90_put_var(ncid_o, varid_o6, dTdt))
  call check( nf90_close(ncid_o) )

  contains
 
  subroutine check(errorcode)
      implicit none
      integer, intent(in) :: errorcode
      if (errorcode /= NF90_NOERR) then
          write(*, '(A)') nf90_strerror(errorcode)
          stop 1
      end if
  end subroutine check


  function der2_center(y, dx) result(r)
    ! ***  This calculates the  2nd order central difference of the second
    ! *** derivative of our function
    real, intent(in) :: y(:)
    real, intent(in) :: dx
    real :: r(size(y))
    integer :: i
    
    associate(n => size(y))
      r(2:n-1) = ( y(3:n) - y(1:n-2) )/(2*dx)
      r(1) = ( y(2) - y(1) )/dx
      r(n) = ( y(n) - y(n-1) )/dx
    
    do i=2, n-1
      if (y(i+1) .eq. defval .or. y(i-1) .eq. defval) r(i) = defval
    end do
    if (y(2) .eq. defval .or. y(1) .eq. defval) r(1) = defval
    if (y(n) .eq. defval .or. y(n-1) .eq. defval) r(n) = defval
    end associate

  end function der2_center

  function der2_forward(y, dx) result(r)
    ! ***  This calculates the  2nd order central difference of the second
    ! *** derivative of our function
    real, intent(in) :: y(:)
    real, intent(in) :: dx
    real :: r(size(y))
    integer :: i
    
    associate(n => size(y))
      r(1:n-1) = ( y(2:n) - y(1:n-1) )/dx
    end associate

  end function der2_forward


  function integrate(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real, intent(in)  :: x(:)         !! Variable x
    real, intent(in)  :: y(size(x))   !! Function y(x)
    real              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(2:n) + y(1:n-1))*(x(2:n) - x(1:n-1)))/2
    end associate
    if ( any(y.eq.defval) ) r = defval
  end function integrate

end program main
