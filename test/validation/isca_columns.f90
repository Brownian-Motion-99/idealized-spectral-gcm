! Reads fixed input columns; uses the original public Isca SBM API.
program compare_columns
  use qe_moist_convection_mod
  use sat_vapor_pres_mod, only: escomp
  use constants_mod, only: rdgas, rvgas
  implicit none
  integer :: n, count, unit_in, unit_out, k, i, enabled
  real, allocatable :: t(:,:,:), q(:,:,:), pf(:,:,:), ph(:,:,:)
  real, allocatable :: dt(:,:,:), dq(:,:,:), tr(:,:,:), qr(:,:,:)
  real :: rain(1,1), snow(1,1), cape(1,1), cin(1,1), iq(1,1), it(1,1)
  integer :: flag(1,1), lzb(1,1), lcl(1,1)
  real :: origin_es, origin_r, origin_rs
  logical :: cold(1,1)
  character(len=4096) :: input_path, output_path
  call get_command_argument(1, input_path)
  call get_command_argument(2, output_path)
  open(newunit=unit_in, file=trim(input_path), status='old')
  open(newunit=unit_out, file=trim(output_path), status='replace')
  read(unit_in,*) count, n
  allocate(t(1,1,n),q(1,1,n),pf(1,1,n),ph(1,1,n+1))
  allocate(dt(1,1,n),dq(1,1,n),tr(1,1,n),qr(1,1,n))
  cold=.false.
  call qe_moist_convection_init()
  do i=1,count
    read(unit_in,*) enabled
    do k=1,n
      read(unit_in,*) pf(1,1,k),ph(1,1,k),ph(1,1,k+1),t(1,1,k),q(1,1,k)
    enddo
    if (enabled == 0) cycle
    call qe_moist_convection(1.,t,q,pf,ph,cold,rain,snow,dt,dq,qr,flag,lzb,cape,cin,iq,it,tr,lcl)
    call escomp(t(1,1,n),origin_es)
    origin_r=q(1,1,n)/(1.-q(1,1,n))
    origin_rs=rdgas*origin_es/rvgas/(pf(1,1,n)-origin_es)
    write(unit_out,*) i,flag(1,1),lcl(1,1),lzb(1,1),cape(1,1),cin(1,1),rain(1,1),origin_r-origin_rs
    do k=1,n
      write(unit_out,*) dt(1,1,k),dq(1,1,k),tr(1,1,k),qr(1,1,k)
    enddo
  enddo
  call qe_moist_convection_end()
  close(unit_in)
  close(unit_out)
end program
