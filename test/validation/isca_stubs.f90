! Minimal isolated harness support; constants match JGCM defaults.
module constants_mod
  implicit none
  real, parameter :: rdgas=287.04, rvgas=461.50, kappa=2.0/7.0
  real, parameter :: Cp_air=rdgas/kappa, HLv=2.5e6, HLs=2.834e6, Grav=9.80
end module

module fms_mod
  implicit none
  integer, parameter :: FATAL=2
contains
  logical function file_exist(file)
    character(len=*),intent(in)::file
    inquire(file=file,exist=file_exist)
  end function
  subroutine error_mesg(module_name,message,severity)
    character(len=*),intent(in)::module_name,message
    integer,intent(in)::severity
    print *,trim(module_name),trim(message),severity
    error stop
  end subroutine
  integer function open_file(file,action)
    character(len=*),intent(in)::file,action
    if(action=='append')then
      open(newunit=open_file,file=file,position='append')
    else
      open(newunit=open_file,file=file,action=action)
    endif
  end function
  integer function check_nml_error(io,name)
    integer,intent(in)::io
    character(len=*),intent(in)::name
    check_nml_error=io
  end function
  integer function mpp_pe()
    mpp_pe=0
  end function
  subroutine close_file(unit)
    integer,intent(in)::unit
    close(unit)
  end subroutine
end module
