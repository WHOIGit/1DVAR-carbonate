module types
  use netcdf
  implicit none


  integer, parameter :: char_len=16, char_len_long=32

  type ncvar_attrib
     character(len=char_len) :: short_name
     character(len=char_len_long) :: long_name
     character(len=char_len_long) :: units
     integer :: id        
     integer :: type      
  end type ncvar_attrib


  ! NetCDF types
  integer, parameter :: &
       ncf_float=NF90_FLOAT

  real, parameter :: &
       ncf_fillval_float=NF90_FILL_FLOAT


end module types
