program workflow_input_probe
  ! Parse/validate without evaluating the known unsupported all-warm parallel kernel.
  use functions, only: read_in_params,read_scan_input,read_map_input
  use vars, only: nscan,option,use_map
  implicit none
  call read_in_params
  if (nscan>0) call read_scan_input
  if (option==0 .or. option==4 .or. use_map) call read_map_input
  write(*,'(a)') 'INPUT_OK'
end program
