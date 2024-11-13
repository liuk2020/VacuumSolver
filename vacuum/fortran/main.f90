program test_dofs

  implicit none
  
  integer :: Mpol, Ntor, Lrad, mn, NAdof

  print *, 'Mpol:'
  read *, Mpol
  print *, 'Ntor:'
  read *, Ntor
  print *, 'Lrad:'
  read *, Lrad

  mn = 1 + Ntor + Mpol*(2*Ntor+1)

  call get_NAdof(Mpol, Ntor, mn, Lrad, .true., NAdof)

  print *, 'NAdof:', NAdof

end program test_dofs