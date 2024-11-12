subroutine set_resolution( Mpol, Ntor, Nfp, mn, im, in )

  implicit none

  integer, intent(in)  :: Mpol, Ntor, Nfp, mn
  integer, intent(out) :: im(mn), in(mn)

  integer              :: imn, mm, nn

  imn = 0

  ;  mm = 0
  ;do nn = 0, Ntor
  ; imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
  ;enddo
  ;

  do mm = 1, Mpol
   do nn = -Ntor, Ntor
    imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
   enddo
  enddo

end subroutine set_resolution