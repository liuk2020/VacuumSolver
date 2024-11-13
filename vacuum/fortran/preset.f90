subroutine get_resolution( Mpol, Ntor, Nfp, mn, im, in )

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

end subroutine get_resolution


subroutine get_NAdof(Mpol, Ntor, mn, Lrad, Istellsym, NAdof)
 
  implicit none
  integer, intent(in)  :: Mpol, Ntor, mn, Lrad, Istellsym
  integer, intent(out) :: NAdof 
  logical              :: YESstellsym, NOTstellsym
  integer              :: zerdof
  integer              :: ii, jj

  select case( Istellsym )
  case( 0 )    ; YESstellsym = .false. ; NOTstellsym = .true.
  case( 1 )    ; YESstellsym = .true.  ; NOTstellsym = .false.
  case default ;
  end select

  zerdof = 0                                       ! count Zernike degree of freedom 30 Jun 19
  do ii = 2, Mpol                                  ! for m>1
    do jj = ii, Lrad, 2
    zerdof = zerdof + 2 * Ntor + 1                 ! plus and minus sign for n>1, unique for n==0
    if( NOTstellsym ) zerdof = zerdof + 2*Ntor + 1 ! plus and minus sign for n
    enddo
  enddo
  zerdof = zerdof * 2                              ! we have one for At and one for Az

  do jj = 0, Lrad, 2                                ! for m==0
    zerdof = zerdof + Ntor + 1                      ! minus sign for n, Aze
    if (jj .ge. 2) zerdof = zerdof + Ntor + 1       ! minus sign for n, Ate, without l=0 due to recombination

    if( NOTstellsym ) then
    zerdof = zerdof + Ntor                         ! sin component minus sign for n, Azo
    if (jj .ge. 2) zerdof = zerdof + Ntor          ! minus sign for n, Ato, without l=0 due to recombination
    endif
  enddo

  if (Mpol .ge. 1) then ! for m==1
    do jj = 1, Lrad, 2
      zerdof = zerdof + 2 * Ntor + 1                  ! minus and plus sign for n, Aze
      if (jj .ge. 2) zerdof = zerdof + 2 * Ntor + 1   ! minus sign for n, Ate, without l=0 due to recombination

      if( NOTstellsym ) then
        zerdof = zerdof + 2 * Ntor + 1                 ! sin component minus and plus sign for n, Azo
        if (jj .ge. 2) zerdof = zerdof + 2 * Ntor + 1  ! minus and plus sign for n, Ato, without l=0 due to recombination
      endif
    enddo
  endif
                                    !                                     a    c      b        d      e      f      g   h
  if( YESstellsym ) NAdof = zerdof                               + mn        + Ntor+1        + mn-1        + 1 + 0
  if( NOTstellsym ) NAdof = zerdof                               + mn + mn-1 + Ntor+1 + Ntor + mn-1 + mn-1 + 1 + 0 ! this is broken at the moment

  ! due to basis recombination, Lma will not have the m=0 and m=1 harmonics. We substract them now
  ! m = 0
  NAdof = NAdof - (Ntor + 1)
  if (NOTstellsym) NAdof = NAdof - Ntor

  ! m = 1
  if (Mpol .ge. 1) then
    NAdof = NAdof - (2 * Ntor + 1)
    if (NOTstellsym) NAdof = NAdof - (2 * Ntor + 1)
  endif

end subroutine