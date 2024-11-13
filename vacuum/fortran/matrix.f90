subroutine matrix( lvol, mn, lrad )
  use mod_kinds, only: wp => dp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  use constants, only : zero, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wmatrix, mpol

  use cputiming, only : Tmatrix

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        YESstellsym, NOTstellsym, &
                        im, in, &
                        NAdof, &
                        dMA, dMD, dMB, dMG, &
                        Ate, Ato, Aze, Azo, &
                        iVns, iBns, iVnc, iBnc, &
                        Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                        Lcoordinatesingularity, TT, RTT, RTM, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        dBdX

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


#ifdef OPENMP
  USE OMP_LIB
#endif
  use mpi
  implicit none
  integer   :: ierr, astat, ios, nthreads, ithread
  real(wp)      :: cput, cpui, cpuo=0 ! cpu time; cpu initial; cpu old; 31 Jan 13;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  integer, intent(in)  :: lvol, mn, lrad

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  integer              :: NN, ii, jj, ll, kk, pp, ll1, pp1, mi, ni, mj, nj, mimj, minj, nimj, ninj, mjmi, mjni, njmi, njni, id, jd

  real(wp)                 :: Wtete, Wteto, Wtote, Wtoto
  real(wp)                 :: Wteze, Wtezo, Wtoze, Wtozo
  real(wp)                 :: Wzete, Wzeto, Wzote, Wzoto
  real(wp)                 :: Wzeze, Wzezo, Wzoze, Wzozo

  real(wp)                 :: Htete, Hteto, Htote, Htoto
  real(wp)                 :: Hteze, Htezo, Htoze, Htozo
  real(wp)                 :: Hzete, Hzeto, Hzote, Hzoto
  real(wp)                 :: Hzeze, Hzezo, Hzoze, Hzozo

  real(wp),allocatable     :: TTdata(:,:,:), TTMdata(:,:) ! queues to construct sparse matrices


  cpui = MPI_WTIME()
  cpuo = cpui
#ifdef OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

   if( .not.allocated(dMA) ) then
     write(6,'("matrix :      fatal : myid=",i3," ; .not.allocated(dMA) ; error ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "matrix : .not.allocated(dMA) : error  ;"
    endif


   if( .not.allocated(dMD) ) then
     write(6,'("matrix :      fatal : myid=",i3," ; .not.allocated(dMD) ; error ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "matrix : .not.allocated(dMD) : error  ;"
    endif


   if( .not.allocated(dMB) ) then
     write(6,'("matrix :      fatal : myid=",i3," ; .not.allocated(dMB) ; error ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "matrix : .not.allocated(dMB) : error  ;"
    endif


   if( .not.allocated(dMG) ) then
     write(6,'("matrix :      fatal : myid=",i3," ; .not.allocated(dMG) ; error ;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "matrix : .not.allocated(dMG) : error  ;"
    endif

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  NN = NAdof(lvol) ! shorthand;

  dMA(0:NN,0:NN) = zero
  dMD(0:NN,0:NN) = zero


   allocate( TTdata(0:lrad, 0:mpol, 0:1), stat=astat )
   TTdata(0:lrad, 0:mpol, 0:1) = zero


   allocate( TTMdata(0:lrad, 0:mpol), stat=astat )
   TTMdata(0:lrad, 0:mpol) = zero


  ! fill in Zernike/Chebyshev polynomials depending on Lcooridnatesingularity
  if (Lcoordinatesingularity) then
    TTdata(0:lrad,0:mpol,0:1) = RTT(0:lrad,0:mpol,0:1,0)
    TTMdata(0:lrad,0:mpol) = RTM(0:lrad,0:mpol)
  else
    do ii = 0, mpol
      TTdata(0:lrad,ii,0:1) = TT(0:lrad,0:1,0)
      TTMdata(0:lrad,ii) = TT(0:lrad,0,0)
    enddo
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( YESstellsym ) then
!$OMP PARALLEL DO PRIVATE(ii,jj,ll,pp,mi,ni,mj,nj,mimj,minj,nimj,ninj,ll1,pp1,Wtete,Wzete,Wteze,Wzeze,Htete,Hzete,Hteze,Hzeze,id,jd,kk) SHARED(lvol,lrad)
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mimj = mi * mj ; minj = mi * nj ; nimj = ni * mj ; ninj = ni * nj

     do ll = 0, lrad

      do pp = 0, lrad

       if (Lcoordinatesingularity) then
        if (ll < mi .or. pp < mj) cycle ! rule out zero components of Zernike; 02 Jul 19
        if (mod(ll+mi,2)+mod(pp+mj,2)>0) cycle ! rule out zero components of Zernike; 02 Jul 19
        ll1 = (ll - mod(ll, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
        pp1 = (pp - mod(pp, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
       else
        ll1 = ll
        pp1 = pp
       end if

       Wtete = + 2 * ninj * TTssss(ll1,pp1,ii,jj) - 2 * ni      * TDszsc(ll1,pp1,ii,jj) - 2      * nj * TDszsc(pp1,ll1,jj,ii) + 2 * DDzzcc(ll1,pp1,ii,jj)
       Wzete = + 2 * nimj * TTssss(ll1,pp1,ii,jj) + 2 * ni      * TDstsc(ll1,pp1,ii,jj) - 2      * mj * TDszsc(pp1,ll1,jj,ii) - 2 * DDtzcc(pp1,ll1,jj,ii)
       Wteze = + 2 * minj * TTssss(ll1,pp1,ii,jj) + 2      * nj * TDstsc(pp1,ll1,jj,ii) - 2 * mi      * TDszsc(ll1,pp1,ii,jj) - 2 * DDtzcc(ll1,pp1,ii,jj)
       Wzeze = + 2 * mimj * TTssss(ll1,pp1,ii,jj) + 2 * mi      * TDstsc(ll1,pp1,ii,jj) + 2      * mj * TDstsc(pp1,ll1,jj,ii) + 2 * DDttcc(ll1,pp1,ii,jj)

       Htete =   zero
       Hzete = - DToocc(pp1,ll1,jj,ii) + DToocc(ll1,pp1,ii,jj)
       Hteze = + DToocc(pp1,ll1,jj,ii) - DToocc(ll1,pp1,ii,jj)
       Hzeze =   zero

       id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtete ; dMD(id,jd) = Htete
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzete ; dMD(id,jd) = Hzete
       id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteze ; dMD(id,jd) = Hteze
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeze ; dMD(id,jd) = Hzeze

      enddo ! end of do pp ;

     enddo ! end of do ll ;

    enddo ! end of do jj ;

    if (dBdX%L) cycle

    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif

    do ll = 0, lrad
    ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,  ii)       ; dMA(id,jd) = +      TTMdata(ll, mi)
    ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,  ii)       ; dMA(id,jd) = +      TTdata(ll, mi,kk)
    if( ii.gt.1 ) then ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - ni * TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,  ii)       ; dMA(id,jd) = - mi * TTdata(ll, mi, 1)
    else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,  ii)       ; dMA(id,jd) = +      TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,  ii)       ; dMA(id,jd) = +      TTdata(ll, mi, 1)
    endif

    enddo ! end of do ll ;

    do pp = 0, lrad     ; id = Lma(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTMdata(pp, mi)
    ;                  ; id = Lmb(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi,kk)
    if( ii.gt.1 ) then ; id = Lme(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TTdata(pp, mi, 1)
      ;                 ; id = Lme(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TTdata(pp, mi, 1)
    else               ; id = Lmg(lvol,  ii)       ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
      ;                 ; id = Lmh(lvol,  ii)       ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
    endif
    enddo ! end of do pp ;

   enddo ! end of do ii ;
!$OMP END PARALLEL DO

  else ! NOTstellsym ;

!$OMP PARALLEL DO PRIVATE(ii,jj,ll,pp,mi,ni,mj,nj,mjmi,mjni,njmi,njni,ll1,pp1,Wtete,Wzete,Wteze,Wzeze,Htete,Hzete,Hteze,Hzeze,Wteto,Wzeto,Wtezo,Wzezo,Hteto,Hzeto,Htezo,Hzezo,Wtote,Wzote,Wtoze,Wzoze,Htote,Hzote,Htoze,Hzoze,Wtoto,Wzoto,Wtozo,Wzozo,Htoto,Hzoto,Htozo,Hzozo,id,jd,kk) SHARED(lvol,lrad)
   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

    do jj = 1, mn ; mj = im(jj) ; nj = in(jj) ; mjmi = mi * mj ; mjni = ni * mj ; njmi = mi * nj ; njni = ni * nj

     do ll = 0, lrad

      do pp = 0, lrad

       if (Lcoordinatesingularity) then
        if (ll < mi .or. pp < mj) cycle ! rule out zero components of Zernike; 02 Jul 19
        if (mod(ll+mi,2)+mod(pp+mj,2)>0) cycle ! rule out zero components of Zernike; 02 Jul 19
        ll1 = (ll - mod(ll, 2)) / 2! shrinked dof for Zernike; 02 Jul 19
        pp1 = (pp - mod(pp, 2)) / 2 ! shrinked dof for Zernike; 02 Jul 19
       else
        ll1 = ll
        pp1 = pp
       end if

       Wtete = 2 * ( + njni * TTssss(pp1,ll1,jj,ii) - nj * TDszsc(pp1,ll1,jj,ii) - ni * TDszsc(ll1,pp1,ii,jj) + DDzzcc(pp1,ll1,jj,ii) )
       Wtote = 2 * ( - njni * TTsscs(pp1,ll1,jj,ii) + nj * TDszcc(pp1,ll1,jj,ii) - ni * TDszss(ll1,pp1,ii,jj) + DDzzsc(pp1,ll1,jj,ii) )
       Wzete = 2 * ( + mjni * TTssss(pp1,ll1,jj,ii) - mj * TDszsc(pp1,ll1,jj,ii) + ni * TDstsc(ll1,pp1,ii,jj) - DDtzcc(pp1,ll1,jj,ii) )
       Wzote = 2 * ( - mjni * TTsscs(pp1,ll1,jj,ii) + mj * TDszcc(pp1,ll1,jj,ii) + ni * TDstss(ll1,pp1,ii,jj) - DDtzsc(pp1,ll1,jj,ii) )

       Wteto = 2 * ( - njni * TTsssc(pp1,ll1,jj,ii) - nj * TDszss(pp1,ll1,jj,ii) + ni * TDszcc(ll1,pp1,ii,jj) + DDzzcs(pp1,ll1,jj,ii) )
       Wtoto = 2 * ( + njni * TTsscc(pp1,ll1,jj,ii) + nj * TDszcs(pp1,ll1,jj,ii) + ni * TDszcs(ll1,pp1,ii,jj) + DDzzss(pp1,ll1,jj,ii) )
       Wzeto = 2 * ( - mjni * TTsssc(pp1,ll1,jj,ii) - mj * TDszss(pp1,ll1,jj,ii) - ni * TDstcc(ll1,pp1,ii,jj) - DDtzcs(pp1,ll1,jj,ii) )
       Wzoto = 2 * ( + mjni * TTsscc(pp1,ll1,jj,ii) + mj * TDszcs(pp1,ll1,jj,ii) - ni * TDstcs(ll1,pp1,ii,jj) - DDtzss(pp1,ll1,jj,ii) )

       Wteze = 2 * ( + njmi * TTssss(pp1,ll1,jj,ii) + nj * TDstsc(pp1,ll1,jj,ii) - mi * TDszsc(ll1,pp1,ii,jj) - DDtzcc(pp1,ll1,jj,ii) )
       Wtoze = 2 * ( - njmi * TTsscs(pp1,ll1,jj,ii) - nj * TDstcc(pp1,ll1,jj,ii) - mi * TDszss(ll1,pp1,ii,jj) - DDtzsc(pp1,ll1,jj,ii) )
       Wzeze = 2 * ( + mjmi * TTssss(pp1,ll1,jj,ii) + mj * TDstsc(pp1,ll1,jj,ii) + mi * TDstsc(ll1,pp1,ii,jj) + DDttcc(pp1,ll1,jj,ii) )
       Wzoze = 2 * ( - mjmi * TTsscs(pp1,ll1,jj,ii) - mj * TDstcc(pp1,ll1,jj,ii) + mi * TDstss(ll1,pp1,ii,jj) + DDttsc(pp1,ll1,jj,ii) )

       Wtezo = 2 * ( - njmi * TTsssc(pp1,ll1,jj,ii) + nj * TDstss(pp1,ll1,jj,ii) + mi * TDszcc(ll1,pp1,ii,jj) - DDtzcs(pp1,ll1,jj,ii) )
       Wtozo = 2 * ( + njmi * TTsscc(pp1,ll1,jj,ii) - nj * TDstcs(pp1,ll1,jj,ii) + mi * TDszcs(ll1,pp1,ii,jj) - DDtzss(pp1,ll1,jj,ii) )
       Wzezo = 2 * ( - mjmi * TTsssc(pp1,ll1,jj,ii) + mj * TDstss(pp1,ll1,jj,ii) - mi * TDstcc(ll1,pp1,ii,jj) + DDttcs(pp1,ll1,jj,ii) )
       Wzozo = 2 * ( + mjmi * TTsscc(pp1,ll1,jj,ii) - mj * TDstcs(pp1,ll1,jj,ii) - mi * TDstcs(ll1,pp1,ii,jj) + DDttss(pp1,ll1,jj,ii) )

       Htete =   zero
       Htote =   zero
       Hzete = - DToocc(pp1,ll1,jj,ii) + DToocc(ll1,pp1,ii,jj)
       Hzote = - DToosc(pp1,ll1,jj,ii) + DToocs(ll1,pp1,ii,jj)

       Hteto =   zero
       Htoto =   zero
       Hzeto = - DToocs(pp1,ll1,jj,ii) + DToosc(ll1,pp1,ii,jj)
       Hzoto = - DTooss(pp1,ll1,jj,ii) + DTooss(ll1,pp1,ii,jj)

       Hteze = + DToocc(pp1,ll1,jj,ii) - DToocc(ll1,pp1,ii,jj)
       Htoze = + DToosc(pp1,ll1,jj,ii) - DToocs(ll1,pp1,ii,jj)
       Hzeze =   zero
       Hzoze =   zero

       Htezo = + DToocs(pp1,ll1,jj,ii) - DToosc(ll1,pp1,ii,jj)
       Htozo = + DTooss(pp1,ll1,jj,ii) - DTooss(ll1,pp1,ii,jj)
       Hzezo =   zero
       Hzozo =   zero

       id = Ate(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtete ; dMD(id,jd) = Htete
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtote ; dMD(id,jd) = Htote
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzete ; dMD(id,jd) = Hzete
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzote ; dMD(id,jd) = Hzote
       id = Ato(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteto ; dMD(id,jd) = Hteto
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtoto ; dMD(id,jd) = Htoto
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeto ; dMD(id,jd) = Hzeto
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzoto ; dMD(id,jd) = Hzoto
       id = Aze(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wteze ; dMD(id,jd) = Hteze
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtoze ; dMD(id,jd) = Htoze
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzeze ; dMD(id,jd) = Hzeze
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzoze ; dMD(id,jd) = Hzoze
       id = Azo(lvol,0,ii)%i(ll) ; jd = Ate(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtezo ; dMD(id,jd) = Htezo
       ;                         ; jd = Ato(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wtozo ; dMD(id,jd) = Htozo
       ;                         ; jd = Aze(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzezo ; dMD(id,jd) = Hzezo
       ;                         ; jd = Azo(lvol,0,jj)%i(pp) ; dMA(id,jd) = Wzozo ; dMD(id,jd) = Hzozo

      enddo ! end of do pp ;

     enddo ! end of do jj ;

    enddo ! end of do ll ;

    if (dBdX%L) cycle

    if( Lcoordinatesingularity .and. ii.eq.1 ) then ; kk = 1
    else                                            ; kk = 0
    endif

    do ll = 0, lrad ! Chebyshev polynomial ;

    ;                  ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lma(lvol,ii)         ; dMA(id,jd) = TTMdata(ll, mi)
    ;                  ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmb(lvol,ii)         ; dMA(id,jd) = TTdata(ll, mi,kk)
    if( ii.gt.1 ) then ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmc(lvol,ii)         ; dMA(id,jd) = TTMdata(ll, mi)
      ;                 ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmd(lvol,ii)         ; dMA(id,jd) = TTdata(ll, mi,0)
      ;                 ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lme(lvol,ii)         ; dMA(id,jd) = - ni * TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lme(lvol,ii)         ; dMA(id,jd) = - mi * TTdata(ll, mi, 1)
      ;                 ; id = Ato(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,ii)         ; dMA(id,jd) = + ni * TTdata(ll, mi, 1)
      ;                 ; id = Azo(lvol,0,ii)%i(ll) ; jd = Lmf(lvol,ii)         ; dMA(id,jd) = + mi * TTdata(ll, mi, 1)
    else               ; id = Ate(lvol,0,ii)%i(ll) ; jd = Lmg(lvol,ii)         ; dMA(id,jd) = +      TTdata(ll, mi, 1)
      ;                 ; id = Aze(lvol,0,ii)%i(ll) ; jd = Lmh(lvol,ii)         ; dMA(id,jd) = +      TTdata(ll, mi, 1)
    endif

    enddo ! end of do ll;

    do pp = 0, lrad
    ;                  ; id = Lma(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTMdata(pp, mi)
    ;                  ; id = Lmb(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTdata(pp, mi,kk)
    if( ii.gt.1 ) then ; id = Lmc(lvol,ii)         ; jd = Ato(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTMdata(pp, mi)
      ;                 ; id = Lmd(lvol,ii)         ; jd = Azo(lvol,0,ii)%i(pp) ; dMA(id,jd) = TTdata(pp, mi,0)
      ;                 ; id = Lme(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = - ni * TTdata(pp, mi, 1)
      ;                 ; id = Lme(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = - mi * TTdata(pp, mi, 1)
      ;                 ; id = Lmf(lvol,ii)         ; jd = Ato(lvol,0,ii)%i(pp) ; dMA(id,jd) = + ni * TTdata(pp, mi, 1)
      ;                 ; id = Lmf(lvol,ii)         ; jd = Azo(lvol,0,ii)%i(pp) ; dMA(id,jd) = + mi * TTdata(pp, mi, 1)
    else               ; id = Lmg(lvol,ii)         ; jd = Ate(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
      ;                 ; id = Lmh(lvol,ii)         ; jd = Aze(lvol,0,ii)%i(pp) ; dMA(id,jd) = +      TTdata(pp, mi, 1)
    endif
    enddo ! end of do pp ;

   enddo ! end of do ii ;
!$OMP END PARALLEL DO
  endif ! end of if( YESstellsym ) ;

  ! call subroutine matrixBG to construct dMB and dMG

   cput = MPI_WTIME()
   Tmatrix = Tmatrix + ( cput-cpuo )
   call matrixBG( lvol, mn, lrad )
   cpuo = MPI_WTIME()



   deallocate(TTdata ,stat=astat)


   deallocate(TTMdata ,stat=astat)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef DEBUG

  if( Wmatrix ) then ! check symmetry of matrices ;

   do ii = 1, NN

    do jj = 1, NN
     if( abs(dMA(ii,jj)-dMA(jj,ii)) .gt. small*abs(dMA(ii,jj)+dMA(jj,ii)) ) then
      write(ounit,1000) myid, dMA(ii,jj), dMA(jj,ii), dMA(ii,jj)-dMA(jj,ii)
     endif

     if( abs(dMD(ii,jj)-dMD(jj,ii)) .gt. small*abs(dMD(ii,jj)+dMD(jj,ii)) ) then
      write(ounit,1001) myid, dMD(ii,jj), dMD(jj,ii), dMD(ii,jj)-dMD(jj,ii)
     endif

    enddo ! end of do jj;

   enddo ! end of do ii;

  endif ! end of if( Wmatrix ) ;

1000 format("matrix : " 10x " : myid="i3" : dMA(ii,jj)="es23.15", dMA(jj,ii)="es23.15", dMA(ii,jj)-dMA(jj,ii)="es13.5" ;")
1001 format("matrix : " 10x " : myid="i3" : dMD(ii,jj)="es23.15", dMD(jj,ii)="es23.15", dMD(ii,jj)-dMD(jj,ii)="es13.5" ;")

#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


9999 continue
  cput = MPI_WTIME()
  Tmatrix = Tmatrix + ( cput-cpuo )
  return


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine matrix