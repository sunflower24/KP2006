! This subroutine compute the Initial Steady State

subroutine FinSS

use params
implicit none
real(prec),dimension(ns*ni*pa*ntypes)::lorincx,lorincy,lorconsx,lorconsy
real(prec) :: giniinc,ginicons

external residss


! Initialize the time counter

tc=nt
R(tc)=1.04

! Redefine Grid

call dgrid0

! Compute Equilibrium Interest rate

call guess

call Newtonss(residss,R(tc))

! Save Final Distribution 
finomega=omega(nt,1:ntypes,1:ns,1:ni,1:pa)

call gini2(finomega,lorincx,lorincy,lorconsx,lorconsy,giniinc,ginicons)
giniyvec(tc)=giniinc
ginicvec(tc)=ginicons

!pause

end subroutine FinSS