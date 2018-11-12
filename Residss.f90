subroutine residss(x,fv)
use params

implicit none

real(prec),intent(in) :: x
real(prec),intent(out):: fv 

R(tc)=x

KA(tc)=LS(tc)*(A*alpha/(R(tc)-1.0+delta))**(1.0/(1.0-alpha))
w(tc)=MPL(KA(tc),LS(tc))

! Compute value of Autarky

if (choice==1) then
	call autval
endif

call dgrid0

! Compute Equilibrium Interest rate

call guess

! Solve for the planners decision rule, given the intertemporal price
call enforce

! compute stationary distribution for given intertemporal price
call diss

print*,'Stationary distribution computed'
! compute excess demand from given stationary distribution
call resource

fv=resd1(tc)

end subroutine residss