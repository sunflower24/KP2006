subroutine valuecm(x ,fv, n)
! this system deals with the case in which highest DC is binding

use params


implicit none
integer , intent(in) :: n
real(prec) , intent(in) ::x(n)
real(prec) , intent(out) :: fv(n)
real(prec)::consum,val

! Consumption in complete markets in first state

consum=(1.0- ( (R(tc)*bta)**(1.0/sigma) )/R(tc) )*( x(1)+expinc(tc,tyc,sc,ic) )

if (abs(sigma-1.0) > 0.0001) then
	val=( (1.0-bta)/(1.0-bta*(bta*R(tc))**((1.0-sigma)/sigma) ) )*U(consum)
else
	val=U(consum)+(bta/(1.0-bta))*log(bta*R(tc))
end if

fv(1)=val-vmaut(tc,tyc,sc,ic)

end subroutine valuecm
