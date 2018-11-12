subroutine sysextr(x ,fv, n)
! this system deals with the exogenous incomplete markets

use params


implicit none
integer , intent(in) :: n
real(prec) , intent(in) ::x(n)
real(prec) , intent(out) :: fv(n)
real(prec) :: erhs,sprime,con
real(prec) :: vals(2)
integer inds(2)


 sprime=x(1)


 con=gridm(tyc,rc)+w(tc)*yat(tc,tyc,sc,ic)-sprime/Rm(tc)

	! compute  euler equation

	erhs=0.0
	call basefun (gridm(tyc,:),pm,sprime,vals,inds)
	do scc=1,ns
		do icc=1,ni
			erhs=erhs+probs(tc,sc,scc)*pi(tc+1,icc)*( vals(1)*vpmfun(tc+1,tyc,scc,icc,inds(1))+vals(2)*vpmfun(tc+1,tyc,scc,icc,inds(2)) )
		end do
	end do

	fv(1)=(1.0-bta)*MU(con)-bta*Rm(tc)*erhs

end subroutine sysextr
