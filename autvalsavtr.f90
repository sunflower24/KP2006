subroutine autvalsavtr

use params
implicit none

real(prec) ::vmn(ns,ni)
real(prec),dimension(ns)::sol1,gue1
real(prec) :: gue(1),sol(1),vals(2),evn
integer :: inds(2)

! Choose at what rate to save in autarky

external sysextr

select case(choice2)

case(2)

Rm(tc)=Rex

case(3)

Rm(tc)=R(tc)

end select

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do rc=1,pm
				gue=sfun(tc+1,tyc,sc,ic,rc)
       			call dneqnf(sysextr,errel,1,itmax,gue,sol,fnorm)
				if (sol(1) < 0.0 ) then
					sol(1) = 0.0
				end if
				sfun(tc,tyc,sc,ic,rc)=sol(1)
				cmfun(tc,tyc,sc,ic,rc)=gridm(tyc,rc)+w(tc)*yat(tc,tyc,sc,ic)-sfun(tc,tyc,sc,ic,rc)/Rm(tc)
				vpmfun(tc,tyc,sc,ic,rc)=(1.0-bta)*MU(cmfun(tc,tyc,sc,ic,rc))
				
				vmn=0.0
				evn=0.0
				do scc=1,ns
					do icc=1,ni
						call basefun (gridm(tyc,:),pm,sfun(tc,tyc,sc,ic,rc),vals,inds)
						vmn(scc,icc)=probs(tc,sc,scc)*pi(tc+1,icc)*( vals(1)*vmfun(tc+1,tyc,scc,icc,inds(1))+vals(2)*vmfun(tc+1,tyc,scc,icc,inds(2)) )
					end do
				end do
				evn=sum(vmn)
				vmfun(tc,tyc,sc,ic,rc)=(1.0-bta)*U(cmfun(tc,tyc,sc,ic,rc))+bta*evn
			end do
		end do
	end do
end do

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			vmaut(tc,tyc,sc,ic)=vmfun(tc,tyc,sc,ic,1)
		end do
	end do
end do

end subroutine autvalsavtr