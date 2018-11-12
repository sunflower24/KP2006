subroutine autvaltr

use params
implicit none

real(prec),dimension(ns,ni)::vmn,emn

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do scc=1,ns
				do icc=1,ni
					vmn(scc,icc)=bta*probs(tc,sc,scc)*pi(tc+1,icc)*vmaut(tc+1,tyc,scc,icc)
					emn(scc,icc)=(1.0/R(tc))*probs(tc,sc,scc)*pi(tc+1,icc)*expinc(tc+1,tyc,scc,icc)
				end do
			end do
			vmaut(tc,tyc,sc,ic)=(1.0-bta)*U( w(tc)*yat(tc,tyc,sc,ic) )+sum(vmn)
			expinc(tc,tyc,sc,ic)=w(tc)*yat(tc,tyc,sc,ic)+sum(emn)
		end do
	end do
end do

if (choice2 > 1) then
	call autvalsavtr
end if

end subroutine autvaltr