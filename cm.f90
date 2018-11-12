subroutine cm

use params
implicit none

real(prec)::cin(ntypes),cint(ntypes,ns,ni),ain(ntypes,ns,ni)
real(prec)::eye(ns*ni,ns*ni),put(ns*ni),mat(ns*ni,ns*ni),matinv(ns*ni,ns*ni),cmexpin(nt,ntypes,ns,ni)
real(prec)::inimean,finmean,inivar,varfin,inistd,stdfin,meanst,meanstt(ntypes),varfint,stdfint
real(prec),dimension(ns*ni)::helper
real(prec),dimension(ns,ni)::vmn

! Compute Initial Wealth Distribution in Complete Markets

tc=1
cin=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			cin(tyc)=cin(tyc)+diste(tc,sc)*pi(tc,ic)*yat(tc,tyc,sc,ic)*(w(tc)+(R(tc)-1.0)*KA(tc))
		end do
	end do
end do

do tyc=1,ntypes
	eye=0.0
	do sc=1,ns
		do ic=1,ni
			put((sc-1)*ni+ic)=cin(tyc)-w(tc)*yat(tc,tyc,sc,ic)
			eye((sc-1)*ni+ic,(sc-1)*ni+ic)=1.0
			do scc=1,ns
				do icc=1,ni
					mat((sc-1)*ni+ic,(scc-1)*ni+icc)=eye((sc-1)*ni+ic,(scc-1)*ni+icc)-bta*probs(tc,sc,scc)*pi(tc,icc)
				end do
			end do
		end do
	end do
	call dlinrg(ni*ns,mat,ni*ns,matinv,ni*ns)
	helper=matmul(matinv,put)
	do sc=1,ns
		ain(tyc,sc,1:ni)=helper((sc-1)*ni+1:sc*ni)
	end do
end do

! Compute Present Discounted Value of Income from Period 2 onwards

cmexpin=0.0

tc=nt
do tyc=1,ntypes
	eye=0.0
	do sc=1,ns
		do ic=1,ni
			put((sc-1)*ni+ic)=w(tc)*yat(tc,tyc,sc,ic)
			eye((sc-1)*ni+ic,(sc-1)*ni+ic)=1.0
			do scc=1,ns
				do icc=1,ni
					mat((sc-1)*ni+ic,(scc-1)*ni+icc)=eye((sc-1)*ni+ic,(scc-1)*ni+icc)-bta*probs(tc,sc,scc)*pi(tc,icc)
				end do
			end do
		end do
	end do
	call dlinrg(ns*ni,mat,ns*ni,matinv,ns*ni)
	helper=matmul(matinv,put)
	do sc=1,ns
		cmexpin(tc,tyc,sc,1:ni)=helper((sc-1)*ni+1:sc*ni)
	end do
end do

do tc=2,nt-1
	ttc=nt+1-tc
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				do scc=1,ns
					do icc=1,ni
						vmn(scc,icc)=bta*probs(ttc,sc,scc)*pi(ttc*1,icc)*cmexpin(ttc+1,tyc,scc,icc)
						cmexpin(ttc,tyc,sc,ic)=w(tc)*yat(ttc,tyc,sc,ic)+sum(vmn)
					end do
				end do
			end do
		end do
	end do
end do

! Compute new consumption levels

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			cint(tyc,sc,ic)=(1.0-bta)*( ain(tyc,sc,ic)+cmexpin(2,tyc,sc,ic) )
		end do
	end do
end do

! Compute Inequality Statistics

meanst=0.0
meanstt=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			meanst=meanst+pitypes(tyc)*cint(tyc,sc,ic)*diste(2,sc)*pi(2,ic)
			meanstt(tyc)=meanstt(tyc)+log(cint(tyc,sc,ic))*diste(2,sc)*pi(2,ic)	
		end do
	end do
end do

inimean=sum(log(cin)*pitypes)
inivar=sum(pitypes*(log(cin)-inimean)**2.0)
inistd=sqrt(inivar)

finmean=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			finmean=finmean+pitypes(tyc)*diste(2,sc)*pi(2,ic)*log(cint(tyc,sc,ic))
		end do
	end do
end do

varfin=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			varfin=varfin+pitypes(tyc)*diste(2,sc)*pi(2,ic)*( (log(cint(tyc,sc,ic))-finmean)**2.0 )
		end do
	end do
end do
stdfin=sqrt(varfin)

varfint=0.0
do tyc=1,ntypes	
	varfint=varfint+pitypes(tyc)*(meanstt(tyc)-finmean)**2.0
end do
stdfint=sqrt(varfint)

end subroutine cm