subroutine resource

use params
implicit none

real(prec),dimension(ntypes,ns,ni,pa)::distr1
real(prec)::privint
real(prec),dimension(ntypes,ns,ni,pa)::enw
real(prec)::testdis
real(prec)::giniinc,ginicons
real(prec),dimension(ni*ns*pa*ntypes)::lorincx,lorincy,lorconsx,lorconsy
privint=0.0

print*, "resource started"

! unstack distribution over (y,w)

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			distr1(tyc,sc,ic,1:pa)=distr(tc,tyc,((sc-1)*ni+ic-1)*pa+1:((sc-1)*ni+ic)*pa)
		end do
	end do
end do

! construct total resources available in the economy

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa	
				enw(tyc,sc,ic,ac)=yat(tc,tyc,sc,ic)*distr1(tyc,sc,ic,ac)
				privint=privint+abs(cfun(tc,tyc,sc,ic,ac)-w(tc)*yat(tc,tyc,sc,ic))*distr1(tyc,sc,ic,ac)
			end do
		end do
	end do
end do

tote(tc)=sum(enw(1:ntypes,1:ns,1:ni,1:pa))

cons(tc)=sum(distr1*cfun(tc,1:ntypes,1:ns,1:ni,1:pa))  ! total consumption
omega(tc,1:ntypes,1:ns,1:ni,1:pa)=distr1

testdis=sum(omega(tc,1:ntypes,1:ns,1:ni,1:pa))
if ( (testdis < 0.9999) .or. (testdis > 1.0001 ) ) then
	print*,'Distribution does not sum to 1 '
	print*,tc, testdis
	pause
end if

resd(tc)=(cons(tc)+delta*KA(tc)-F(KA(tc),LS(tc)))/KA(tc)

assdemand(tc)=0.0
assdemand(tc)= sum(omega(tc,1:ntypes,1:ns,1:ni,1:pa)*grida(tc,1:ntypes,1:ns,1:ni,1:pa))

assdemand1(tc)=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				do scc=1,ns
					do icc=1,ni
						assdemand1(tc)=assdemand1(tc)+probs(tc,sc,scc)*pi(tc,icc)*omega(tc,tyc,sc,ic,ac)*afun(tc,tyc,sc,ic,ac,scc,icc)
					end do
				end do
			end do
		end do
	end do
end do

if ( abs( assdemand(tc)-assdemand1(tc) )>0.0001 ) then
		print*,'Oh, oh, aggregation fails',tc,assdemand(tc),assdemand1(tc) 
		pause
endif

resd1(tc)=(-assdemand(tc)+R(tc)*KA(tc))/KA(tc)

print*,'Excess Demand in the Goods and the Asset Market are '
print*,resd(tc),resd1(tc)
print*,'Interest Rate in Percent is ',100*(R(tc)-1)
print*,'Interest Rate Target is ',100*rtarget
print*,'Wealth to Output Ratio is ',KA(tc)/F(KA(tc),LS(tc))
print*,'Target is ',rat

call condist

call gini2(omega(tc,1:ntypes,1:ns,1:ni,1:pa),lorincx,lorincy,lorconsx,lorconsy,giniinc,ginicons)
!print*,giniinc,ginicons

end subroutine resource