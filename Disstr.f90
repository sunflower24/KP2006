! This is a subroutine to compute the distribution over w,y in different tc

subroutine disstr

use params
IMPLICIT NONE

real(prec), dimension(ns,ni,pa,ns,ni,pa)::TT
real(prec), dimension(ntypes,ns,ni,pa)::enw
integer::acc,inds(2),cc
real(prec),dimension(2)::vals
real(prec),dimension(ni*ns*pa*ntypes)::lorincx,lorincy,lorconsx,lorconsy
real(prec),dimension(nt,ni*ns*pa*ntypes)::lorenzyx,lorenzyy,lorenzcx,lorenzcy
real(prec) :: ginic,giniy
real(prec) :: testdis

! Initialize Distribution from initial steady state. Note that grid has changed between first and second period

omega(2:nt,:,:,:,:)=0.0

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				call basefun(grida(2,tyc,sc,ic,:),pa,grida(1,tyc,sc,ic,ac),vals,inds)
				omega(2,tyc,sc,ic,inds(1))=omega(2,tyc,sc,ic,inds(1))+omega(1,tyc,sc,ic,ac)*vals(1)
				omega(2,tyc,sc,ic,inds(2))=omega(2,tyc,sc,ic,inds(2))+omega(1,tyc,sc,ic,ac)*vals(2)
			end do
		end do
	end do
end do


! Loop to find distributions for t>2

do tc=2,nt-1

do tyc=1,ntypes


TT=0.0
do sc=1,ns
	do ic=1,ni				! Loop over States of productivity
		do ac=1,pa				! Loop over old assets
			do scc=1,ns
				do icc=1,ni
					call basefun(grida(tc+1,tyc,scc,icc,:),pa,afun(tc,tyc,sc,ic,ac,scc,icc),vals,inds)
					TT(sc,ic,ac,scc,icc,inds(1))=vals(1)*probs(tc,sc,scc)*pi(tc+1,icc)
					TT(sc,ic,ac,scc,icc,inds(2))=vals(2)*probs(tc,sc,scc)*pi(tc+1,icc)
				end do
			end do
		end do
	end do
end do

omega(tc+1,tyc,1:ns,1:ni,1:pa)=0.0
do scc=1,ns
	do icc=1,ni
		do acc=1,pa
			do sc=1,ns
				do ic=1,ni
					do ac=1,pa
						omega(tc+1,tyc,scc,icc,acc)=omega(tc+1,tyc,scc,icc,acc)+TT(sc,ic,ac,scc,icc,acc)*omega(tc,tyc,sc,ic,ac)
					end do
				end do
			end do
		end do
	end do
end do

end do

! Test whether distributions sum to 1

testdis=sum(omega(tc,1:ntypes,1:ns,1:ni,1:pa))
if ( (testdis < 0.9999) .or. (testdis > 1.0001 ) ) then
	print*,'Distribution does not sum to 1 '
	print*,tc, testdis
	pause
end if

! Average Consumption

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				enw(tyc,sc,ic,ac)=yat(tc,tyc,sc,ic)*omega(tc,tyc,sc,ic,ac)
			end do
		end do
	end do
end do

cons(tc)=sum(omega(tc,1:ntypes,1:ns,1:ni,1:pa)*cfun(tc,1:ntypes,1:ns,1:ni,1:pa))  ! total consumption
resd(tc)=(cons(tc)+KA(tc+1)-(1.0-delta)*KA(tc)-F(KA(tc),LS(tc)))/KA(tc)

call condist

cons1(tc)=0.0
do tyc=1,ntypes
	do cc=1,pc
		do sc=1,ns
			do ic=1,ni
				cons1(tc)=cons1(tc)+cdist(tc,tyc,cc,sc,ic)*gridc(cc)
			end do
		end do
	end do
end do

end do

do tc=1,nt

assdemand(tc)=0.0
assdemand(tc)= sum(omega(tc,1:ntypes,1:ns,1:ni,1:pa)*grida(tc,1:ntypes,1:ns,1:ni,1:pa))

end do
assdemand1(2)=assdemand1(1)
do tc=1,nt-1
	assdemand1(tc+1)=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				do scc=1,ns
					do icc=1,ni
						assdemand1(tc+1)=assdemand1(tc+1)+probs(tc,sc,scc)*pi(tc+1,icc)*omega(tc,tyc,sc,ic,ac)*afun(tc,tyc,sc,ic,ac,scc,icc)
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
end do

resd1(2)=resd1(1)
do tc=3,nt
resd1(tc)=(-assdemand(tc)+R(tc-1)*KA(tc))/KA(tc)
end do


! Restore final distribution
omega(nt,1:ntypes,1:ns,1:ni,1:pa)=finomega(1:ntypes,1:ns,1:ni,1:pa)

open(unit=80,file='residual.txt')

do tc=1,nt
	write(80,'(6f9.6)') resd1(tc),resd(tc),assdemand(tc),KA(tc),tote(tc),cons1(tc),cons(tc)
end do

close(80)

do tc=1,nt
	call gini2(omega(tc,1:ntypes,1:ns,1:ni,1:pa),lorincx,lorincy,lorconsx,lorconsy,giniy,ginic)
	giniyvec(tc)=giniy
	ginicvec(tc)=ginic
	lorenzyx(tc,1:pa*ns*ni*ntypes)=lorincx
	lorenzyy(tc,1:pa*ns*ni*ntypes)=lorincy
	lorenzcx(tc,1:pa*ns*ni*ntypes)=lorconsx
	lorenzcy(tc,1:pa*ns*ni*ntypes)=lorconsy
end do

totint=0.0
do tc=1,nt
  do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				if (grida(tc,tyc,sc,ic,ac)<=0.0) then
					totint(tc)=totint(tc)+(omega(tc,tyc,sc,ic,ac)*abs(grida(tc,tyc,sc,ic,ac)))
				end if
			end do
		end do
	end do
  end do
end do 

end Subroutine disstr
