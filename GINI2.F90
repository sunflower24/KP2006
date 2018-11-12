subroutine gini2(distrib,lorincx,lorincy,lorconsx,lorconsy,giniinc,ginicons)

use params
use dflib
implicit none


integer i,j,index,hh
real(prec)::totinc,totcons,mass
real(prec),dimension(ntypes,ns,ni,pa),intent(in) :: distrib
real(prec),dimension(ntypes,ns,ni,pa)::inc,distribu
real(prec),dimension(ni*ns*pa*ntypes,2)::incmat,incmat1,consmat,consmat1
real(prec),dimension(ni*ns*pa*ntypes)::incv,consv
real(prec),dimension(ni*ns*pa*ntypes),intent(out)::lorincx,lorincy,lorconsx,lorconsy
real(prec),intent(out)::giniinc,ginicons
real(prec),parameter::epsil=0.000000001
integer,dimension(ni*ns*pa*ntypes)::perm
integer::indexh


! type inizialization  (only valid for two types)

hh=ni*ns*pa*ntypes
do i=1,hh
	perm(i)=i
end do

! Rescale the distribution to allow for types  (there is equal mass of each type)
distribu=distrib !/real(ntypes)

! Generate pre-tax and after tax income

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				inc(tyc,sc,ic,ac)=w(tc)*yat(tc,tyc,sc,ic)
			end do
		end do
	end do
end do


! Stack matrix as vector
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			indexh=ns*ni*(tyc-1)+(sc-1)*ni+ic
			incmat((indexh-1)*pa+1:indexh*pa,1)=inc(tyc,sc,ic,1:pa)
			incmat((indexh-1)*pa+1:indexh*pa,2)=distribu(tyc,sc,ic,1:pa)
			consmat((indexh-1)*pa+1:indexh*pa,1)=cfun(tc,tyc,sc,ic,1:pa)
			consmat((indexh-1)*pa+1:indexh*pa,2)=distribu(tyc,sc,ic,1:pa)
		end do
	end do
end do


! first sort matrices according to first dimension


do i=1,hh
	perm(i)=i
end do
call dsvrgp(hh,incmat(1:hh,1),incmat1(1:hh,1),perm)
call dpermu(hh,incmat(1:hh,2),perm,1,incmat1(1:hh,2))

do i=1,hh
	perm(i)=i
end do
call dsvrgp(hh,consmat(1:hh,1),consmat1(1:hh,1),perm)
call dpermu(hh,consmat(1:hh,2),perm,1,consmat1(1:hh,2))


! compute the lorenz curves


totinc=sum(incmat1(1:hh,1)*incmat1(1:hh,2))
totcons=sum(consmat1(1:hh,1)*consmat1(1:hh,2))


do i=1,hh
	mass=sum(incmat1(1:i,1)*incmat1(1:i,2))
	lorincy(i)=mass/totinc
	lorincx(i)=sum(incmat1(1:i,2))

	mass=sum(consmat1(1:i,1)*consmat1(1:i,2))
	lorconsy(i)=mass/totcons
	lorconsx(i)=sum(consmat1(1:i,2))

end do

giniinc=0.0
ginicons=0.0

do i=1,hh
	do j=1,hh
		giniinc=giniinc+incmat1(i,2)*incmat1(j,2)*abs(incmat1(i,1)-incmat1(j,1))
		ginicons=ginicons+consmat1(i,2)*consmat1(j,2)*abs(consmat1(i,1)-consmat1(j,1))
	end do
end do



giniinc=giniinc/(2.0*totinc)
ginicons=ginicons/(2.0*totcons)

end subroutine gini2