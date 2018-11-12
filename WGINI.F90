subroutine cginiw(vector,dim,weights,gindex,quint)

use dflib
use params
implicit none


integer i,dim
real(prec),dimension(dim) :: vector,weights,rweights
real(prec),dimension(dim) :: lorx,lory
real(prec) gindex,mass,tcon,quint(5)


! first sort vec

CALL sortqq (loc(vector),dim,srt$real8) 

rweights=weights/sum(weights)
tcon=sum(vector*rweights)

! compute the lorenz curve

do i=1,dim
	mass=sum(vector(1:i)*rweights(1:i))
	lory(i)=mass/tcon
	lorx(i)=sum(rweights(1:i))/sum(rweights(1:dim))
	if(i>1) then

		if ((lorx(i-1)<.2).and.(lorx(i)>.2)) then
			quint(1)=lory(i)
		endif

		if ((lorx(i-1)<.4).and.(lorx(i)>.4)) then
			quint(2)=lory(i)-quint(1)
		endif

		if ((lorx(i-1)<.6).and.(lorx(i)>.6)) then
			quint(3)=lory(i)-quint(2)-quint(1)
		endif

		if ((lorx(i-1)<.8).and.(lorx(i)>.8)) then
			quint(4)=lory(i)-quint(3)-quint(2)-quint(1)
		endif
	endif

end do

quint(5)=1.0-quint(4)-quint(3)-quint(2)-quint(1)



gindex=sum(lorx-lory)/sum(lorx)

end subroutine cginiw