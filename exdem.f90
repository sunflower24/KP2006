subroutine exdem

use params
implicit none

integer::irc
integer,parameter::nr=20
real(prec)::minr,maxr
real(prec),dimension(nr)::gridr,exd,dem,sup

minr=1.02
maxr=1.0/bta-0.0001

do irc=1,nr
	gridr(irc)=minr+(maxr-minr)*((irc-1.0)/(nr-1.0))
end do

do irc=1,nr
	call residss(gridr(irc),exd(irc))
	dem(irc)=assdemand(tc)
	sup(irc)=R(tc)*KA(tc)
end do

open(unit=17,file='exdem.txt')
do irc=1,nr
	write(17,'(4f13.6)') gridr(irc),dem(irc),sup(irc),exd(irc)
end do	
close(17)

end subroutine exdem

