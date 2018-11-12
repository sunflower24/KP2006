subroutine dgrid0  ! Set the grid for Alvarez-Jermann economy
use params

implicit none
real(prec) :: alb,ustepsize,lstepsize,center,wwc,ascale
real(prec)::gue1(1),sol1(1),ahb
real(prec),parameter::curv=2.0

integer i

external valuecm

do tyc=1,ntypes
amin(tc,tyc,:,:)=0.0
do sc=1,ns
	do ic=1,ni
		gue1(1)=0.0
		call dneqnf(valuecm,errel,1,itmax,gue1,sol1,fnorm)
		amin(tc,tyc,sc,ic)=sol1(1)
		alb=sol1(1)
		ahb=10*w(tc)*yat(tc,tyc,ns,ni)
		ascale=ahb-alb
		grida(tc,tyc,sc,ic,1)=alb
		do ac=1,pa
			grida(tc,tyc,sc,ic,ac)=alb+ascale*((ac-1.0)/(pa-1.0))**curv
		end do
	end do
end do
end do

end subroutine dgrid0
