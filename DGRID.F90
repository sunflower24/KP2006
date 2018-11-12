subroutine dgrid  ! Set the grid for Alvarez-Jermann economy
use params

implicit none
real(prec) :: alb,ustepsize,lstepsize,center,wwc,ascale
real(prec)::gue1(1),sol1(1),ahb
real(prec),parameter::curv=2.0

integer i

external borrow

do tyc=1,ntypes

	if ( (tc>1) .and. (tc<nt-1) ) then

		do sc=1,ns
			do ic=1,ni
				alb=amin(tc,tyc,sc,ic)
				ahb=40*w(tc)*yat(tc,tyc,ns,ni) 
				ascale=ahb-alb
				grida(tc+1,tyc,sc,ic,1)=alb
				do ac=1,pa
					grida(tc+1,tyc,sc,ic,ac)=alb+ascale*((ac-1.0)/(pa-1.0))**curv
				end do
			end do
		end do

	else

		do sc=1,ns
			do ic=1,ni
				alb=amin(tc,tyc,sc,ic)
				ahb=40.0*w(tc)*yat(tc,tyc,ns,ni)   
				ascale=ahb-alb
				grida(tc,tyc,sc,ic,1)=alb
				do ac=1,pa
					grida(tc,tyc,sc,ic,ac)=alb+ascale*((ac-1.0)/(pa-1.0))**curv
				end do
			end do
		end do

	endif

end do

end subroutine dgrid



  