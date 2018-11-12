subroutine dgrid2  ! Set the grid for the enforcement and exogenous incomplete mkts
use params

implicit none
real(prec) :: clb,chb,cstepsize,blb,bhb,bstepsize
real(prec) :: ascale,curv
integer i 
 
 
 
 ! set grid for consumption

clb=0.0001
chb=1.5*maxval(yat(1:nt,1:ntypes,ns,ni))
cstepsize=(chb-clb)/(pc-1)
gridc=(/(i*cstepsize,i=0,pc-1)/)
gridc=gridc+clb

 ! set grid for savings in autarky (same spacing as uncontingent bonds grid) 


do tyc=1,ntypes

ascale=10.0*yat(1,tyc,ns,ni)
curv=1.5

gridm(tyc,1)=0.0
do rc=2,pm
	gridm(tyc,rc)=ascale*((rc-1.0)/(pm-1.0))**curv
end do

end do 

end subroutine dgrid2