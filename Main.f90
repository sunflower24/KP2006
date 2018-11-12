program main

use params
use MSIMSL
implicit none

real(prec),dimension(nt,2):: helps
real(prec),dimension(2,nt):: helpst
real(prec):: concm

external residss

call erset(4,-1,0)          !Change stop-setting for fatal errors
							!if the 3rd element =0 then do not stop 
							!if the 3rd element =1 then stop

call erset(0,0,-1)

! Set Constant in the Production Function

A=(rat**(-alpha))*(1.0-alpha)**(alpha-1.0)
delta=alpha/rat-rtarget

! Initialize Equilibrium Factor Prices

do tc=1,nt
	KA(tc)=rat/(1.0-alpha)
	LS(tc)=1.0
	R(tc)=1.0 + MPC(KA(tc),LS(tc))
	w(tc)=MPL(KA(tc),LS(tc))
end do

! parametrize markov chains
call mctransnew

! design grids for ex incomplete, consumption and autarky with savings

call dgrid2 


 ! Choice of Debt Constraint

 print *, 'Which debt constraint do you want?'
 print *, '1: Adjusting to Changes in Income Process'
 print *, '2: No Adjustment of Constraint'

 read *, choice
 print *, ' '

! Choice of saving option in autarky
 
 print *,'Which Saving Option do you want?'
 print *,'1: No Saving in Autarky'
 print *,'2: Saving at exogenous Rate ',Rex
 print *,'3: Saving at equilibrium risk free rate'

 read *,choice2
 print *,' '

  
! Compute values of autarky for the definition of the grid

select case(choice)

case(1)
tc=1
call autval

tc=nt
call autval


do ttc=2,nt-1
	tc=nt+1-ttc
	call autvaltr          
	minvaut(tc)=vmaut(tc,1,1,1) 
	maxvaut(tc)=vmaut(tc,ntypes,ns,ni)
end do

case(2)


tc=1
call autval
tc=nt
call autval
do ttc=2,nt-1
	tc=nt+1-ttc
	call autvaltr          
	minvaut(tc)=vmaut(tc,1,1,1) 
	maxvaut(tc)=vmaut(tc,ntypes,ns,ni)
end do

end select

! Compute Complete Markets Ginis

call cm

tc=1
vcm=0.0
do tyc=1,ntypes
	concm=0.0
	do sc=1,ns
		do ic=1,ni
			concm=concm+diste(tc,sc)*pi(tc,ic)*yat(tc,tyc,sc,ic)*(w(tc)+(R(tc)-1.0)*KA(tc))
		end do
	end do
	vcm(tyc)=U(concm)		 
	print*,'Value of Complete Markets ',vcm(tyc)
	print*,'Value of Autarky ',vmaut(tc,tyc,1:ns,1:ni)
end do

tc=nt
vcm=0.0
do tyc=1,ntypes
	concm=0.0
	do sc=1,ns
		do ic=1,ni
			concm=concm+diste(tc,sc)*pi(tc,ic)*yat(tc,tyc,sc,ic)*(w(tc)+(R(tc)-1.0)*KA(tc))
		end do
	end do
	vcm(tyc)=U(concm)		 
	print*,'Value of Complete Markets ',vcm(tyc)
	print*,'Value of Autarky ',vmaut(tc,tyc,1:ns,1:ni)
end do
!pause


!call exdem

! Compute Initial Steady State
print*,'Initial Steady State Started'
call IniSS
print*,R(1)
pause

! Compute Final Steady State
print*,'Final Steady State Started'
call FinSS
print*,R(nt)
!pause

! Compute the Transition

call transition

! Save our results
call saving

end program main