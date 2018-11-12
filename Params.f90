		  	  	    
 module params
 implicit none

 integer,parameter :: prec=selected_real_kind(9,300)

 integer,parameter::ns=7	     ! Number of individual states for the persistent part
 integer,parameter::ni=2		 ! Number of individual states for the iid part
 integer :: NVAR,nlag=1
 integer ::nval(10)
 integer :: nsm
 integer,parameter::nt=100 			! Number of time periods for transtion to be completed
 integer::tc,ttc						! Time counter
						

 ! Parameters for type definition
 integer,parameter :: ntypes=2 
 integer :: tyc   
 real(prec),dimension(ntypes) :: typevec
 real(prec),dimension(nt) :: alfavec
 real(prec)::alfa
 real(prec),dimension(ntypes) :: pitypes

 real(prec),parameter :: sigma=1.0  ! Risk aversion
 real(prec),parameter :: bta=0.95945   ! subjective discount factor
 real(prec),parameter :: betainv=1.0/bta  	! subjective discount factor
 real(prec),parameter :: minc=0.0001, minu=0.0001  ! mimimum cons and utility
 real(prec),parameter :: Rex=1.04

 real(prec),parameter :: rat=2.6	! Capital-Output Ratio (Wealth-Output-Ratio)
 real(prec),parameter :: rtarget=0.04 ! Target for the Interest Rate
 real(prec)::A						! Constant in the Production Function
 real(prec),parameter :: alpha=0.3	! Capital Share
 real(prec):: delta					! Depreciation rate


 real(prec) :: theta(100)		   ! 1:nvar = constant
								   ! nvar*nvar = autoregressive matrix
								   ! nvar*nvar = epsilons var cov matrix


 ! Parameters derived from Tauchen procedure (Markov)
 real(prec),allocatable :: mstae(:)			 ! Stationary distribution of e
 real(prec),allocatable :: mstates(:)		 ! Matrix of states e1,e2
 real(prec),allocatable :: mprobs(:,:)      ! Transition probabilities


 ! Parameters used in the main program (in the VAR case they are tranformed using the sty routine)
 
 real(prec) :: states(nt,ns)	 ! Matrix of states for the persistent process
 real(prec) :: probs(nt,ns,ns)       ! Transition probabilities for the persistent process
 real(prec) :: diste(nt,ns)			 ! Distribution of persistent part 
 real(prec) :: giniyvec(nt),ginicvec(nt)           ! vector of income gini coefficients 

 ! Parameters governing the iid shocks

 real(prec) :: pi(nt,ni)			! Distribution of iid shocks
 real(prec) :: iidstates(nt,ni)		! Transitory States


 integer,parameter :: pa=90
 real(prec),dimension(nt,ntypes,ns,ni,pa) :: grida

 real(prec) :: R(nt),Rnew(nt)					! Intertemporal price
 real(prec) :: w(nt)							! Wages
 real(prec) :: KA(nt)							! Capital Stock
 real(prec) :: LS(nt)							! Labor Supply
 real(prec) :: yat(nt,ntypes,ns,ni)				! after tax income
 
 integer, parameter	:: pc=201			! grid points for c
 real(prec),dimension(nt,ntypes,pa*ns*ni) :: distr ! Distribution over w and y


 real(prec),dimension(nt) :: govint2,privint2,totint2,avglc,stdlc,stdly,stdle,varlc
 real(prec),dimension(nt,4) :: stdsvec,varvec
 real(prec),dimension(nt,ntypes)::tavglc,tvarlc,tstdlc
 real(prec),dimension(nt)::bvarlc,bstdlc

 
 real(prec)::vcm(ntypes)

! Parameters for computing autarky with savings
 integer,parameter :: pm=90
 real(prec),dimension(ntypes,pm) :: gridm
 real(prec),dimension(nt,ntypes,ns,ni,pm)::cmfun,cmfunn,sfun,sfunn,vpmfun,vpmfunn,vmfun,vmfunn
 integer rc,fsav
 real(prec) :: Rm(nt),minvaut(nt),maxvaut(nt)				 ! Intertemporal price

 integer sc,scc,ac,ic,icc,hsc,hic

 
 real(prec) :: errel=0.00000001, fnorm	! parameters used by the non lin eq solver
 integer :: itmax=50000, choice, choice2
 real(prec) :: tol=0.0001            ! tolerance parameter


 real(prec),dimension(nt,ntypes,ns,ni,pa) ::lafun,lafunn,omega,vfun,vfunn,cfun,cfunn
 real(prec),dimension(nt,ntypes,ns,ni,pa,ns,ni) ::afun,afunn
 real(prec),dimension(ntypes,ns,ni,pa)::iniomega,finomega
 
 real(prec),dimension(nt,ntypes,ns,ni)::vmaut,expinc,amin
 real(prec),dimension(nt,ntypes,pc,ns,ni)::cdist
 real(prec),dimension(pc)::gridc

 real(prec),dimension(nt)::tote,cons,cons1,resd,resd1,assdemand,assdemand1,exput1,exput2,ginivec,totint,totint1
 real(prec),dimension(nt,ntypes,ns,ni)::help



 integer porl

 contains

 

function U(c)
implicit none
real(prec) U,c																		                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
 
if (abs(sigma-1.0)<0.0001) then
	if (c<minc) then
		U=log(minc)+(c-minc)/minc	
	else
		U=log(c)
	endif
else
	if (c<minc) then
		U=(minc**(1.0-sigma))/(1.0-sigma) + (c-minc)*minc**(-sigma)	
	else
		U=(c**(1.0-sigma))/(1.0-sigma)
	endif
endif
  
end function U


function MU(c)
implicit none
real(prec) MU,c

if (abs(sigma-1.0)<0.0001) then
	if (c<minc) then
		MU=1.0/minc - (c-minc)/(minc**2.0)
	else
		MU=1.0/c
	endif
else
	if (c<minc) then
		MU=minc**(-sigma) - (c-minc)*sigma*minc**(-1.0-sigma)	
	else
		MU=c**(-sigma)
	endif
endif
   
end function MU


 function C(u)
 implicit none
 real(prec) C,u

 
 

 if ( (u*(1.0-sigma)) <= 0.0)  then  
	u=minu*(1.0-sigma)			   
 end if

 if (abs(sigma-1.0)<.0001) then
 C = exp(u)
 else
 C = ( (1.0-sigma)*u  )** ( 1.0/(1.0-sigma) )
 endif
  
 end function C															 


 function MC(u)
 implicit none
 real(prec) MC,u

 if (abs(sigma-1.0)<.0001) then
 MC = exp(u)
 else
 MC = C(u)**sigma
 endif

 end function MC

 function F(cap,lab)
 implicit none
 real(prec) F,cap,lab
 
 F=A*(cap**alpha)*lab**(1.0-alpha)
 
 end function F
 
 function MPL(cap,lab)
 implicit none
 real(prec) MPL,cap,lab
 
 MPL=(1.0-alpha)*A*(cap**alpha)*lab**(-alpha)						

 end function MPL

 function MPC(cap,lab)
 implicit none
 real(prec) MPC,cap,lab
 
 MPC=alpha*A*(cap**(alpha-1.0))*lab**(1.0-alpha) - delta					

 end function MPC								



	subroutine basefun (grid_x,npx,x,vals,inds) 
	implicit none
	! this subroutine returns the values and the indices of the two basis
	! functions that are positive on a given x in the grid_x


	real(prec),intent(in) :: x
	integer , intent(in):: npx
	real(prec), intent(in) :: grid_x (npx)
	real(prec), intent(out) ::vals(2)
	integer ,intent(out) ::inds(2)
	integer :: i,ju,jl,jm

	


	jl=1     ! 
	ju=npx   !	


	do

	if (ju-jl<=1) exit
		
	jm=(ju+jl)/2
	if (x>=grid_x(jm)) then

		jl=jm
	else
		ju=jm
	endif

	end do


	i=jl+1



		vals(2)=( x-grid_x(i-1) )/(grid_x(i)-grid_x(i-1))
		vals(1)=( grid_x(i)-x )/(grid_x(i)-grid_x(i-1))
		inds(2)=i
		inds(1)=i-1

	

	end subroutine basefun



 end module params