 subroutine mctransnew
 use params
 use MSIMSL

 implicit none

 integer ini


real(prec),dimension(ns*pa*ntypes)::lorincx,lorincy,lorconsx,lorconsy
real(prec) :: giniinc,ginicons
real(prec):: quint(5),giniy
real(prec) :: stdlyper,stdlytra,stdlytot,meanp,meant,meaniid
real(prec) :: perminc(nt,ntypes)
real(prec) :: cperm(2),ginicm(2)
real(prec) :: test
real(prec) :: avgy(nt)
 
 
 complex(8),dimension(ns,ns):: Evec
 complex(8),dimension(ns) :: Eval
 real,dimension(ns,ns):: revec
 real,dimension(ns) :: reval



  ! modified Heaton and Lucas
 real(prec) :: rho
 real(prec) :: sigmaeps0  ! Initial value of the standard deviation of the epsilons
 real(prec) :: sigmaeps
 real(prec) :: toll=.00001



 real(prec),dimension(ns) :: esttax
 real(prec) acoeff,bcoeff,bhat
 integer,parameter :: npyc=24		 ! number of periods the income process is changing
 real(prec),dimension(nt,2) :: gdamod
 real(prec),dimension(nt,8) :: varinc
 real(prec),dimension(ntypes)::sumhelp


! Income Gini from CEX, household equivalent scales

! gdamod(1:npyc,1)=(/ 0.306, 0.297, 0.305, 0.317, 0.327, 0.332, 0.336, 0.342, 0.348, 0.358, 0.371, 0.383, 0.395, 0.399, 0.392, 0.394, 0.384, 0.384, 0.384, 0.381, 0.390, 0.384, 0.384, 0.393, 0.407, 0.416, 0.408/)

! Variance of Log Income from CEX, after taking out age and sex
! NOW WE NEED TO SPECIFY 4 COMPONENTS: permanent, persistent, transitory, and total as the sum of them

! Permanent Part
varinc(1:npyc,1)=(/0.0658,  0.0676,  0.0696,  0.0716,  0.0737,  0.0758,  0.0781,  0.0805,  0.0829,  0.0855,  0.0881,  0.0908,  0.0936,  0.0965,  0.0995,  0.1026,  0.1057,  0.1090,  0.1123,  0.1158,  0.1193,  0.1229,  0.1266,  0.1304 /)

! Persistent Part
varinc(1:npyc,2)=(/0.1928, 0.1969,0.2010, 0.2049, 0.2087, 0.2124, 0.2160, 0.2195, 0.2229, 0.2261, 0.2293, 0.2323, 0.2352, 0.2380, 0.2407, 0.2433, 0.2457, 0.2481, 0.2503, 0.2524, 0.2544, 0.2563, 0.2581,0.2598  /)

! Transitory Part
varinc(1:npyc,3)=(/0.0744,0.0769,0.0794,0.0818,0.0842,0.0865,0.0888,0.0910,0.0931,0.0952,0.0972,0.0991,0.1010,0.1028,0.1046,0.1063,0.1080,0.1095,0.1111,0.1125,0.1139,0.1153,0.1166,0.1178/)

! Total Variance

varinc(1:npyc,4)=(/0.3329,    0.3415,       0.3499,       0.3583,       0.3666,       0.3748,       0.3829,       0.3909,       0.3989,       0.4067,       0.4145,       0.4222,       0.4298,       0.4374,       0.4448,       0.4521,       0.4594,       0.4666,       0.4737,       0.4807,       0.4877,       0.4945,       0.5013,       0.5079/)

! income gini from CPS

! gdamod(1:npyc,1)=(/ 0.290, 0.295, 0.299, 0.305, 0.31, 0.315, 0.32, 0.325, 0.3267,0.3347,0.3455,0.3545,0.3579,0.3582,0.3613,0.3661,0.3704,0.3680,0.3665,0.3886,0.3912,0.3975,0.3953,0.4141,0.4162,0.4185/)

! gdamod(npyc+1:nt,1)=gdamod(npyc,1)

varinc(npyc+1:nt,1)=varinc(npyc,1) 
varinc(npyc+1:nt,2)=varinc(npyc,2) 
varinc(npyc+1:nt,3)=varinc(npyc,3) 
varinc(npyc+1:nt,4)=varinc(npyc,4) 

! Specify the iid part

! Specify the iid part

do tc=1,npyc
	do ic=1,ni
		pi(tc,ic)=1.0/ni			
		iidstates(tc,ic)=1.0 + ( 2.0*ic-3.0 ) * sqrt(varinc(tc,3)) 
!		iidstates(tc,ic)=1.0 + ( ic-2.0 ) * sqrt(1.5*varinc(tc,3))
	end do
end do

do tc=npyc+1,nt
	do ic=1,ni
		pi(tc,ic)=1.0/ni			
		iidstates(tc,ic)=iidstates(npyc,ic)
	end do
end do

iidstates=exp(iidstates)

! Specify the Persistent Part
     
nvar=1
nval(1)=ns
nsm=ns

allocate(mstae(ns**nvar))	
allocate(mstates(ns**nvar))
allocate(mprobs(ns**nvar,ns**nvar)) 


do tc=1,nt

	rho=0.9989

	sigmaeps=8.57* sqrt(varinc(tc,2))*sqrt(1.0-rho**2)
		
	theta(1)=0.0
	theta(2)=rho
	theta(3)=sigmaeps**2
	call markov
	probs(tc,1:ns,1:ns)=mprobs(1:ns,1:ns)

	if (tc==1) then
		diste(tc,1:ns)=mstae(1:ns)
		CALL DEVCRG (ns, mprobs, ns, EVAL, EVEC, ns)
		reval=dreal(eval)
		revec=dreal(evec)
		print*,'Eigen-values computed '
		print*,reval
			pause
	else
		diste(tc,1:ns)=matmul(diste(tc-1,1:ns),probs(tc-1,1:ns,1:ns))
	end if

	states(tc,1:ns)=exp(mstates(1:ns))

end do

! Specify the Variability of the permanent component and proportion of types

do tc=1,nt
	do tyc=1,ntypes
		perminc(tc,tyc)=exp(  (2*tyc-3)*sqrt(varinc(tc,1))   )
	!	perminc(tc,tyc)=exp(  (2*tyc-3)*sqrt(varinc(1,1))   )
		pitypes(tyc)=1.0/ntypes
	end do
end do

! Specify income 

do tc=1,nt
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				yat(tc,tyc,sc,ic)=perminc(tc,tyc)*states(tc,sc)*iidstates(tc,ic)
			end do
		end do
	end do
end do

! Normalize income so that total income sums to 1

do tc=1,nt
avgy(tc)=0.0
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				avgy(tc)=avgy(tc)+yat(tc,tyc,sc,ic)*pitypes(tyc)*pi(tc,ic)*diste(tc,sc)
			end do
		end do
	end do
end do

do tc=1,nt
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				yat(tc,tyc,sc,ic)=yat(tc,tyc,sc,ic)/avgy(tc)
			end do
		end do
	end do
end do


! Check whether total endowments sum to 1

do tc=1,nt
	test=0.0
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				test=test+yat(tc,tyc,sc,ic)*pitypes(tyc)*pi(tc,ic)*diste(tc,sc)
			end do
		end do
	end do
	if ( abs(test-1.0)>0.0001 ) then
		print*,'Total Endowments not equal to 1'
		pause
	endif
end do

! Compute Statistics for the Income Process

do tc=1,nt

	! Permanent Part

	meanp=sum( pitypes(1:ntypes)*log(perminc(tc,1:ntypes)) )
	varvec(tc,1)=sum( pitypes(1:ntypes)* ( log(perminc(tc,1:ntypes))-meanp)**2 )
	stdsvec(tc,1)=sqrt(varvec(tc,1))

	! Persistent Part
	
	meant=sum( diste(tc,1:ns)*log(states(tc,1:ns)) )
	varvec(tc,2)=sum( diste(tc,1:ns)* ( log(states(tc,1:ns))-meant)**2 )	
	stdsvec(tc,2)=sqrt(varvec(tc,2))
	
	! Transitory Part

	meaniid=sum(pi(tc,1:ni)*log(iidstates(tc,1:ni)))
	varvec(tc,3)=sum( pi(tc,1:ni)* ( log(iidstates(tc,1:ni))-meaniid)**2 )
	stdsvec(tc,3)=sqrt(varvec(tc,3))

	varvec(tc,4)=varvec(tc,1)+varvec(tc,2)+varvec(tc,3)
	stdsvec(tc,4)=sqrt(varvec(tc,4))


end do

! Compute Income Gini

do tc=1,nt
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				omega(1,tyc,sc,ic,1)=pitypes(tyc)*diste(tc,sc)*pi(tc,ic)
				iniomega=omega(1,1:ntypes,1:ns,1:ni,1:pa)
			end do
		end do
	end do
!	call gini2(iniomega,lorincx,lorincy,lorconsx,lorconsy,giniinc,ginicons)		
!	giniyvec(tc)=giniinc
	omega=0.0
	iniomega=0.0
end do

! Construct Consumption under Perfect Insurance within Types

cperm=0.0
do tyc=1,ntypes
	cperm(tyc)=perminc(npyc+1,tyc)*(bta**(npyc+1.0))
	do tc=1,npyc
		cperm(tyc)=cperm(tyc)+(1.0-bta)*perminc(tc,tyc)*(bta**(tc-1))
	end do
end do

! Compute Complete Markets Gini

ginicm(1)=0.5*(1.0-perminc(1,1)/perminc(1,2))/(1.0+perminc(1,1)/perminc(1,2))
ginicm(2)=0.5*(1.0-cperm(1)/cperm(2))/(1.0+cperm(1)/cperm(2))

gdamod(1:nt,2)=giniyvec(1:nt)  ! look at this matrix to see how the model gini compares to the data gini

varinc(1:nt,5:8)=varvec(1:nt,1:4)


end subroutine mctransnew
