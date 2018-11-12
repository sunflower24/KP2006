subroutine diss

use params

implicit none


integer inds(2),i
integer unevalc,uneval(ni*ns*pa)
real(prec), dimension(ni*ns*pa,ni*ns*pa)::Amat,revec 
real(prec) :: vals(2),sdtol=0.0001,reval(ni*ns*pa),test,test1
complex(8),dimension(ni*ns*pa,ni*ns*pa):: Evec
complex(8),dimension(ni*ns*pa) :: Eval

	open(unit=13,file='diss.txt')
	open(unit=15,file='amat.txt')



do tyc=1,ntypes

Amat=0.0

do sc=1,ns
	do ic=1,ni
		do ac=1,pa
			do scc=1,ns
				do icc=1,ni
					call basefun(grida(tc,tyc,scc,icc,:),pa,afun(tc,tyc,sc,ic,ac,scc,icc),vals,inds)

!		if (vals(1)<0.0) then
!		vals(1)=0.0
!		endif

!		if (vals(2)>1.0) then
!		vals(2)=1.0
!		endif

!		if (vals(2)<0.0) then
!		vals(2)=0.0
!		endif

!		if (vals(1)>1.0) then
!		vals(1)=1.0
!		endif

					AMat(((scc-1)*ni+icc-1)*pa+inds(1),((sc-1)*ni+ic-1)*pa+ac)=vals(1)*probs(tc,sc,scc)*pi(tc,icc)
					AMat(((scc-1)*ni+icc-1)*pa+inds(2),((sc-1)*ni+ic-1)*pa+ac)=vals(2)*probs(tc,sc,scc)*pi(tc,icc)

				end do
			end do
		end do
	end do
end do

! Test the Matrix

do i=1,ns*ni*pa
	test=sum(AMAT(:,i))
	test1=sum(AMAT(i,:))
	if ( ( test> 1.0001 ) .or. (test < 0.999) ) then
		print*,'Something is wrong with transition matrix ',i,test
		pause
	end if	
!	if ( ( test1> 1.0001 ) .or. (test1 < 0.999) ) then
!		print*,'Something is wrong with transition matrix ',i,test1
!		pause
!	end if	
end do


CALL DEVCRG (ni*ns*pa, Amat, ni*ns*pa, EVAL, EVEC, ni*ns*pa)

reval=dreal(eval)
revec=dreal(evec)

print*,'Eigen-values computed '




if (reval(1)<1.0-sdtol) then
 print*,' No stationary distribution here !!'
end if



unevalc=0
do i=1,pa*ns*ni

	if (reval(i)<1.0-sdtol)  then
		exit
	end if 

	
	if (reval(i)>1.0+sdtol) then
		cycle
	end if

	unevalc=unevalc+1
	uneval(unevalc)=i
	print*,eval(i)


  end do						

if ( unevalc==1) then 
	print*,'Unique stationary distribution '

	distr(tc,tyc,1:pa*ns*ni)=pitypes(tyc)*revec(1:pa*ns*ni,uneval(1))/sum(revec(1:pa*ns*ni,uneval(1)))

	rewind(13)

	do sc=1,ns
		do ic=1,ni
			write(13,fmt=*) distr(tc,tyc,(((sc-1)*ni+ic-1))*pa+1:((sc-1)*ni+ic)*pa)
		end do
	end do 
	write(13,fmt=*)	grida

	print*,' Distr. saved '
	rewind(13)


end if

if (unevalc>1) then
	print*,'Warning...'
	print*,'No unique stationary distribution..'
	print*,'Number of unit Eigen-values ',unevalc
	pause
end if

if ( ( sum(distr(tc,tyc,1:pa*ns*ni)) < pitypes(tyc) - 0.0001 ) .or. (sum(distr(tc,tyc,1:pa*ns*ni)) > pitypes(tyc) + 0.0001 ) ) then
	print*,'Distribution does not sum to 1'
	print*,sum(distr(tc,tyc,1:pa*ns*ni)),tyc
	pause
end if

end do

end subroutine diss