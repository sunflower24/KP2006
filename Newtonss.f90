subroutine Newtonss(fun,gues)

use params

implicit none

integer :: i, maxit=1000
real(prec) ::eps=0.0001,del, gues,guesp,ngues,nngues,fval1,fval2


do i=1,maxit


	ngues=gues
	call fun(gues,fval2)

!	pause

	if (abs(fval2)<50.0*tol) then
		exit
	end if

	guesp=gues+gues*eps

	call fun(guesp,fval1)

	del=( fval1-fval2 )/( ngues*eps )
	if (abs(del)<0.0000005) then
		del=0.0000005*abs(del)/del
	end if
	print*,'del is',del

	print*,'Old Guess is ',ngues
	
	nngues=ngues-fval2/del


if (nngues>1.0/bta) then
	nngues=1.0/bta-0.000001
end if

gues=0.6*ngues+0.4*nngues

print*,'Nguess is ',gues

end do


end subroutine Newtonss