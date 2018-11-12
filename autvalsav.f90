subroutine autvalsav

use params
implicit none

integer ::maxiter=2500,iter,inds(2),neq,i
real(prec) :: cdif,sdif,enfv(ns),vdif,npla,vals(2),vpdif
real(prec) ::vmn(ns,ni)
real(prec),dimension(ns)::sol1,gue1
real(prec) :: gue(1),sol(1)
real(prec) :: debtcon,adj=1.0,frac=1.0

	external sysex

! Choose at what rate to save in autarky

select case(choice2)

case(2)

Rm(tc)=Rex

case(3)

Rm(tc)=R(tc)

end select

! Guess savings functions

call guess2

! Compute optimal Savings

do iter=1,maxiter
!	print*,'Iteration number ',iter
	sfunn(tc,1:ntypes,1:ns,1:ni,1:pm)=0.0
	cmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)=0.0
	vpmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)=0.0
	vmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)=0.0
	
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni		
       			do rc=1,pm

					gue=sfun(tc,tyc,sc,ic,rc) !gridm(1)
       				call dneqnf(sysex,errel,1,itmax,gue,sol,fnorm)

					if (sol(1) < 0.0 ) then
						sol(1) = 0.0
					end if
					sfunn(tc,tyc,sc,ic,rc)=sol(1)
					cmfunn(tc,tyc,sc,ic,rc)=gridm(tyc,rc)+w(tc)*yat(tc,tyc,sc,ic)-sfunn(tc,tyc,sc,ic,rc)/Rm(tc)
					vpmfunn(tc,tyc,sc,ic,rc)=(1.0-bta)*MU(cmfunn(tc,tyc,sc,ic,rc))
					vmn=0.0
					do scc=1,ns
						do icc=1,ni
							call basefun (gridm(tyc,:),pm,sfunn(tc,tyc,sc,ic,rc),vals,inds)
							vmn(scc,icc)=bta*probs(tc,sc,scc)*pi(tc,icc)*( vals(1)*vmfun(tc,tyc,scc,icc,inds(1))+vals(2)*vmfun(tc,tyc,scc,icc,inds(2)) )
						end do
					end do
					vmfunn(tc,tyc,sc,ic,rc)=(1.0-bta)*U(cmfunn(tc,tyc,sc,ic,rc))+sum(vmn) 
				end do
			end do
		end do
	end do
       
    sdif=100*maxval( abs( (sfunn(tc,1:ntypes,1:ns,1:ni,1:pm)-sfun(tc,1:ntypes,1:ns,1:ni,1:pm)) ) )
 	cdif=100*maxval( abs( (cmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)-cmfun(tc,1:ntypes,1:ns,1:ni,1:pm)) ) )
	vpdif=maxval( abs( (vpmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)-vpmfun(tc,1:ntypes,1:ns,1:ni,1:pm)) ) )
	vdif=maxval( abs( (vmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)-vmfun(tc,1:ntypes,1:ns,1:ni,1:pm)) ) )


	cmfun(tc,1:ntypes,1:ns,1:ni,1:pm)=(1.0-adj)*cmfun(tc,1:ntypes,1:ns,1:ni,1:pm)+adj*cmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)
	sfun(tc,1:ntypes,1:ns,1:ni,1:pm)=(1.0-adj)*sfun(tc,1:ntypes,1:ns,1:ni,1:pm)+adj*sfunn(tc,1:ntypes,1:ns,1:ni,1:pm)
	vpmfun(tc,1:ntypes,1:ns,1:ni,1:pm)=(1.0-adj)*vpmfun(tc,1:ntypes,1:ns,1:ni,1:pm)+adj*vpmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)
	vmfun(tc,1:ntypes,1:ns,1:ni,1:pm)=(1.0-adj)*vmfun(tc,1:ntypes,1:ns,1:ni,1:pm)+adj*vmfunn(tc,1:ntypes,1:ns,1:ni,1:pm)
       
!	print*,max(sdif,cdif,vpdif,vdif)
    
    
    if ( max(sdif,cdif,vpdif,vdif) < tol ) then
    exit
    endif
       
       
end do

print*,'convergence achieved after '
print*,iter
print*,'iterations'

vmaut(tc,1:ntypes,1:ns,1:ni)=vmfun(tc,1:ntypes,1:ns,1:ni,1)

! Test: how are things changing if we take as punishment eating a fraction of your endowment in the period of 
! default, and saving from then on

!do tyc=1,ntypes
!	do sc=1,ns
!		do ic=1,ni
!			do scc=1,ns
!				do icc=1,ni
!					vmn(scc,icc)=bta*probs(tc,sc,scc)*pi(tc,icc)*vmfun(tc,tyc,scc,icc,1)
!				end do
!			end do
!			vmaut(tc,tyc,sc,ic)=(1.0-bta)*U(frac*yat(tc,tyc,sc,ic))+sum(vmn)
!		end do
!	end do
!end do

end subroutine autvalsav