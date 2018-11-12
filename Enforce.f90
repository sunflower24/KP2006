
subroutine enforce

    use params


    
    implicit none
	integer ::maxiter=300,iter,inds(2),neq,neqi,i,hc,helpc,helpcc
	real(prec) :: cdif,adif,enfv(ns,ni),ladif,vdif,npla,vals(2)
	real(prec) ::vmn(ns,ni),adj
	real(prec), dimension(nt,ntypes,ns,ni,pa)::oldgrida,oldcfun,oldvfun,oldlafun
	real(prec), dimension(nt,ntypes,ns,ni,pa,ns,ni)::oldafun
	real(prec), dimension(ns,ni) ::valport


	real(prec),dimension(ns*ni)::sol1,gue1


	external sysbin
	external borrow

	open(unit=11,file='drules.txt')


adj=0.4

print*,'R is ',R(tc)

do iter=1,maxiter

!print*,iter

! Solve for the state-contingent borrowing constraints

select case(choice)

case(1)

amin(tc,:,:,:)=0.0

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			gue1(1)=grida(tc,tyc,sc,ic,1)
			fnorm=0.0
			call dneqnf(borrow,errel,1,itmax,gue1,sol1,fnorm)
			amin(tc,tyc,sc,ic)=sol1(1)
		end do
	end do
end do


case(2)
	if (tc==1) then
		do tyc=1,ntypes
			amin(tc,tyc,:,:)=0.0
			do sc=1,ns
				do ic=1,ni
					gue1(1)=0.0
					fnorm=0.0
					call dneqnf(borrow,errel,1,itmax,gue1,sol1,fnorm)
					amin(tc,tyc,sc,ic)=sol1(1)
				end do
			end do
		end do
	elseif (tc==nt) then
		amin(tc,:,:,:)=amin(1,:,:,:)
	endif

end select

! Redefine Grid

oldgrida(tc,:,:,:,:)=grida(tc,:,:,:,:)
oldcfun(tc,:,:,:,:)=cfun(tc,:,:,:,:)
oldafun(tc,:,:,:,:,:,:)=afun(tc,:,:,:,:,:,:)
oldvfun(tc,:,:,:,:)=vfun(tc,:,:,:,:)
oldlafun(tc,:,:,:,:)=lafun(tc,:,:,:,:)

call dgrid

! Redefine all functions to lie on the new grid

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				call basefun(oldgrida(tc,tyc,sc,ic,:),pa,grida(tc,tyc,sc,ic,ac),vals,inds)
				cfun(tc,tyc,sc,ic,ac)=vals(1)*oldcfun(tc,tyc,sc,ic,inds(1))+vals(2)*oldcfun(tc,tyc,sc,ic,inds(2))
				vfun(tc,tyc,sc,ic,ac)=vals(1)*oldvfun(tc,tyc,sc,ic,inds(1))+vals(2)*oldvfun(tc,tyc,sc,ic,inds(2))
				lafun(tc,tyc,sc,ic,ac)=vals(1)*oldlafun(tc,tyc,sc,ic,inds(1))+vals(2)*oldlafun(tc,tyc,sc,ic,inds(2))
				do scc=1,ns
					do icc=1,ni
						afun(tc,tyc,sc,ic,ac,scc,icc)=vals(1)*oldafun(tc,tyc,sc,ic,inds(1),scc,icc)+vals(2)*oldafun(tc,tyc,sc,ic,inds(2),scc,icc)
					end do
				end do
			end do
		end do
	end do
end do

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
		if ( amin(tc,tyc,sc,ic)  < grida(tc,tyc,sc,ic,1) -0.001 )  then
			print*,'Grid too small'
			pause
		endif
		end do
	end do
end do

!	print*,'Iteration number ',iter

cfunn(tc,1:ntypes,1:ns,1:ni,1:pa)=0.0
afunn(tc,1:ntypes,1:ns,1:ni,1:pa,1:ns,1:ni)=0.0

do tyc=1,ntypes
	do sc=1,ns		
		do ic=1,ni
       		do ac=1,pa

				do neq=ns*ni,1,-1
					
					
						! Determine which constraints are binding

						do scc=1,ns
							if ( ( neq>(scc-1)*ni) .and. (neq<=scc*ni) ) then

								hsc=scc
								hic=neq-(scc-1)*ni
								do hc=1,hsc-1
									gue1((hc-1)*ni+1:hc*ni)=afun(tc,tyc,sc,ic,ac,hc,:)
								end do
								gue1((hsc-1)*ni+1:neq)=afun(tc,tyc,sc,ic,ac,hsc,1:hic)
								fnorm=0.0
								call dneqnf(sysbin,errel,neq,itmax,gue1(1:neq),sol1(1:neq),fnorm)
								afunn(tc,tyc,sc,ic,ac,1:ns,1:ni)=amin(tc,tyc,1:ns,1:ni)
								do hc=1,hsc-1				
									afunn(tc,tyc,sc,ic,ac,hc,:)=sol1((hc-1)*ni+1:hc*ni)
								end do
								afunn(tc,tyc,sc,ic,ac,hsc,1:hic)=sol1((hsc-1)*ni+1:neq)
								
								valport=0.0
								do helpc=1,ns
									do helpcc=1,ni
										valport(helpc,helpcc)=(1.0/R(tc))*probs(tc,sc,helpc)*pi(tc,helpcc)*afunn(tc,tyc,sc,ic,ac,helpc,helpcc) 
									end do
								end do
								cfunn(tc,tyc,sc,ic,ac)= grida(tc,tyc,sc,ic,ac)+w(tc)*yat(tc,tyc,sc,ic)-sum(valport) 
								enfv(1:ns,1:ni)=afunn(tc,tyc,sc,ic,ac,1:ns,1:ni)-amin(tc,tyc,1:ns,1:ni)	
							end if
						end do		
			
						if  (minval(enfv)>-.000000001 ) then 
							goto 2001	
						end if
				
				end do

	! If program gets here all constraints are binding
	afunn(tc,tyc,sc,ic,ac,1:ns,1:ni)=amin(tc,tyc,1:ns,1:ni)
	valport=0.0
	do helpc=1,ns
		do helpcc=1,ni
			valport(helpc,helpcc)=(1.0/R(tc))*probs(tc,sc,helpc)*pi(tc,helpcc)*afunn(tc,tyc,sc,ic,ac,helpc,helpcc) 
		end do
	end do
	cfunn(tc,tyc,sc,ic,ac)=grida(tc,tyc,sc,ic,ac)+w(tc)*yat(tc,tyc,sc,ic)-sum(valport ) 

   
2001	   lafunn(tc,tyc,sc,ic,ac)=(1.0-bta)*MU(cfunn(tc,tyc,sc,ic,ac))

		   do scc=1,ns
				do icc=1,ni
					call basefun(grida(tc,tyc,scc,icc,:),pa,afunn(tc,tyc,sc,ic,ac,scc,icc),vals,inds)
					vmn(scc,icc)= bta*probs(tc,sc,scc)*pi(tc,icc)*(vals(1)*vfun(tc,tyc,scc,icc,inds(1))+vals(2)*vfun(tc,tyc,scc,icc,inds(2)))
				end do
		   end do

		   vfunn(tc,tyc,sc,ic,ac)=(1.0-bta)*U(cfunn(tc,tyc,sc,ic,ac))+sum(vmn) 
		   end do
       end do
	end do
end do

    cdif=maxval( abs( (cfunn(tc,1:ntypes,1:ns,1:ni,1:pa)-cfun(tc,1:ntypes,1:ns,1:ni,1:pa))/cfun(tc,1:ntypes,1:ns,1:ni,1:pa) ) )
 	adif=maxval( abs( (afunn(tc,1:ntypes,1:ns,1:ni,1:pa,1:ns,1:ni)-afun(tc,1:ntypes,1:ns,1:ni,1:pa,1:ns,1:ni)) ) )
	ladif=maxval( abs( (lafunn(tc,1:ntypes,1:ns,1:ni,1:pa)-lafun(tc,1:ntypes,1:ns,1:ni,1:pa)) ) )
	vdif=maxval( abs( (vfunn(tc,1:ntypes,1:ns,1:ni,1:pa)-vfun(tc,1:ntypes,1:ns,1:ni,1:pa)) ) )

	cfun(tc,1:ntypes,1:ns,1:ni,1:pa)=adj*cfunn(tc,1:ntypes,1:ns,1:ni,1:pa)+(1.0-adj)*cfun(tc,1:ntypes,1:ns,1:ni,1:pa)
	afun(tc,1:ntypes,1:ns,1:ni,1:pa,1:ns,1:ni)=adj*afunn(tc,1:ntypes,1:ns,1:ni,1:pa,1:ns,1:ni)+(1.0-adj)*afun(tc,1:ntypes,1:ns,1:ni,1:pa,1:ns,1:ni)
	lafun(tc,1:ntypes,1:ns,1:ni,1:pa)=adj*lafunn(tc,1:ntypes,1:ns,1:ni,1:pa)+(1.0-adj)*lafun(tc,1:ntypes,1:ns,1:ni,1:pa)
	vfun(tc,1:ntypes,1:ns,1:ni,1:pa)=adj*vfunn(tc,1:ntypes,1:ns,1:ni,1:pa)+(1.0-adj)*vfun(tc,1:ntypes,1:ns,1:ni,1:pa)

       
!	print*,cdif,adif,ladif,vdif,iter	  
    if ( max(cdif,adif,ladif,vdif) < 100.0*tol ) then
    exit
    endif
       
       
end do

print*,'convergence achieved after '
print*,iter
print*,'iterations'

end subroutine enforce

