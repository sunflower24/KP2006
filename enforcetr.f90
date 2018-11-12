
subroutine enforcetr

    use params


    
    implicit none
	integer ::maxiter=1000,iter,inds(2),neq,neqi,i,hc,helpc,helpcc
	real(prec) :: cdif,adif,enfv(ns,ni),ladif,vdif,npla,vals(2)
	real(prec) ::vmn(ns,ni),valport(ns,ni)



	real(prec),dimension(ns*ni)::sol1,gue1


	external sysbintr
	external borrowtr

	open(unit=12,file='drulestr.txt')


	select case(choice)

	case(1)

		do tyc=1,ntypes
			do sc=1,ns
				do ic=1,ni
					gue1(1)=0.0
					fnorm=0.0
					call dneqnf(borrowtr,errel,1,itmax,gue1,sol1,fnorm)
					help(tc,tyc,sc,ic)=sol1(1)
				end do
			end do
		end do

	case(2)

		do tyc=1,ntypes
			do sc=1,ns
				do ic=1,ni
					help(tc,tyc,sc,ic)=amin(1,tyc,sc,ic)
				end do
			end do
		end do

	end select


!	print*,'Iteration number ',iter
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
									gue1((hc-1)*ni+1:hc*ni)=afun(tc+1,tyc,sc,ic,ac,hc,:)
								end do
								gue1((hsc-1)*ni+1:neq)=afun(tc+1,tyc,sc,ic,ac,hsc,1:hic)
								fnorm=0.0
								call dneqnf(sysbintr,errel,neq,itmax,gue1(1:neq),sol1(1:neq),fnorm)
								afun(tc,tyc,sc,ic,ac,1:ns,1:ni)=amin(tc,tyc,1:ns,1:ni)
								do hc=1,hsc-1				
									afun(tc,tyc,sc,ic,ac,hc,:)=sol1((hc-1)*ni+1:hc*ni)
								end do
								afun(tc,tyc,sc,ic,ac,hsc,1:hic)=sol1((hsc-1)*ni+1:neq)

								valport=0.0
								do helpc=1,ns
									do helpcc=1,ni
										valport(helpc,helpcc)=(1.0/R(tc))*probs(tc,sc,helpc)*pi(tc+1,helpcc)*afun(tc,tyc,sc,ic,ac,helpc,helpcc) 
									end do
								end do
								cfun(tc,tyc,sc,ic,ac)= grida(tc,tyc,sc,ic,ac)+w(tc)*yat(tc,tyc,sc,ic)-sum(valport) 
								enfv(1:ns,1:ni)=afun(tc,tyc,sc,ic,ac,1:ns,1:ni)-amin(tc,tyc,1:ns,1:ni)	
							end if
						end do	

					if  (minval(enfv)>-.000000001 ) then 
						goto 2001	
					end if
				end do

	! If program gets here all constraints are binding
	afun(tc,tyc,sc,ic,ac,1:ns,1:ni)=amin(tc,tyc,1:ns,1:ni)
	valport=0.0
	do helpc=1,ns
		do helpcc=1,ni
			valport(helpc,helpcc)=(1.0/R(tc))*probs(tc,sc,helpc)*pi(tc+1,helpcc)*afun(tc,tyc,sc,ic,ac,helpc,helpcc) 
		end do
	end do
	cfun(tc,tyc,sc,ic,ac)=grida(tc,tyc,sc,ic,ac)+w(tc)*yat(tc,tyc,sc,ic)-sum(valport) 

	   
	if ( cfun(tc,tyc,sc,ic,ac) < 0.0 ) then
		print*,'Negative Consumption'
		print*,cfun(tc,tyc,sc,ic,ac),tc,tyc,sc,ic,ac
!		pause
	endif

2001	   lafun(tc,tyc,sc,ic,ac)=(1.0-bta)*MU(cfun(tc,tyc,sc,ic,ac))

		   do scc=1,ns
			do icc=1,ni
				call basefun (grida(tc+1,tyc,scc,icc,:),pa,afun(tc,tyc,sc,ic,ac,scc,icc),vals,inds)
				vmn(scc,icc)=bta*probs(tc,sc,scc)*pi(tc+1,icc)* (vals(1)*vfun(tc+1,tyc,scc,icc,inds(1))+vals(2)*vfun(tc+1,tyc,scc,icc,inds(2)))
			end do
		   end do

		   vfun(tc,tyc,sc,ic,ac)=(1.0-bta)*U(cfun(tc,tyc,sc,ic,ac))+sum(vmn ) 

		end do
	  end do
    end do
end do


end subroutine enforcetr

