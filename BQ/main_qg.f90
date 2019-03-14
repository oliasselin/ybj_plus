PROGRAM main

  USE parameters
  USE mpi
  USE fft
  USE init
  USE derivatives
  USE elliptic
  USE diagnostics
  USE files

  !********************** Declaring variables *****************************!

  double precision, dimension(n1d,n2d,n3h2)   :: ur,vr,wr,br                     !Velocity and potential temperature fields (r-space)
  double complex,   dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk                     !Velocity and potential temperature fields (k-space)

  double precision, dimension(n1d,n2d,n3h1)   :: war       !1st order vertical velocity as computed in QG (w^1)
  double complex,   dimension(iktx,ikty,n3h1) :: wak

  double precision, dimension(n1d,n2d,n3h1)   :: qr                  
  double complex,   dimension(iktx,ikty,n3h1) :: qk        
  double complex,   dimension(iktx,ikty,n3h1) :: qok          !Old q 
  double complex,   dimension(iktx,ikty,n3h1) :: qtempk        

  double complex,   dimension(iktx,ikty,n3h0) :: dqk         !dissipation
  double complex,   dimension(iktx,ikty,n3h1) :: psik                     !pressure, and rhs of pressure equation!
  double precision, dimension(n1d,n2d,n3h1)   :: psir
  double complex,   dimension(iktx,ikty,n3h1) :: psi_old      !For computing w...

  double complex,   dimension(iktx,ikty,n3h0) :: rhs         !RHS of elliptic equation (n3h0 version of q at n+1)
  double precision, dimension(n1d,n2d,n3h0)   :: rhsr

  double complex, dimension(iktx,n3, iktyp) :: qt             !Transposed (ky-parallelization) right-hand side   

  double precision, dimension(n1d,n2d,n3h0)   :: nqr                  
  double complex,   dimension(iktx,ikty,n3h0) :: nqk        

  equivalence(ur,uk)
  equivalence(vr,vk)
  equivalence(wr,wk)
  equivalence(br,bk)
  equivalence(war,wak)

  equivalence(rhsr,rhs)

  equivalence(psir,psik)
  equivalence(qr,qk)
  equivalence(nqr,nqk)



  double precision, dimension(n1d,n2d) :: array2dr
  double complex,   dimension(iktx,ikty) :: array2di

  double precision, dimension(n3)   :: fr_even,fk_even
  double precision, dimension(n3-1) :: fr_odd ,fk_odd

  equivalence(fr_even,fk_even)
  equivalence(fr_odd ,fk_odd )
  equivalence(array2dr,array2di)

  !For implicit dissipation
  double precision :: diss             ! nu_H * kH**(2*ilap) delt

  !Rotational part of u for slice...                                                                                                                                                                                                         
  double complex, dimension(iktx,ikty,n3h1) :: u_rot
  double precision, dimension(n1d,n2d,n3h1) :: u_rotr

  equivalence(u_rotr,u_rot)




  !********************** Initializing... *******************************!


  iter=0

  call initialize_mpi
  call init_files
  call initialize_fftw(array2dr,array2di,fr_even,fk_even,fr_odd,fk_odd)
  call init_arrays
  call init_base_state
  if(mype==0)  call validate_run

  if(init_vertical_structure==generic) then 
     call init_psi_generic(uk,vk,wk,bk,psik,psir)
     call init_q(qk,psik)
  elseif(init_vertical_structure==smith_bernard) then
     call init_psi_generic(uk,vk,wk,bk,psik,psir)   !Init a fully general psi
     call init_q_sb(rhs,psik)       !Compute the jump-dominated PV as in Smith & Bernard (2013)
     call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
     call psi_solver(psik,qt)
     call compute_velo(uk,vk,wk,bk,psik)

    !Set q from normalized rhs                                                                                                                                                                                                      
    do izh0=1,n3h0
       izh1=izh0+1
       do iky=1,ikty
          do ikx=1,iktx
             if (L(ikx,iky).eq.1) then
                 qk(ikx,iky,izh1) =  rhs(ikx,iky,izh0)
              else
                 qk(ikx,iky,izh1) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo
  end if

 if(norm_trop==1) call normalize_trop(uk,vk,wk,bk,psik,qk,wak)

 call generate_halo(uk,vk,wk,bk)
 call generate_halo_q(qk) 

 qok=qk


 !Initial diagnostics!
 !*******************!

 !Compute war/wak if desired                                                                                                                                                     
 if(out_omega==1)  then
    call omega_eqn_rhs(rhs,rhsr,psik)
    call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
    call omega_equation(wak,qt)
 end if

 if(out_etot ==1) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)


 do id_field=1,nfields                                            
    if(out_slice ==1)  call slices(uk,vk,wk,bk,wak,u_rot,ur,vr,wr,br,war,u_rotr,id_field)
 end do
 
! if(out_slab == 1 ) then
!    if(mype==slab_mype) call print_slab(uk,vk)
!    if(mype==slab_mype) call slab_klist
! end if

 if(out_slab == 1) call set_klist
 if(out_slab == 1 .and. mype==slab_mype) call print_time_series(uk,vk)


 if(out_eta == 1 ) call tropopause(uk,vk,wk,bk,ur,vr,wr,br)



 !************************************************************************!
 !*** 1st time timestep using the projection method with Forward Euler ***!
 !************************************************************************!
 
 time=delt
 if(itermax>0) then
! if(mype==0) write(*,*) "First time step"

 iter=1
 
 call convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)                                                                          

 if(linear==1) then
    nqk=(0.D0,0.D0)
 end if

 !Compute dissipation 
 call dissipation_q_nv(dqk,qok)
 
 if(inviscid==1) then
    dqk=(0.D0,0.D0)
 end if

 !Compute q^1 with Forward Euler 
 !First compute u* and b^1
 !u* = u^0 + dt Fu^0 + dt Du^0
 !b^1= b^0 + dt Fb^0 + dt Db^0
 
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss = nuh*delt*(1.*kh2)**(1.*ilap)              !This does not work !!!!! diss = nuh*(kh2**ilap)*delt 
              if (L(ikx,iky).eq.1) then
                 qk(ikx,iky,izh1) = ( qok(ikx,iky,izh1) - delt*nqk(ikx,iky,izh0)  + delt*dqk(ikx,iky,izh0) )*exp(-diss)
                 rhs(ikx,iky,izh0)= qk(ikx,iky,izh1)
              else
                 qk(ikx,iky,izh1) = (0.D0,0.D0)
                 rhs(ikx,iky,izh0)= (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo


 !Generate halo for q
 call generate_halo_q(qk)

 !Transpose rhs -> ft                                                                                                                                                                                                          
 call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)

 !Solve the pressure equation laplacian(phi)=f                                                                                                                                                                                      
 call psi_solver(psik,qt)

 !Compute the corresponding u,v,w and t (u and v to be used in convol)                                                                                                                                    
 call compute_velo(uk,vk,wk,bk,psik)
 call generate_halo(uk,vk,wk,bk)


! if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. mype==slab_mype) call print_slab(uk,vk)                                                                                                                                               
 if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. mype==slab_mype) call print_time_series(uk,vk)




end if



 !********************************************************************************!
 !*** Subsequent timesteps using the projection method + leapfrog timestepping ***!
 !********************************************************************************!


!  if(mype==0) write(*,*) "Subsequent time steps"
  do iter=2,itermax

!     if(mype==0)  cputime=etime(tarray1)
     
     time=iter*delt

     call convol_q(nqk,nqr,uk,vk,qk,ur,vr,qr)                                                                          

     if(linear==1) then
        nqk=(0.D0,0.D0)
     end if

     !Compute dissipation from qok
     call dissipation_q_nv(dqk,qok)

     if(inviscid==1) then
        dqk=(0.D0,0.D0)
     end if


     !Compute q^n+1
     
     do izh0=1,n3h0
        izh1=izh0+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss = nuh*delt*(1.*kh2)**(1.*ilap)              !This does not work !!!!! diss = nuh*(kh2**ilap)*delt                                                                                                          
              if (L(ikx,iky).eq.1) then
                 qtempk(ikx,iky,izh1) =  qok(ikx,iky,izh1)*exp(-2*diss) - 2*delt*nqk(ikx,iky,izh0)*exp(-diss)  + 2*delt*dqk(ikx,iky,izh0)*exp(-2*diss)
                 rhs(ikx,iky,izh0) = qtempk(ikx,iky,izh1)
              else
                 qtempk(ikx,iky,izh1) = (0.D0,0.D0)
                 rhs(ikx,iky,izh0) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo


 !Filter shitz out.
 do izh0=1,n3h0
    izh1=izh0+1
    do iky=1,ikty
       do ikx=1,iktx
          if (L(ikx,iky).eq.1) then
             qok(ikx,iky,izh1) =  qk(ikx,iky,izh1) + gamma * ( qok(ikx,iky,izh1) - 2 * qk(ikx,iky,izh1) + qtempk(ikx,iky,izh1) )
          else
             qok(ikx,iky,izh1) = (0.D0,0.D0)
           endif
       enddo
    enddo
 enddo

!Overwrite the new field uk with u^{n+1} 
 qk = qtempk

 !Generate halo for q
 call generate_halo_q(qk)
 call generate_halo_q(qok)
 

 !Transpose rhs -> ft                                                                                                                                                               
 call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)

 !Solve the pressure equation laplacian(phi)=f                                                                                                                                                                                      
 call psi_solver(psik,qt)

 !Compute the corresponding u,v,w and t 
 call compute_velo(uk,vk,wk,bk,psik)
 call generate_halo(uk,vk,wk,bk) 


 !*** Diagnostics ***!
 !-------------------!

 !Compute w if desired
 if(out_omega==1 .and. (mod(iter,freq_omega) ==0))  then
    call omega_eqn_rhs(rhs,rhsr,psik)
    call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
    call omega_equation(wak,qt)
    call generate_halo_q(wak)
 end if
 


if(out_etot ==1 .and. mod(iter,freq_etot )==0) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)

 do id_field=1,nfields
    if(out_slice ==1 .and. mod(iter,freq_slice)==0 .and. count_slice(id_field)<max_slices)  call slices(uk,vk,wk,bk,wak,u_rot,ur,vr,wr,br,war,u_rotr,id_field)
 end do

 if(out_eta == 1 .and. mod(iter,freq_eta)==0 ) call tropopause(uk,vk,wk,bk,ur,vr,wr,br)


! if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. mype==slab_mype) call print_slab(uk,vk)
 if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. mype==slab_mype) call print_time_series(uk,vk)


 if(time>maxtime) EXIT
end do !End loop

 if(mype==0)  write(*,*) n1,ave_cpu/(1.*itermax-1.)

!************ Terminating processes **********************!                                                                                                                         

  call kill_fftw                                                                                                                                              
  call kill_mpi                                                                                                                                  
 
END PROGRAM main
