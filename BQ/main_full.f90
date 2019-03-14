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
  double precision, dimension(n1d,n2d,n3h1)   :: nur,nvr,nwr,nbr     !Vorticity and nonlinear fields (r-space)
  double complex,   dimension(iktx,ikty,n3h2) :: uk,vk,wk,bk                     !Velocity and potential temperature fields (k-space)
  double complex,   dimension(iktx,ikty,n3h1) :: nuk,nvk,nwk,nbk     !Vorticity and nonlinear fields (k-space)
  double complex,   dimension(iktx,ikty,n3h1) :: psik,qk                !streamfunction
  double precision, dimension(n1d,n2d,n3h1)   :: psir
  double precision, dimension(n1d,n2d,n3h1)   :: qr
  double complex, dimension(iktx,n3, iktyp) :: qt             !Transposed (ky-parallelization) right-hand side         (For TS init only)

  double complex,   dimension(iktx,ikty,n3h1) :: wak
  double precision, dimension(n1d,n2d,n3h1)   :: war

  equivalence(ur,uk)
  equivalence(vr,vk)
  equivalence(wr,wk)
  equivalence(br,bk)

  equivalence(qr,qk)
  equivalence(psir,psik)
  equivalence(war,wak)

  equivalence(nur,nuk)
  equivalence(nvr,nvk)
  equivalence(nwr,nwk)
  equivalence(nbr,nbk)

  double complex,   dimension(iktx,ikty,n3h2) :: uok,vok,wok,bok        !Old fields 

  double complex,   dimension(iktx,ikty,n3h2) :: utempk,vtempk,wtempk,btempk 
  double precision, dimension(n1d,n2d,n3h2)   :: utempr,vtempr,wtempr,btempr



  equivalence(utempr,utempk)
  equivalence(vtempr,vtempk)
  equivalence(wtempr,wtempk)
  equivalence(btempr,btempk)


  double complex,   dimension(iktx,ikty,n3h1) :: duk,dvk,dwk,dbk         !dissipation
  double complex,   dimension(iktx,ikty,n3h1) :: usk,vsk,wsk             !u*!
  double complex,   dimension(iktx,ikty,n3h1) :: phi                     !pressure, and rhs of pressure equation!
  double complex,   dimension(iktx,ikty,n3h0) :: phi_x,phi_y,phi_z,rhs   !grad(phi)
  double precision, dimension(n1d,n2d,n3h0)   :: rhsr

  equivalence(rhsr,rhs)

  double complex, dimension(iktx,n3, iktyp) :: ft             !Transposed (ky-parallelization) right-hand side (HMMM?)  

  double precision :: diss_u,diss_t             ! nu_H * kH**(2*ilap) delt  (one for nuh, the other for nuth)                                                                                                                                                               
  double complex,   dimension(iktx,ikty) :: ws_bot,ws_top      !Top and bottom w* are the Neumann boundary conditions for pressure


  !********************** Test-only declarations **************************!

  
  double precision, dimension(n1d,n2d,n3h2)   :: ftest1 ,ftest2 ,ftest3
  double complex,   dimension(iktx,ikty,n3h2) :: ftest1k,ftest2k,ftest3k

  double precision, dimension(n1d,n2d,n3h2)   :: ftest1s ,ftest2s ,ftest3s
  double complex,   dimension(iktx,ikty,n3h2) :: ftest1sk,ftest2sk,ftest3sk

  double precision, dimension(n1d,n2d) :: array2dr 
  double complex,   dimension(iktx,ikty) :: array2di

  double precision, dimension(n3)   :: fr_even,fk_even
  double precision, dimension(n3-1) :: fr_odd ,fk_odd

  equivalence(fr_even,fk_even)
  equivalence(fr_odd ,fk_odd )

  equivalence(ftest1,ftest1k)
  equivalence(ftest2,ftest2k)
  equivalence(ftest3,ftest3k)

  equivalence(ftest1s,ftest1sk)
  equivalence(ftest2s,ftest2sk)
  equivalence(ftest3s,ftest3sk)

  equivalence(array2dr,array2di)
  
  !Ageostrophic field variables
!  double complex, dimension(iktx,ikty,n3h1) :: u_a,v_a,w_a,b_a
!  double precision, dimension(n1d,n2d,n3h1) :: u_ar,v_ar,w_ar,b_ar

!  equivalence(u_ar,u_a)
!  equivalence(v_ar,v_a)
!  equivalence(w_ar,w_a)
!  equivalence(b_ar,b_a)

  !Rotational part of u for slice...
  double complex, dimension(iktx,ikty,n3h1) :: u_rot
  double precision, dimension(n1d,n2d,n3h1) :: u_rotr

  equivalence(u_rotr,u_rot)

  integer :: unit_u = 123451234
  integer :: unit_v = 123451235
  integer :: unit_wke = 123451236

  double precision, dimension(n1d,n2d,n3h2)   :: uwr,vwr     !(approximate wave part of the flow)
  double complex,   dimension(iktx,ikty,n3h2) :: uwk,vwk
  
  equivalence(uwr,uwk)
  equivalence(vwr,vwk)

  double precision, dimension(n1d,n2d,1)   :: ubr,vbr     !barotropic part of the flow
  double complex,   dimension(iktx,ikty,1) :: ubk,vbk

  equivalence(ubr,ubk)
  equivalence(vbr,vbk)

  !To calculate the error between the steady solution and the barotropic component
  real :: error_drift_u,error_drift_v,error_drift_E
  real :: norm_u,norm_v,norm_E 

  !********************** Initializing... *******************************!


  iter=0


  call initialize_mpi
  call init_files
  call initialize_fftw(array2dr,array2di,fr_even,fk_even,fr_odd,fk_odd)
  call init_arrays
  call init_base_state
  if(mype==0)  call validate_run

  call generate_fields_stag(ur,n3h2,vr,n3h2,wr,n3h2)
  call fft_r2c(ur,uk,n3h2)
  call fft_r2c(vr,vk,n3h2)
  wk=(0.D0,0.D0)
  bk=(0.D0,0.D0)

  uok=uk
  vok=vk
  wok=wk
  bok=bk

  open (unit=unit_u,file='u.dat',action="write",status="replace")
  open (unit=unit_v,file='v.dat',action="write",status="replace")
  open (unit=unit_wke,file='wke.dat',action="write",status="replace")

  !Initial diagnostics!
  !-------------------!

     !**********************************************************************************************************!
     !*** Get the approximate wave part of the flow by subtracting the otherwise steady part of the solution ***!
     !**********************************************************************************************************!


     call fft_c2r(uk,ur,n3h2)
     call fft_c2r(vk,vr,n3h2)
     
     do ix=1,n1
        x = xa(ix)
        do iy=1,n2
           y = ya(iy)
           do izh2=1,n3h2
              
              uwr(ix,iy,izh2) = ur(ix,iy,izh2) + sin(x)*cos(y)
              vwr(ix,iy,izh2) = vr(ix,iy,izh2) - cos(x)*sin(y)
              
           end do
        end do
     end do
     
!     if(mype==(npe-1)) then
!        write(unit_u,fmt=*)   time,time/(twopi*Ro),uwr(n1/4,n2/4,iztop2)*U_scale
!        write(unit_v,fmt=*)   time,time/(twopi*Ro),vwr(n1/4,n2/4,iztop2)*U_scale
!        write(unit_wke,fmt=*) time,time/(twopi*Ro),0.5*(uwr(n1/4,n2/4,iztop2)**2+vwr(n1/4,n2/4,iztop2)**2)*U_scale*U_scale
!     end if
     
     
     call fft_r2c(ur,uk,n3h2)
     call fft_r2c(vr,vk,n3h2)
     

     !**********************************************************************************************************!



  !Compute war/wak if desired                                                                                                                                                                                                  
  if(out_omega==1 .and. init_wageo/=1)  then
     call omega_eqn_rhs(rhs,rhsr,psik)
     call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
     call omega_equation(wak,qt)
  end if

  if(out_etot ==1) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)
  if(out_cont ==1) call continuity_anelastic(uk,vk,wk)

 ! Compute q (linear PV)                                                                                                                                                                                                                    
 call init_q(qk,psik)
 do id_field=1,nfields
    if(out_slice ==1)  call slices(uk,vk,wk,bk,wak,u_rot,ur,vr,wr,br,war,u_rotr,uwr,vwr,id_field)
 end do

 if(out_slab == 1) call set_klist
 if(out_slab == 1 .and. mype==slab_mype .and. many_slab == 0) call print_time_series(uk,vk)
 if(out_slab == 1 .and. many_slab == 1) then
    do iz=1,num_spec
       call print_many_time_series(uk,vk,iz)   
    end do
 end if

 if(out_3d ==1) call output3d(uk,vk,wk,bk,ur,vr,wr,br)


 if(out_eta == 1) call tropopause(uk,vk,wk,bk,ur,vr,wr,br)

 !************************************************************************!                                                                                                                                            
 !*** 1st time timestep using the projection method with Forward Euler ***!                                                                                                                                         
 !************************************************************************!                                                                                                                                               

 time=delt
 if(itermax>0) then

 iter=1

 call convol(nuk,nvk,nwk,nbk,nur,nvr,nwr,nbr,uk,vk,wk,bk,ur,vr,wr,br)

 if(linear==1) then
    nuk=(0.D0,0.D0)
    nvk=(0.D0,0.D0)
    nwk=(0.D0,0.D0)
    nbk=(0.D0,0.D0)
 end if

 !Compute dissipation                                                                                                                                                                                                                
 call dissipation_nv(duk,dvk,dwk,dbk,uk,vk,wk,bk)

 if(inviscid==1) then
    duk=(0.D0,0.D0)
    dvk=(0.D0,0.D0)
    dwk=(0.D0,0.D0)
    dbk=(0.D0,0.D0)
 end if


 !First compute u* and b^1                                                                                                                                                                                                    
 !u* = u^0 + dt Fu^0 + dt Du^0                                                                                                                                                                                                          
 !b^1= b^0 + dt Fb^0 + dt Db^0                                                                                                                                                                                                           

     do izh1=1,n3h1
        izh2=izh1+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss_u = nuh *delt*(1.*kh2)**(1.*ilap)
              diss_t = nuth*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
                 usk(ikx,iky,izh1) =  (uok(ikx,iky,izh2) - delt*nuk(ikx,iky,izh1) + delt*(1./Ro)                                              *vok(ikx,iky,izh2) + delt*duk(ikx,iky,izh1))*exp(-diss_u)
                 vsk(ikx,iky,izh1) =  (vok(ikx,iky,izh2) - delt*nvk(ikx,iky,izh1) - delt*(1./Ro)                                              *uok(ikx,iky,izh2) + delt*dvk(ikx,iky,izh1))*exp(-diss_u)
                 wsk(ikx,iky,izh1) =  (wok(ikx,iky,izh2) - delt*nwk(ikx,iky,izh1) + delt*(1./(Ar2*Ro))*r_1(izh2)*0.5*( bok(ikx,iky,izh2+1) + bok(ikx,iky,izh2) ) + delt*dwk(ikx,iky,izh1))*exp(-diss_u)
                  bk(ikx,iky,izh2) =  (bok(ikx,iky,izh2) - delt*nbk(ikx,iky,izh1) - delt*(Ro/(Fr*Fr)) *r_2(izh2)*0.5*( wok(ikx,iky,izh2) + wok(ikx,iky,izh2-1) ) + delt*dbk(ikx,iky,izh1))*exp(-diss_t)
              else
                 usk(ikx,iky,izh1) = (0.D0,0.D0)
                 vsk(ikx,iky,izh1) = (0.D0,0.D0)
                 wsk(ikx,iky,izh1) = (0.D0,0.D0)
                  bk(ikx,iky,izh2) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo


     !Boundaries!
     !----------!

     if(mype==0) then   !Watch out: wsk is defined at z=0 (izbot1-1) too
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss_u = nuh *delt*(1.*kh2)**(1.*ilap)
              diss_t = nuth*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
                  bk(ikx,iky,izbot2)   =  (bok(ikx,iky,izbot2) - delt*nbk(ikx,iky,izbot1) - delt*(Ro/(Fr*Fr)) *r_2(izbot2)*0.5*( wok(ikx,iky,izbot2) ) + delt*dbk(ikx,iky,izbot1))*exp(-diss_t)
                 wsk(ikx,iky,izbot1-1) =  ( delt*(1./(Ar2*Ro))*r_1(izbot2-1)*bok(ikx,iky,izbot2) )*exp(-diss_u)    !This requires that we define r_1 at z=0 (generally in the X-marked regions.
              else
                  bk(ikx,iky,izbot2) = (0.D0,0.D0)
                 wsk(ikx,iky,izbot1-1) = (0.D0,0.D0)
              endif

              ws_bot(ikx,iky) = Ar2*wsk(ikx,iky,izbot1-1)
              
           enddo
        enddo
     elseif(mype==(npe-1)) then   !Modify top w* (WE MUST KEEP HORIZONTAL DIFFUSION)
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss_u = nuh *delt*(1.*kh2)**(1.*ilap)
              diss_t = nuth*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
                 wsk(ikx,iky,iztop1) =  ( delt*(1./(Ar2*Ro))*r_1(iztop2)*bok(ikx,iky,iztop2) )*exp(-diss_u)
              else
                 wsk(ikx,iky,iztop1) =  (0.D0,0.D0)
              endif

              ws_top(ikx,iky) = Ar2*wsk(ikx,iky,iztop1)              

           enddo
        enddo
     end if

     !Broadcast top/bottom w* (they are the boundary conditions on dphi/dz) 
     call mpi_bcast(ws_bot,iktx*ikty,MPI_DOUBLE_COMPLEX,0    ,MPI_COMM_WORLD,ierror)
     call mpi_bcast(ws_top,iktx*ikty,MPI_DOUBLE_COMPLEX,npe-1,MPI_COMM_WORLD,ierror)


     !Compute the gradient of u*, store in rhs                                                                                                                                                                  
     call divergence(usk,vsk,wsk,rhs)

     !Transpose rhs -> ft                                                                                                                                                                                                  
     call mpitranspose(rhs,iktx,ikty,n3h0,ft,n3,iktyp)

     !Solve the pressure equation laplacian(phi)=f                                                                                                                                                                  
     call helmholtzdouble(phi,ft,ws_bot,ws_top)

     !Compute the gradient of phi                                                                                                                                                                                     
     call gradient(phi,phi_x,phi_y,phi_z)


 !Get u^n+1 = u* - grad phi                                                                                                                                        
                                                                                                                                                                                                                   
 do izh0=1,n3h0
    izh1=izh0+1
    izh2=izh0+2
    do iky=1,ikty
       ky = kya(iky)
       do ikx=1,iktx
          kx = kxa(ikx)
          if (L(ikx,iky).eq.1) then
             uk(ikx,iky,izh2) =  usk(ikx,iky,izh1) - phi_x(ikx,iky,izh0)
             vk(ikx,iky,izh2) =  vsk(ikx,iky,izh1) - phi_y(ikx,iky,izh0)
             wk(ikx,iky,izh2) =  wsk(ikx,iky,izh1) - phi_z(ikx,iky,izh0)/Ar2
          else
             uk(ikx,iky,izh2) = (0.D0,0.D0)
             vk(ikx,iky,izh2) = (0.D0,0.D0)
             wk(ikx,iky,izh2) = (0.D0,0.D0)
          endif
       enddo
    enddo
 enddo

 !Set w to 0 at the bottom                                                                                                                                                                                                          
 if(mype==(npe-1)) then
    do iky=1,ikty
       do ikx=1,iktx

          wk(ikx,iky,iztop2) = (0.D0,0.D0)

       enddo
    enddo
 end if

 !Set mode kh=0 to 0 at all levels                                                                                                                                                                          
 do izh2=1,n3h2
    wk(1,1,izh2) = (0.D0,0.D0)
 end do

 !Generate halo for the new variables                                                                                                                                                                                       
 call generate_halo(uk,vk,wk,bk)


 if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. mype==slab_mype .and. many_slab == 0) call print_time_series(uk,vk)
 if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. many_slab == 1) then
    do iz=1,num_spec
       call print_many_time_series(uk,vk,iz)
    end do
 end if

end if




 !********************************************************************************!                                                                                                                                            
 !*** Subsequent timesteps using the projection method + leapfrog timestepping ***!                                                                                                                                       
 !********************************************************************************!                                                                                                                                        


!  if(mype==0) write(*,*) "Subsequent time steps"                                                                                                                                                                                
  do iter=2,itermax

!     if(mype==0)  cputime=etime(tarray1)                                                                                                                                                                                       
     time=iter*delt

     !Compute the nlt^n                                                                                                                                                                                                        
     call convol(nuk,nvk,nwk,nbk,nur,nvr,nwr,nbr,uk,vk,wk,bk,ur,vr,wr,br)

     if(linear==1) then
        nuk=(0.D0,0.D0)
        nvk=(0.D0,0.D0)
        nwk=(0.D0,0.D0)
        nbk=(0.D0,0.D0)
     end if

     !Compute dissipation                                                                                                                                                                                                    
     call dissipation_nv(duk,dvk,dwk,dbk,uok,vok,wok,bok)

     if(inviscid==1) then
        duk=(0.D0,0.D0)
        dvk=(0.D0,0.D0)
        dwk=(0.D0,0.D0)
        dbk=(0.D0,0.D0)
     end if


     !Get u*    = uold + 2*dt*F^n + 2*dt*D^n-1                                                                                                                                                                      
     !Get b_n+1 = bold + 2*dt*F^n + 2*dt*D^n-1                                                                                                                                                                                              

     do izh1=1,n3h1
        izh2=izh1+1
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss_u = nuh *delt*(1.*kh2)**(1.*ilap)
              diss_t = nuth*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
                 usk(ikx,iky,izh1) =  uok(ikx,iky,izh2)*exp(-2*diss_u) + 2*delt*(- nuk(ikx,iky,izh1) + (1./Ro)                                             *vk(ikx,iky,izh2))*exp(-diss_u) + 2*delt*duk(ikx,iky,izh1)*exp(-2*diss_u)
                 vsk(ikx,iky,izh1) =  vok(ikx,iky,izh2)*exp(-2*diss_u) + 2*delt*(- nvk(ikx,iky,izh1) - (1./Ro)                                             *uk(ikx,iky,izh2))*exp(-diss_u) + 2*delt*dvk(ikx,iky,izh1)*exp(-2*diss_u)
                 wsk(ikx,iky,izh1) =  wok(ikx,iky,izh2)*exp(-2*diss_u) + 2*delt*(- nwk(ikx,iky,izh1) + (1./(Ar2*Ro))*r_1(izh2)*0.5*( bk(ikx,iky,izh2+1) + bk(ikx,iky,izh2) ))*exp(-diss_u) + 2*delt*dwk(ikx,iky,izh1)*exp(-2*diss_u)
              btempk(ikx,iky,izh2) =  bok(ikx,iky,izh2)*exp(-2*diss_t) + 2*delt*(- nbk(ikx,iky,izh1) - (Ro/(Fr*Fr)) *r_2(izh2)*0.5*( wk(ikx,iky,izh2) + wk(ikx,iky,izh2-1) ))*exp(-diss_t) + 2*delt*dbk(ikx,iky,izh1)*exp(-2*diss_t)
              else
                 usk(ikx,iky,izh1) = (0.D0,0.D0)
                 vsk(ikx,iky,izh1) = (0.D0,0.D0)
                 wsk(ikx,iky,izh1) = (0.D0,0.D0)
              btempk(ikx,iky,izh2) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo


     !Boundaries!                                                                                                                                                                                                                
     !----------!                                                                                                                                                                                                                          

     if(mype==0) then   !Watch out: wsk is defined at z=0 (izbot1-1) too                                                                                                                                                      
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss_u = nuh *delt*(1.*kh2)**(1.*ilap)
              diss_t = nuth*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
              btempk(ikx,iky,izbot2)   =  bok(ikx,iky,izbot2)*exp(-2*diss_t) + 2*delt*(- nbk(ikx,iky,izbot1) - (Ro/(Fr*Fr)) *r_2(izbot2)*0.5*( wk(ikx,iky,izbot2)  ))*exp(-diss_t) + 2*delt*dbk(ikx,iky,izbot1)*exp(-2*diss_t)
                 wsk(ikx,iky,izbot1-1) =  ( 2*delt*(1./(Ar2*Ro))*r_1(izbot2-1)*bk(ikx,iky,izbot2) )*exp(-diss_u)    !This requires that we define r_1 at z=0 (generally in the X-marked regions.                                           
              else
              btempk(ikx,iky,izbot2) = (0.D0,0.D0)
                 wsk(ikx,iky,izbot1-1) = (0.D0,0.D0)
              endif

              ws_bot(ikx,iky) = Ar2*wsk(ikx,iky,izbot1-1)

           enddo
        enddo
     elseif(mype==(npe-1)) then   !Modify top w* (WE MUST KEEP HORIZONTAL DIFFUSION)                                                                                                                                                  
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              kh2=kx*kx+ky*ky
              diss_u = nuh *delt*(1.*kh2)**(1.*ilap)
              diss_t = nuth*delt*(1.*kh2)**(1.*ilap)
              if (L(ikx,iky).eq.1) then
                 wsk(ikx,iky,iztop1) =  ( 2*delt*(1./(Ar2*Ro))*r_1(iztop2)*bk(ikx,iky,iztop2) )*exp(-diss_u)
              else
                 wsk(ikx,iky,iztop1) =  (0.D0,0.D0)
              endif

              ws_top(ikx,iky) = Ar2*wsk(ikx,iky,iztop1)

           enddo
        enddo
     end if

     !Broadcast top/bottom w* (they are the boundary conditions on dphi/dz)                                                                                                                                                              
     call mpi_bcast(ws_bot,iktx*ikty,MPI_DOUBLE_COMPLEX,0    ,MPI_COMM_WORLD,ierror)
     call mpi_bcast(ws_top,iktx*ikty,MPI_DOUBLE_COMPLEX,npe-1,MPI_COMM_WORLD,ierror)

     !Compute the gradient of u*, store in rhs                                                                                                                                                                            
     call divergence(usk,vsk,wsk,rhs)

     !Transpose rhs -> ft                                                                                                                                                                                                    
     call mpitranspose(rhs,iktx,ikty,n3h0,ft,n3,iktyp)

     !Solve the pressure equation laplacian(phi)=f                                                                                                                                                                             
     call helmholtzdouble(phi,ft,ws_bot,ws_top)

     !Compute the gradient of phi                                                                                                                                                                                                 
     call gradient(phi,phi_x,phi_y,phi_z)

     !Get u^n+1 = u* - grad phi                                                                                                                                                                                                  
     do izh0=1,n3h0
        izh1=izh0+1
        izh2=izh0+2
        do iky=1,ikty
           ky = kya(iky)
           do ikx=1,iktx
              kx = kxa(ikx)
              if (L(ikx,iky).eq.1) then
                 utempk(ikx,iky,izh2) =  usk(ikx,iky,izh1) - phi_x(ikx,iky,izh0)
                 vtempk(ikx,iky,izh2) =  vsk(ikx,iky,izh1) - phi_y(ikx,iky,izh0)
                 wtempk(ikx,iky,izh2) =  wsk(ikx,iky,izh1) - phi_z(ikx,iky,izh0)/Ar2
              else
                 utempk(ikx,iky,izh2) = (0.D0,0.D0)
                 vtempk(ikx,iky,izh2) = (0.D0,0.D0)
                 wtempk(ikx,iky,izh2) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo


     !Apply Robert-Asselin time filter
     do izh0=1,n3h0
        izh2=izh0+2
        do iky=1,ikty
           do ikx=1,iktx
              if (L(ikx,iky).eq.1) then
                 uok(ikx,iky,izh2) =  uk(ikx,iky,izh2) + gamma * ( uok(ikx,iky,izh2) - 2 * uk(ikx,iky,izh2) + utempk(ikx,iky,izh2) )
                 vok(ikx,iky,izh2) =  vk(ikx,iky,izh2) + gamma * ( vok(ikx,iky,izh2) - 2 * vk(ikx,iky,izh2) + vtempk(ikx,iky,izh2) )
                 wok(ikx,iky,izh2) =  wk(ikx,iky,izh2) + gamma * ( wok(ikx,iky,izh2) - 2 * wk(ikx,iky,izh2) + wtempk(ikx,iky,izh2) )
                 bok(ikx,iky,izh2) =  bk(ikx,iky,izh2) + gamma * ( bok(ikx,iky,izh2) - 2 * bk(ikx,iky,izh2) + btempk(ikx,iky,izh2) )
              else
                 uok(ikx,iky,izh2) = (0.D0,0.D0)
                 vok(ikx,iky,izh2) = (0.D0,0.D0)
                 wok(ikx,iky,izh2) = (0.D0,0.D0)
                 bok(ikx,iky,izh2) = (0.D0,0.D0)
              endif
           enddo
        enddo
     enddo
     
     !Overwrite the new field uk with u^{n+1}                                                                                                                                                                                
     uk = utempk
     vk = vtempk
     wk = wtempk
     bk = btempk

     !Set w to at the top                                                                                                                                                                                                                 
     if(mype==(npe-1)) then
        do iky=1,ikty
           do ikx=1,iktx
              
              wk(ikx,iky,iztop2) = (0.D0,0.D0)
              wok(ikx,iky,iztop2) = (0.D0,0.D0)
              
           enddo
        enddo
     end if
     
     
     !Set mode kh=0 to 0 at all levels                                                                                                                                                                                            
     do izh2=1,n3h2
        wk(1,1,izh2) = (0.D0,0.D0)
        wok(1,1,izh2) = (0.D0,0.D0)
     end do


     !Generate halo for the new variables                                                                                                                                                                                    
     call generate_halo(uk,vk,wk,bk)
     call generate_halo(uok,vok,wok,bok)


     !Diagnostics!                                                                                                                                                                                                               
     !-----------!                                                                                                                                                                                                                       

     !**********************************************************************************************************!
     !*** Get the approximate wave part of the flow by subtracting the otherwise steady part of the solution ***!
     !**********************************************************************************************************!

     call barotropic_flow(uk,vk,ubk,vbk)

     call fft_c2r(uk,ur,n3h2)
     call fft_c2r(vk,vr,n3h2)

     !Compute the wave part of the flow
     if(wave_from_steady==1) then   !Calculate the wave component by subtracting the steady solution

        do ix=1,n1
           x = xa(ix)
           do iy=1,n2
              y = ya(iy)
              do izh2=1,n3h2
                 
                 uwr(ix,iy,izh2) = ur(ix,iy,izh2) + sin(x)*cos(y)
                 vwr(ix,iy,izh2) = vr(ix,iy,izh2) - cos(x)*sin(y)
                 
              end do
           end do
        end do

     else                          !Calculate the wave component by subtracting the barotropic part of the flow

        call fft_c2r(ubk,ubr,1)
        call fft_c2r(vbk,vbr,1)

        do ix=1,n1
           do iy=1,n2
              do izh2=1,n3h2
                 
                 uwr(ix,iy,izh2) = ur(ix,iy,izh2) - ubr(ix,iy,1)
                 vwr(ix,iy,izh2) = vr(ix,iy,izh2) - vbr(ix,iy,1)
                 
              end do
           end do
        end do


        !Compute the drift of the barotropic flow from steady solution (L2-norm error)
        error_drift_u = 0.
        error_drift_v = 0.
        error_drift_E = 0.

        norm_u = 0.
        norm_v = 0.
        norm_E = 0.

        do ix=1,n1
           x = xa(ix)
           do iy=1,n2
              y = ya(iy)

              error_drift_u = error_drift_u + (ubr(ix,iy,1) + sin(x)*cos(y))**2
              error_drift_v = error_drift_v + (vbr(ix,iy,1) - cos(x)*sin(y))**2
              error_drift_E = error_drift_E + (  0.5*(ubr(ix,iy,1)**2 + vbr(ix,iy,1)**2) - 0.5*( (sin(x)*cos(y))**2 +  (cos(x)*sin(y))**2 ) )**2  

              norm_u = norm_u + (sin(x)*cos(y))**2
              norm_v = norm_v + (cos(x)*sin(y))**2
              norm_E = norm_E + 0.5*( (sin(x)*cos(y))**2 +  (cos(x)*sin(y))**2 )

           end do
        end do

        error_drift_u = sqrt(error_drift_u/norm_u)
        error_drift_v = sqrt(error_drift_v/norm_v)
        error_drift_E = sqrt(error_drift_E/norm_E)

        if(mype==0) then
           write(unit_u,fmt=*)   time,time/(twopi*Ro),error_drift_u
           write(unit_v,fmt=*)   time,time/(twopi*Ro),error_drift_v
           write(unit_wke,fmt=*) time,time/(twopi*Ro),error_drift_E
        end if

     end if


           



     
     
     call fft_r2c(ur,uk,n3h2)
     call fft_r2c(vr,vk,n3h2)
     

     !**********************************************************************************************************!


     call compute_streamfunction(uk,vk,psik)

     !Compute wa if desired                                                                                                                                                                                                     
     if(out_omega==1 .and. (mod(iter,freq_omega) ==0))  then
        call omega_eqn_rhs(rhs,rhsr,psik)
        call mpitranspose(rhs,iktx,ikty,n3h0,qt,n3,iktyp)
        call omega_equation(wak,qt)
        !    call generate_halo_q(wak)
     end if

     if(out_etot ==1 .and. mod(iter,freq_etot )==0) call diag_zentrum(uk,vk,wk,bk,wak,psik,u_rot)
     if(out_cont ==1 .and. mod(iter,freq_cont )==0) call continuity_anelastic(uk,vk,wk)
     
     if( out_slice ==1  .and. mod(iter,freq_slice)==0 .and. count_slice(id_field)<max_slices) call init_q(qk,psik)
     do id_field=1,nfields
        if(out_slice ==1 .and. mod(iter,freq_slice)==0 .and. count_slice(id_field)<max_slices)  call slices(uk,vk,wk,bk,wak,u_rot,ur,vr,wr,br,war,u_rotr,uwr,vwr,id_field)
     end do

     !Calculate the wave kinetic energy spectrum
     call fft_r2c(uwr,uwk,n3h2)
     call fft_r2c(vwr,vwk,n3h2)
     do iz=1,num_spec
        if(out_hspecw ==1 .and. mod(iter,freq_hspecw) ==0) call hspec_waves(uwk,vwk,iz)
     end do


     if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. mype==slab_mype .and. many_slab == 0) call print_time_series(uk,vk)
     if(out_slab == 1 .and. mod(iter,freq_slab)==0 .and. many_slab == 1) then
        do iz=1,num_spec
           call print_many_time_series(uk,vk,iz)
        end do
     end if
     

     if(out_eta == 1 .and. mod(iter,freq_eta)==0 ) call tropopause(uk,vk,wk,bk,ur,vr,wr,br)

     if(out_3d ==1 .and. mod(iter,freq_3d)==0 ) call output3d(uk,vk,wk,bk,ur,vr,wr,br)


 if(time>maxtime) EXIT
end do !End loop         

!************ Terminating processes **********************!                                                                                                                         

  call kill_fftw                                                                                                                                              
  call kill_mpi                                                                                                                                  
 
END PROGRAM main
