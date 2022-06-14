subroutine ll_init_model(N_params, params)
use example_heisenberg_params_mod
implicit none
   integer :: N_params
   double precision :: params(N_params)

   ! PARAMS: Jmag (exchange parameter), cutoff
   if (N_params /= 2) then
      print *, "example_heisenberg_model.F90 got 2 /= N_params ", N_params
      print *, "Need: Jmag, cutoff"
      call exit(1)
   endif 

   jmag = params(1)
   cutoff = params(2)

   cutoff_sq = cutoff*cutoff

end subroutine ll_init_model

subroutine ll_init_config(N, pos, cell, Emax)
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   double precision :: Emax
   return
end subroutine ll_init_config

subroutine ll_init_config_mag(N, pos, mom, cell, Emax)
implicit none
   integer :: N
   double precision :: pos(3,N), mom(3,N), cell(3,3)
   double precision :: Emax
   return
end subroutine ll_init_config_mag

double precision function ll_eval_energy_mag(N, pos, mom, n_extra_data, extra_data, cell)
use example_mat_mod
use example_heisenberg_params_mod
implicit none
   integer :: N
   double precision :: pos(3,N), pos_2D(2,N), mom(3,N), cell(3,3), cell_2D(2,2)

   integer :: i, j
   double precision :: dr(2), dr_mag, dr_mag_sq, dr_l(2), dr_l0(2), pos_l(2,N), drr(2), dot_moms

   double precision :: cell_inv(2,2), E_term
   integer :: dj1, dj2

   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)

   integer :: n_images
   double precision cell_height(2)

   do i=1, N
      pos_2D(1,i) = pos(1,i)
      pos_2D(2,i) = pos(2,i)
   end do
   do i=1, 2
      cell_2D(1,i) = cell(1,i)
      cell_2D(2,i) = cell(2,i)
   end do

   call matrix2x2_inverse(cell_2D, cell_inv)
   ! into lattice coodinates 

   pos_l = matmul(cell_inv, pos_2D)

   do i=1, 2
      cell_height(i) = sqrt(sum(cell(:,i)**2))
   end do

   n_images = ceiling(cutoff/minval(cell_height))

   ll_eval_energy_mag = 0.0
   do i=1, N
   do j=i, N
      dr_l0 = pos_l(:,i)-pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-n_images,n_images
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-n_images,n_images
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      if (i == j .and. dj1 == 0 .and. dj2 == 0) cycle !  no self-interaction
      dr(1) = cell(1,1) * dr_l(1) + cell(1,2) * dr_l(2) ! sum(cell(1,:)*dr_l)
      dr(2) = cell(2,1) * dr_l(1) + cell(2,2) * dr_l(2) ! sum(cell(2,:)*dr_l)
      dr_mag_sq = dr(1)*dr(1) + dr(2)*dr(2) ! sum(dr*dr)
      if (dr_mag_sq < cutoff_sq) then
          dr_mag = sqrt(dr_mag_sq)
          dot_moms = mom(1,i) * mom(1,j) + mom(2,i) * mom(2,j) + mom(3,i) * mom(3,j) 
          E_term = Jmag * (dot_moms)
          if (i == j) E_term = E_term * 0.5
          ll_eval_energy_mag = ll_eval_energy_mag + E_term
      endif

      end do
      end do
   end do
   end do
end function ll_eval_energy_mag

double precision function ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)
use example_mat_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)

   ll_eval_energy = 0.0

end function ll_eval_energy

integer function ll_move_atom_1(N, Z, pos, n_extra_data, extra_data, cell, d_i, d_pos, dEmax, dE)
use example_mat_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   integer :: d_i
   double precision :: d_pos(3)
   double precision :: dEmax, dE

   ll_move_atom_1 = 1
   dE = 0.0

end function ll_move_atom_1

function ll_eval_forces(N, Z, pos, n_extra_data, extra_data, cell, forces) result(energy)
use example_mat_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   double precision :: energy ! result

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_l(3), dr_l0(3), pos_l(3,N)
   double precision :: cell_inv(3,3)
   integer :: dj1, dj2, dj3

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6, E_term

   energy = 0.0
   forces = 0.0

end function ll_eval_forces
