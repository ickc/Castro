/* problem-specific stuff goes here */

extern "C"
{

  void problem_checkpoint(int * int_dir_name, int * len);

  void problem_restart(int * int_dir_name, int * len);

  void flame_width_temp(BL_FORT_FAB_ARG_3D(temp), 
			const int* lo, const int* hi,
			const Real* dx, const Real* time,
			Real* T_max, Real* T_min, Real* grad_T_max);

  void flame_speed_data(BL_FORT_FAB_ARG_3D(omegadot), 
			const int* lo, const int* hi,
			const Real* dx,
			Real* rho_X_dot);

}
