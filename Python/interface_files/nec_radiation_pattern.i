%nodefault;
class nec_radiation_pattern
{
public:
	
	/*! Returns the frequency in Herz. */
	nec_float get_frequency();
	
	
	/*! Returns the associated ground object. */
	nec_ground get_ground();
	
	
	/*! Returns the radial distance in meters ( or the rho cylindrical coordinate in meters if the calculation mode chosen is mode 1 ). */
	nec_float get_range();
		
	
	/*! Returns the array of gains in dB used in the normalization process. */
	real_array get_gain();
	
	
	/*! Returns the array of vertical (or major axis, depending on the output format chosen) gains in dB. */
	real_array get_gain_vert();
	
	
	/*! Returns the array of horizontal (or minor axis, depending on the output format chosen) gains in dB*/
	real_array get_gain_horiz();
	
	
	/*! Returns the array of total gains in dB*/
	real_array get_gain_tot();
	
	
	/*! Returns the array of polarization axial ratios (no units). */
	real_array get_pol_axial_ratio();
	
	
	/*! Returns the array of polarization tilts in degrees. */
	real_array get_pol_tilt();
	
	
	/*! Returns the array of polarization sense indexes (no units). The relationship between the index and the actual sense is the following :
		0 : linear
		1 : right 
		2 : left
	*/
	int_array get_pol_sense_index();
	
	
	/*! Returns the array of complex theta-components of electric field E in Volt/meter. */
	complex_array get_e_theta();
	
	
	/*! Returns the array of complex phi-components of electric field E in Volt/meter. */
	complex_array get_e_phi();
	
	
	/*! Returns the array of complex radial-components of electric field E in Volt/meter - only available for the calculation mode 1. */
	complex_array get_e_r();
	
	
	/*! Returns the normalization factors in dB provided a normalization has been requested. */
	nec_float get_normalization_factor();
	
	
	/*! Returns the increment for theta in degrees (or for z in meters if the calculation mode chosen is mode 1 ). */
	nec_float get_delta_theta();
	
	
	/*! Returns the first value of theta in degrees (or of z in meters if the calculation mode chosen is mode 1 ). */
	nec_float get_theta_start();
	
	
	/*! Returns the increment for phi in degrees. */
	nec_float get_delta_phi();
	
	
	/*! Returns the first value of phi in degrees. */
	nec_float get_phi_start();
	
	
	/*! Returns the number of theta angles. */
	int get_ntheta() const;
	
	
	/*! Returns the number of phi angles. */
	int get_nphi() const;
	
	
	/*! Returns the array of average power gains in dB, provided its computation has been requested. */
	nec_float get_average_power_gain();
	
	
	/*! Returns the solid angle in steradians used in the averaging process, provided the computation of an average gain has been requested. */ 
	nec_float get_average_power_solid_angle();
	
	
	/*! Returns the flag (no units) which indicates the calculation mode chosen. */
	int get_ifar();
	
	
	/*! Returns the flag (no units) which indicates the target of the normalization process. */ 
	int get_rp_normalization();
	
	
	/*! Returns the flag (no units) which indicates the output format chosen. */
	int get_rp_output_format();
	
	
	/*! Returns the flag (no units) which indicates whether the average gain will be computed or not. */
	int get_rp_power_average();
	
	
	/*! Returns the flag (no units) which indicates the type of gain computed : power or directive gain. */
	int get_rp_ipd();
};
