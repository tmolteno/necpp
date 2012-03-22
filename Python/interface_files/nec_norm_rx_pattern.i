%nodefault;
class nec_norm_rx_pattern
{
public:
	
	/*! Returns the frequency in Herz. */
	nec_float get_frequency();
	
	
	/*! Returns the number of theta angles. */
	int get_n_theta();
	
	
	/*! Returns the number of phi angles. */
	int get_n_phi();
	
	
	/*! Returns the first value of theta in degrees. */
	nec_float get_theta_start();
	
	
	/*! Returns the first value of phi angles in degrees. */
	nec_float get_phi_start();
	
	
	/*! Returns the increment for theta in degrees. */
	nec_float get_delta_theta();
	
	
	/*! Returns the increment for phi in degrees. */
	nec_float get_delta_phi();
	
	
	/*! Returns the value of eta in degrees. */
	nec_float get_eta();
	
	
	/*! Returns the axial ratio (no units). */
	nec_float get_axial_ratio();
	
	
	/*! Returns the segment number. */
	int get_segment_number();
	
	
	/*! Return the polarization type. */
	string get_type();
	
	
	/*! Returns the normalization factor in dB. */
	nec_float get_norm_factor();
	
	
	/*! Returns the array of receiving gains not yet normalized. */
	real_array get_mag();
		
};
