class nec_structure_currents
{
public:

	/*! Returns the frequency in Herz. */
	nec_float get_frequency();
	
	
	/*! Returns the flag which controls the printing of the currents. */
	int get_iptflg();
	
	
	/*! Returns the flag which controls the printing of charge densities. */
	int get_iptflq();
	
	
	/*! Returns the number the wire segments in the geometry. */
	int get_n();
	
	
	/*! Returns the number of patches in the geometry. */
	int get_m();
	
	
	/*! Returns the array of segment numbers for the printing of currrents. */
	vector<int> get_current_segment_number();
	
	
	/*! Returns the array of segment tag numbers for the printing of currents, provided the standard output format has been requested. */
	vector<int> get_current_segment_tag();
	
	
	/*! Returns the array of x-coordinate of segment centers in meters for the printing of currents, provided the standard output format has been requested. */
	vector<nec_float> get_current_segment_center_x();
	
	
	/*! Returns the array of y-coordinate of segment centers in meters for the printing of currents, provided the standard output format has been requested. */
	vector<nec_float> get_current_segment_center_y();
	
	
	/*! Returns the array of z-coordinate of segment centers in meters for the printing of currents, provided the standard output format has been requested. */
	vector<nec_float> get_current_segment_center_z();
	
	
	/*! Returns the array of segment lengths in meters for the printing of currents, provided the standard output format has been requested. */
	vector<nec_float> get_current_segment_length();
	
	
	/*! Returns the array of theta angles in degrees for the printing of currents, provided the format designed for a receiving pattern has been requested. */
	vector<nec_float> get_current_theta();
	
	
	/*! Returns the array of phi angles in degrees for the printing of currents, provided the format designed for a receiving pattern has been requested. */	
	vector<nec_float> get_current_phi();
	
	
	/*! Returns the array of complex currents in Ampere. */
	vector<nec_complex> get_current();
	
	
	/*! Returns the array of segment numbers for the printing of charge densities. */
	vector<int> get_q_density_segment_number();
	
	
	/*! Returns the array of segment tag numbers for the printing of charge densities. */
	vector<int> get_q_density_segment_tag();
	
	
	/*! Returns the array of x-coordinate of segment centers in meters for the printing of charge densities. */
	vector<nec_float> get_q_density_segment_center_x();
	
	
	/*! Returns the array of y-coordinate of segment centers in meters for the printing of charge densities. */
	vector<nec_float> get_q_density_segment_center_y();
	
	
	/*! Returns the array of z-coordinate of segment centers in meters for the printing of charge densities. */
	vector<nec_float> get_q_density_segment_center_z();
	
	
	/*! Returns the array of segment lengths in meters for the printing of charge densities. */
	vector<nec_float> get_q_density_segment_length();
	
	
	/*! Returns the array of complex charge densities in Coulomb/meter. */
	vector<nec_complex> get_q_density();
	
	
	/*! Returns the array of patch numbers. */
	vector<int> get_patch_number();
	
	
	/*! Returns the array of x-coordinate of patch centers. */
	vector<nec_float> get_patch_center_x();
	
	
	/*! Returns the array of y-coordinate of patch centers. */
	vector<nec_float> get_patch_center_y();
	
	
	/*! Returns the array of z-coordinate of patch centers. */
	vector<nec_float> get_patch_center_z();
	
	
	/*! Returns the array of complex tangent vector 1 of the patches. */
	vector<nec_complex> get_patch_tangent_vector1();
	
	
	/*! Returns the array of complex tangent vector 2 of the patches. */
	vector<nec_complex> get_patch_tangent_vector2();
	
	
	/*! Returns the complex x-component of the electric field E. */
	vector<nec_complex> get_patch_e_x();
	
	
	/*! Returns the complex y-component of the electric field E. */
	vector<nec_complex> get_patch_e_y();
	
	
	/*! Returns the complex z-component of the electric field E. */
	vector<nec_complex> get_patch_e_z();
	
	
	
/*this private method won't be wrapped, but allow an error in the compilation
process to be avoided.*/ 	
private:
	
	nec_structure_currents(nec_context * in_context, char * in_pattype,
	int in_nload,
	nec_float in_xpr3, nec_float in_xpr6);
		
};
