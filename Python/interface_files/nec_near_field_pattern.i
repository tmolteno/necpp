class nec_near_field_pattern
{
public:
	
	/*! Returns the frequency in Herz. */
	nec_float get_frequency();
	
	
	/*! Returns the flag indicating whether the result is a near electric or magnetic field pattern. */
	int get_nfeh();
	
	
	/*! Returns the array of x-coordinate in meters of field points. */  
	vector<nec_float> get_x();
	
	
	/*! Returns the array of y-coordinate in meters of field points. */
	vector<nec_float> get_y();
	
	
	/*! Returns the array of z-coordinate in meters of field points. */
	vector<nec_float> get_z();
	
	
	/*! Returns the array of x_components of the electric or magnetic field. */
	vector<nec_complex> get_field_x();
	
	
	/*! Returns the array of y_components of the electric or magnetic field. */
	vector<nec_complex> get_field_y();
	
	
	/*! Returns the array of z_components of the electric or magnetic field. */
	vector<nec_complex> get_field_z();
	
	

/*this private method won't be wrapped, but allow an error in the compilation
process to be avoided.*/ 	
private:
	
	nec_near_field_pattern(int nfeh);
		
};
