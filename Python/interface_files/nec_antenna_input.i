class nec_antenna_input
{
		
public:
	
	/*! Returns the frequency in Herz. */
	nec_float get_frequency();
	
	
	/*! Returns the array of segment tag numbers. */
	vector<int> get_tag();
	
	
	/*! Returns the array of segment numbers. */
	vector<int> get_segment();
	
	
	/*! Returns the array of complex currents in Ampere. */	
	vector<nec_complex> get_current();
	
	
	/*! Returns the array of complex voltages in Volt. */
	vector<nec_complex> get_voltage();
	
	
	/*! Returns the array of power in Watt. */
	vector<nec_float> get_power();	
};
