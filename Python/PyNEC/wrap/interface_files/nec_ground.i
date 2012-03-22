%nodefault;
class nec_ground
{
public:
	%extend
	{
		/*! Returns the relative dielectric constant (no units) of the ground medium 1. */
		nec_float get_relative_dielectric_constant()
		{
			return self->epsr;
		}
		
		
		/*! Returns the conductivity in Siemens/meter of the ground medium 1. */
		nec_float get_conductivity()
		{
			return self->sig;
		}
		
		
		/*! Returns the number of radial wires in the ground screen approximation. If it's zero then this approximation has not been used.*/
		int get_radial_wire_count()
		{
			return self->radial_wire_count;
		}
		
		
		/*! Returns the length of radial wires used in the ground screen approximation - provided this approximation has been used. */
		nec_float get_radial_wire_length()
		{
			return self->radial_wire_count;
		}
		
		
		/*! Returns the radius of radial wires in the ground screen approximation - provided this approximation has been used. */
		nec_float get_radial_wire_radius()
		{
			return self->radial_wire_radius;
		}
		
		
		/*! If there's a cliff problem, returns the distance from the origin of the coordinate system to join between medium 1 and 2.
			This distance is either	the radius of the circle where the two media join or the distance from the X axis to where
			the two media join in a line parallel to the Y axis. Specification of the circular or linear option is on the RP card.
		*/		
		nec_float get_cliff_edge_distance()
		{
			return self->cliff_edge_distance;
		}
		
		
		/*! If there's a cliff problem, returns the distance (positive or zero) by which the surface of medium 2 is below medium 1. */
		nec_float get_cliff_height()
		{
			return self->cliff_height;
		}
		
		
		/*! If there's a cliff problem, returns the relative dielectric constant (no units) of the ground medium 2. */
		nec_float get_relative_dielectric_constant2()
		{
			return self->epsr2;
		}
		
		
		/*! If there's a cliff problem, returns the conductivity in Siemens/meter of the ground medium 2. */
		nec_float get_conductivity2()
		{
			return self->sig2;
		}		
	}
};
