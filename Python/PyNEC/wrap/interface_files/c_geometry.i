class c_geometry
{
public:
		
	/*! Add a wire to the geometry,

	All coordinates are in meters.

		\param tag_id The tag ID.
		\param segment_count The number of segments.
		
		\param xw1 The x coordinate of the wire starting point.
		\param yw1 The y coordinate of the wire starting point.
		\param zw1 The z coordinate of the wire starting point.
		
		\param xw2 The x coordinate of the wire ending point.
		\param yw2 The y coordinate of the wire ending point.
		\param zw2 The z coordinate of the wire ending point.
		
		\param rad The wire radius (meters)
		\param rdel For tapered wires, the. Otherwise set to 1.0
		\param rrad For tapered wires, the. Otherwise set to 1.0
	*/
	void wire( 	int tag_id, int segment_count,
			nec_float xw1, nec_float yw1, nec_float zw1,
			nec_float xw2, nec_float yw2, nec_float zw2,
			nec_float rad,
			nec_float rdel, nec_float rrad);
	
	
			
	/*! Add an arc to the geometry,

	All coordinates are in meters and angles are in degrees.

		\param tag_id The tag ID.
		\param segment_count The number of segments.
		
		\param rada The radius.
		\param ang1 The angle of the arc starting point.
		\param ang2 The angle of the arc end point.
		\param rad The wire radius.
	*/		
	void arc( int tag_id, int segment_count, nec_float rada,
			nec_float ang1, nec_float ang2, nec_float rad );
	
	
	
	/*! Add an helix to the geometry,

	\remark The helix is a versatile m_geometry->element. For example, to generate a spiral printed circuit antenna, use a helix of zero height.	

	All coordinates are in meters.
		
		\param tag_id The tag ID.
		\param segment_count The number of segments.
		\param s The turn spacing.
		\param h1 The total length of the helix (negative for a left-handed helix).
		
		\param a1 x-start radius.
		\param b1 y-start radius.
		
		\param a2 x-end radius.
		\param b2 y-end radius.
		
		\param rad The wire radius.
	*/		
	void helix( nec_float s, nec_float hl, nec_float a1, nec_float b1,
			nec_float a2, nec_float b2, nec_float rad, int segment_count, int tag_id );
	
		
	
	/*! Scale all dimensions of a structure by a constant.
	
		\param xw1 All structure dimensions, including wire radius, are multiplied by xw1.
	*/
	void scale( nec_float xw1);
	
	
	%extend{
	
		/*! Move the structure with respect to its coordinate system or reproduces structure in new positions,
	
		All coordinates are in meters and angles are in degrees.
		
			\param rox_deg The angle in degrees through which the structure is rotated about the X-axis. A positive angle causes a right-hand rotation.  
			\param roy_deg The angle of rotation about Y-axis.
			\param roz_deg The angle of rotation about Z-axis.
			
			\param xs The x component of vector by which the structure is translated with respect to the coordinate system.
			\param ys The y component of vector by which the structure is translated.
			\param zs The z component of vector by which the structure is translated.
			
			\param its The tag number of the segments that will be moved. If its = 0 then the entire structure is moved.
			\param nrpt The number of new Structures to be generated.
			\param itgi The tag number increment. 
		*/	
		void move( nec_float rox_deg, nec_float roy_deg, nec_float roz_deg, nec_float xs,
			nec_float ys, nec_float zs, int its, int nrpt, int itgi )
		{
			nec_float rox_rad = degrees_to_rad(rox_deg);
			nec_float roy_rad = degrees_to_rad(roy_deg);
			nec_float roz_rad = degrees_to_rad(roz_deg);
			
			return self->move( rox_rad, roy_rad, roz_rad, xs, ys, zs, its, nrpt, itgi );
		}
			



		/*! Reflects partial structure along x,y, or z axes.
			
			\param ix If ix = 1 then the structure is reflected along X axis. 
			\param iy If iy = 1 then the structure is reflected along Y axis.
			\param iz If iz = 1 then the structure is reflected along Z axis.
			\param itx The tag number increment.	
		*/
		void reflect(int ix, int iy, int iz, int itx)
		{
			int nop = 100*ix + 10*iy + iz;
			return self->reflect(ix, iy, iz, itx, nop);
		}
		
		
		
		/*! Rotates structure along Z-axis to complete a symmetric structure.
		
			\param itx The tag number increment.
			\param nop The total number of times that the structure is to occur in the cylindrical array.  
		*/
		void generate_cylindrical_structure(int itx, int nop)
		{
			return self->reflect(-1, 0, 0, itx, nop);
		}
		
		
		
		/*! Add a arbitrary-shaped patch to the geometry
		
		All coordinates are in meters, angles are in degrees.
		
			\param ax1 The x-coordinate of patch center.
			\param ay1 The y-coordinate of patch center.
			\param az1 The z-coordinate of patch center.
			
			\param ax2_deg The elevation angle above the X-Y plane of outward normal vector.
			\param ay2_deg The azimuth angle from X-axis of outward normal vector.
			
			\param az2 The patch area if ny=1.			
		*/
		void arbitrary_shaped_patch( nec_float ax1, nec_float ay1, nec_float az1,
						nec_float ax2_deg, nec_float ay2_deg, nec_float az2 )
		{
			nec_float ax2_rad = degrees_to_rad(ax2_deg);
			nec_float ay2_rad = degrees_to_rad(ay2_deg);
			
			return self->patch( 0, 1, ax1, ay1, az1, ax2_rad, ay2_rad, az2, 0, 0, 0, 0, 0, 0 );
		}
		
		
		
		/*! Add a rectangular patch to the geometry
		
		All coordinates are in meters.		
			
			\param ax1 The x-coordinate of corner 1.
			\param ay1 The y-coordinate of corner 1.
			\param az1 The z-coordinate of corner 1.
			
			\param ax2 The x_coordinate of corner 2.
			\param ay2 The y_coordinate of corner 2.
			\param az2 The z-coordinate of corner 2.
			
			\param ax3 The x_coordinate of corner 3.
			\param ay3 The y_coordinate of corner 3.
			\param az3 The z_coordinate of corner 3.			
		*/
		void rectangular_patch( nec_float ax1, nec_float ay1, nec_float az1,
						nec_float ax2, nec_float ay2, nec_float az2,
						nec_float ax3, nec_float ay3, nec_float az3 )
		{
			return self->patch( 0, 2, ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, 0, 0, 0 );
		}
		
		
		
		/*! Add a triangular patch to the geometry
		
		All coordinates are in meters.		
			
			\param ax1 The x-coordinate of corner 1.
			\param ay1 The y-coordinate of corner 1.
			\param az1 The z-coordinate of corner 1.
			
			\param ax2 The x_coordinate of corner 2.
			\param ay2 The y_coordinate of corner 2.
			\param az2 The z-coordinate of corner 2.
			
			\param ax3 The x_coordinate of corner 3.
			\param ay3 The y_coordinate of corner 3.
			\param az3 The z_coordinate of corner 3.			
		*/
		void triangular_patch( nec_float ax1, nec_float ay1, nec_float az1,
						nec_float ax2, nec_float ay2, nec_float az2,
						nec_float ax3, nec_float ay3, nec_float az3 )
		{
			return self->patch( 0, 3, ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, 0, 0, 0 );
		}
		
		
		
		/*! Add a quadrilateral patch to the geometry
		
		All coordinates are in meters.		
			
			\param ax1 The x-coordinate of corner 1.
			\param ay1 The y-coordinate of corner 1.
			\param az1 The z-coordinate of corner 1.
			
			\param ax2 The x_coordinate of corner 2.
			\param ay2 The y_coordinate of corner 2.
			\param az2 The z-coordinate of corner 2.
			
			\param ax3 The x_coordinate of corner 3.
			\param ay3 The y_coordinate of corner 3.
			\param az3 The z_coordinate of corner 3.
			
			\param ax4 The x_coordinate of corner 4.
			\param ay4 The x_coordinate of corner 4.
			\param az4 The x_coordinate of corner 4.
		*/
		void quadrilateral_patch( nec_float ax1, nec_float ay1, nec_float az1,
						nec_float ax2, nec_float ay2, nec_float az2,
						nec_float ax3, nec_float ay3, nec_float az3,
						nec_float ax4, nec_float ay4, nec_float az4 )
		{
			return self->patch( 0, 4, ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4 );
		}
		
				
		
		/*! Add a multiple patch to the geometry.
		
		All coordinates are in meters.
		
			\param nx The rectangular surface is divided into nx patches from corner 1 to corner 2.
			\param ny The rectangular surface is divided into ny patches from corner 2 to corner 3.
			\param ax1 The x-coordinate of corner 1.
			\param ay1 The y-coordinate of corner 1.
			\param az1 The z-coordinate of corner 1.
			
			\param ax2 The x_coordinate of corner 2.
			\param ay2 The y_coordinate of corner 2.
			\param az2 The z-coordinate of corner 2.
			
			\param ax3 The x_coordinate of corner 3.
			\param ay3 The y_coordinate of corner 3.
			\param az3 The z_coordinate of corner 3.
			
			\param ax4 The x_coordinate of corner 4.
			\param ay4 The x_coordinate of corner 4.
			\param az4 The x_coordinate of corner 4.
		*/
		void multiple_patch( int nx, int ny,
			nec_float ax1, nec_float ay1, nec_float az1,
			nec_float ax2, nec_float ay2, nec_float az2,
			nec_float ax3, nec_float ay3, nec_float az3,
			nec_float ax4, nec_float ay4, nec_float az4 )
		{
			return self->patch( nx, ny, ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4 );
		}		
	}
				
};
