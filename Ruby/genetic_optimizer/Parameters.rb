SEGMENTS_PER_WAVELENGTH = 31
FREQUENCY=1575.422 #MHz
BANDWIDTH=10 #MHZ
WAVELENGTH = 3.0e8 / (FREQUENCY * 1.0e6)

ANT_HEIGHT=0.3 #meters
ANT_MAX_X=0.2
ANT_MAX_Y=0.2

THETA_MIN=0.0 #Degrees
THETA_MAX=45.0 #Degrees
GAIN_MEAN_LIMIT=6 #dBi (adjust)

# Used to generate random offsets
RANDOM_RANGE=0.1 # Units are Meters

WIRE_CONDUCTIVITY = 3.72e7 # Copper

# Properties of the animal tissue http://niremf.ifac.cnr.it/tissprop/
GROUND_DIELECTRIC_CONSTANT = 39.27 # relative dielectric constant (dry skin)
GROUND_CONDUCTIVITY = 1.0e-2 # conductivity (mhos per meter) (dry skin)
