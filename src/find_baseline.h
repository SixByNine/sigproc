float find_rms(float * dat, int ndat, float * ave);

void find_baseline(float * dat, int ndat, double tsamp_dat, double tsmooth,
		   float threshold);
