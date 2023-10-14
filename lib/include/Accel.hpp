void accel(int nsp, double dx, double dt, int t, std::vector<double> q, std::vector<double> m, double& ael, std::vector<double>& a, int ng, std::vector<int> N, std::vector <std::vector<double>>& x, std::vector <std::vector<double>>& vx) {
	/// ACCEL - Calculate force and advance velocity

	for (int species = 0; species < nsp; species++) {

		if (t == 0)
			q[species] = -0.5 * q[species];

		double dxdt = dx / dt;
		double ae = (q[species] / m[species]) * (dt / dxdt);

		//  renormalizes acceleration if need be.

		if (ae != ael) {
			double tem = ae / ael;
			for (int j = 0; j < ng + 1; j++) {
				a[j] *= tem;
			}
			ael = ae;
		}

		double v1s = 0;
		double v2s = 0;
		for (int i = 0; i < N[species]; ++i) {
			int j = floor(x[species][i]);
			double vo = vx[species][i];
			double vn = vo + a[j] + (x[species][i] - j) * (a[j + 1] - a[j]);
			v1s += vn;
			v2s += vo * vn;
			vx[species][i] = vn;
		}
	}
}