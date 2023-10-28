void accel(int nsp, double dx, double dt, int t, std::vector<double> q, std::vector<double> m, double& ael, std::vector<double>& a, int ng, std::vector<int> N, std::vector <std::vector<double>>& x, std::vector <std::vector<double>>& vx) {
	/// ACCEL - Calculate force and advance velocity

	const double dxdt = dx / dt;
	for (int species = 0; species < nsp; species++) {

		if (t == 0)
			q[species] = -0.5 * q[species];

		const double ae = (q[species] / m[species]) * (dt / dxdt);

		//  renormalizes acceleration if need be.
		if (ae != ael) {
			const double tem = ae / ael;
			for (int j = 0; j <= ng; j++) {
				a[j] *= tem;
			}
			ael = ae;
		}

		for (int i = 0; i < N[species]; ++i) {
			const int j = floor(x[species][i]);
			const double vo = vx[species][i];
			const double vn = vo + a[j] + (x[species][i] - j) * (a[j + 1] - a[j]);
			vx[species][i] = vn;
		}
	}
}