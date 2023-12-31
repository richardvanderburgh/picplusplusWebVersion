#ifndef UTILS_HPP
#define UTILS_HPP


inline std::vector<double> linspace(double start, double end, int num_points) {
	std::vector<double> result(num_points);
	const double step = (end - start) / (num_points - 1);
	for (int i = 0; i < num_points; ++i) {
		result[i] = start + i * step;
	}
	return result;
}
#endif // !UTILS_HPP
