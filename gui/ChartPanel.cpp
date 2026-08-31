#include "ChartPanel.h"

#include <QChart>
#include <QChartView>
#include <QGridLayout>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QScatterSeries>
#include <QStackedWidget>
#include <QTabWidget>
#include <QValueAxis>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace {

const QColor SPECIES_COLORS[] = {
	QColor("#2563eb"),
	QColor("#dc2626"),
	QColor("#16a34a"),
	QColor("#ea580c"),
};

QChartView* makeChartView(QChart* chart) {
	auto* view = new QChartView(chart);
	view->setRenderHint(QPainter::Antialiasing);
	view->setMinimumHeight(280);
	return view;
}

void setupScatterChart(QChart* chart, QValueAxis*& axisX, QValueAxis*& axisY, const QString& xLabel, const QString& yLabel) {
	chart->setAnimationOptions(QChart::NoAnimation);
	axisX = new QValueAxis();
	axisY = new QValueAxis();
	axisX->setTitleText(xLabel);
	axisY->setTitleText(yLabel);
	chart->addAxis(axisX, Qt::AlignBottom);
	chart->addAxis(axisY, Qt::AlignLeft);
}

} // namespace

ChartPanel::ChartPanel(QWidget* parent)
	: QWidget(parent) {
	setupCharts();

	auto* layout = new QGridLayout(this);
	layout->setContentsMargins(0, 0, 0, 0);
	layout->addWidget(m_tabs, 0, 0);
}

void ChartPanel::setupCharts() {
	m_tabs = new QTabWidget(this);

	m_phaseStack = new QStackedWidget();
	m_phaseChart = new QChart();
	m_phaseChart->setTitle(QStringLiteral("Phase space (x–v)"));
	setupScatterChart(m_phaseChart, m_phaseXAxis, m_phaseYAxis,
		QStringLiteral("position x"), QStringLiteral("velocity v"));
	m_phaseView = makeChartView(m_phaseChart);
	m_phaseStack->addWidget(m_phaseView);

	auto* projectionTabs = new QTabWidget();
	m_projXYChart = new QChart();
	m_projXYChart->setTitle(QStringLiteral("Particle projection XY"));
	QValueAxis *projXYX = nullptr, *projXYY = nullptr;
	setupScatterChart(m_projXYChart, projXYX, projXYY, QStringLiteral("x"), QStringLiteral("y"));
	m_projXYView = makeChartView(m_projXYChart);
	projectionTabs->addTab(m_projXYView, QStringLiteral("XY"));

	m_projXZChart = new QChart();
	m_projXZChart->setTitle(QStringLiteral("Particle projection XZ"));
	QValueAxis *projXZX = nullptr, *projXZY = nullptr;
	setupScatterChart(m_projXZChart, projXZX, projXZY, QStringLiteral("x"), QStringLiteral("z"));
	m_projXZView = makeChartView(m_projXZChart);
	projectionTabs->addTab(m_projXZView, QStringLiteral("XZ"));

	m_projYZChart = new QChart();
	m_projYZChart->setTitle(QStringLiteral("Particle projection YZ"));
	QValueAxis *projYZX = nullptr, *projYZY = nullptr;
	setupScatterChart(m_projYZChart, projYZX, projYZY, QStringLiteral("y"), QStringLiteral("z"));
	m_projYZView = makeChartView(m_projYZChart);
	projectionTabs->addTab(m_projYZView, QStringLiteral("YZ"));

	m_phaseStack->addWidget(projectionTabs);
	m_tabs->addTab(m_phaseStack, QStringLiteral("Particles"));

	m_fieldChart = new QChart();
	m_fieldChart->setTitle(QStringLiteral("Electric field E(x)"));
	m_fieldSeries = new QLineSeries();
	m_fieldSeries->setColor(QColor("#7c3aed"));
	m_fieldSeries->setName(QStringLiteral("E(x)"));
	m_fieldChart->addSeries(m_fieldSeries);
	m_fieldChart->createDefaultAxes();
	m_fieldView = makeChartView(m_fieldChart);
	m_tabs->addTab(m_fieldView, QStringLiteral("E(x)"));

	m_energyChart = new QChart();
	m_energyChart->setTitle(QStringLiteral("Energy budget"));
	m_energyView = makeChartView(m_energyChart);
	m_tabs->addTab(m_energyView, QStringLiteral("Energy"));

	m_modeChart = new QChart();
	m_modeChart->setTitle(QStringLiteral("Fourier amplitude |E_k|"));
	m_modeSeries = new QLineSeries();
	m_modeSeries->setColor(QColor("#b45309"));
	m_modeChart->addSeries(m_modeSeries);
	m_modeChart->createDefaultAxes();
	m_modeView = makeChartView(m_modeChart);
	m_tabs->addTab(m_modeView, QStringLiteral("|E_k|"));
}

void ChartPanel::clearScatterSeries(std::vector<QScatterSeries*>& seriesList, QChart* chart) {
	for (auto* series : seriesList) {
		chart->removeSeries(series);
		delete series;
	}
	seriesList.clear();
}

void ChartPanel::clearResults() {
	m_result = nlohmann::json();
	m_frames.clear();
	m_modeAmplitudes.clear();
	m_modeTimes.clear();
	m_is3D = false;

	clearScatterSeries(m_phaseSeries, m_phaseChart);
	clearScatterSeries(m_projXYSeries, m_projXYChart);
	clearScatterSeries(m_projXZSeries, m_projXZChart);
	clearScatterSeries(m_projYZSeries, m_projYZChart);

	m_fieldSeries->clear();

	for (auto* series : m_keSeries) {
		m_energyChart->removeSeries(series);
		delete series;
	}
	m_keSeries.clear();

	if (m_eseSeries != nullptr) {
		m_energyChart->removeSeries(m_eseSeries);
		delete m_eseSeries;
		m_eseSeries = nullptr;
	}
	if (m_totalSeries != nullptr) {
		m_energyChart->removeSeries(m_totalSeries);
		delete m_totalSeries;
		m_totalSeries = nullptr;
	}

	m_modeSeries->clear();

	for (auto* axis : m_energyChart->axes()) {
		m_energyChart->removeAxis(axis);
		delete axis;
	}
}

void ChartPanel::setSimulationContext(const FormParams& params) {
	m_params = params;
	m_is3D = params.dimension == 3;
}

void ChartPanel::renderResults(const nlohmann::json& result) {
	clearResults();
	m_result = result;
	m_is3D = result.value("dimension", m_params.dimension) == 3;

	if (result.contains("phaseFrames")) {
		m_frames = result["phaseFrames"].get<std::vector<DATA_STRUCTS::Frame>>();
	}

	double vMin = std::numeric_limits<double>::infinity();
	double vMax = -std::numeric_limits<double>::infinity();
	for (const auto& frame : m_frames) {
		for (const auto& particle : frame.particles) {
			vMin = std::min(vMin, particle.velocity);
			vMax = std::max(vMax, particle.velocity);
		}
	}

	const double estimate = std::max((std::abs(m_params.driftVelocity) + 4.0 * std::max(0.0, m_params.thermalVelocity)) * 1.5, 0.5);
	if (!std::isfinite(vMin) || !std::isfinite(vMax)) {
		m_phaseVelocityMin = -estimate;
		m_phaseVelocityMax = estimate;
	} else {
		vMin = std::min(vMin, -estimate);
		vMax = std::max(vMax, estimate);
		const double pad = std::max((vMax - vMin) * 0.15, 0.1);
		m_phaseVelocityMin = vMin - pad;
		m_phaseVelocityMax = vMax + pad;
	}

	m_phaseStack->setCurrentIndex(m_is3D ? 1 : 0);
	m_tabs->setTabVisible(m_tabs->indexOf(m_modeView), !m_is3D);

	renderEnergyPlot(result);
	if (!m_is3D) {
		renderModePlot(result);
	}

	if (!m_frames.empty()) {
		showFrame(0);
	}
}

void ChartPanel::showFrame(int frameIndex) {
	if (m_frames.empty()) {
		return;
	}
	const int clamped = std::clamp(frameIndex, 0, static_cast<int>(m_frames.size()) - 1);
	const auto& frame = m_frames[static_cast<size_t>(clamped)];
	if (m_is3D) {
		renderProjectionPlots(frame);
	} else {
		renderPhasePlot(frame);
	}
	renderFieldPlot(frame);
}

void ChartPanel::renderPhasePlot(const DATA_STRUCTS::Frame& frame) {
	clearScatterSeries(m_phaseSeries, m_phaseChart);

	for (int species = 0; species < m_params.numSpecies; ++species) {
		auto* series = new QScatterSeries();
		series->setName(QString("Species %1").arg(species));
		series->setColor(SPECIES_COLORS[species % 4]);
		series->setMarkerSize(6.0);

		for (const auto& particle : frame.particles) {
			if (particle.species == species) {
				series->append(particle.position, particle.velocity);
			}
		}

		m_phaseChart->addSeries(series);
		series->attachAxis(m_phaseXAxis);
		series->attachAxis(m_phaseYAxis);
		m_phaseSeries.push_back(series);
	}

	m_phaseXAxis->setRange(0.0, m_params.spatialLength);
	m_phaseYAxis->setRange(m_phaseVelocityMin, m_phaseVelocityMax);
}

void ChartPanel::renderProjectionPlots(const DATA_STRUCTS::Frame& frame) {
	auto renderProjection = [&](QChart* chart, std::vector<QScatterSeries*>& seriesList,
		auto xAccessor, auto yAccessor, double xMax, double yMax, const QString& tabName) {
		(void)tabName;
		clearScatterSeries(seriesList, chart);
		QValueAxis* axisX = nullptr;
		QValueAxis* axisY = nullptr;
		for (auto* axis : chart->axes()) {
			if (axis->alignment() == Qt::AlignBottom) {
				axisX = qobject_cast<QValueAxis*>(axis);
			} else if (axis->alignment() == Qt::AlignLeft) {
				axisY = qobject_cast<QValueAxis*>(axis);
			}
		}

		for (int species = 0; species < m_params.numSpecies; ++species) {
			auto* series = new QScatterSeries();
			series->setName(QString("Species %1").arg(species));
			series->setColor(SPECIES_COLORS[species % 4]);
			series->setMarkerSize(5.0);
			for (const auto& particle : frame.particles) {
				if (particle.species == species) {
					series->append(xAccessor(particle), yAccessor(particle));
				}
			}
			chart->addSeries(series);
			if (axisX != nullptr) {
				series->attachAxis(axisX);
			}
			if (axisY != nullptr) {
				series->attachAxis(axisY);
			}
			seriesList.push_back(series);
		}
		if (axisX != nullptr) {
			axisX->setRange(0.0, xMax);
		}
		if (axisY != nullptr) {
			axisY->setRange(0.0, yMax);
		}
	};

	renderProjection(
		m_projXYChart, m_projXYSeries,
		[](const DATA_STRUCTS::Particle& p) { return p.position; },
		[](const DATA_STRUCTS::Particle& p) { return p.positionY; },
		m_params.spatialLength, m_params.spatialLengthY, QStringLiteral("XY"));

	renderProjection(
		m_projXZChart, m_projXZSeries,
		[](const DATA_STRUCTS::Particle& p) { return p.position; },
		[](const DATA_STRUCTS::Particle& p) { return p.positionZ; },
		m_params.spatialLength, m_params.spatialLengthZ, QStringLiteral("XZ"));

	renderProjection(
		m_projYZChart, m_projYZSeries,
		[](const DATA_STRUCTS::Particle& p) { return p.positionY; },
		[](const DATA_STRUCTS::Particle& p) { return p.positionZ; },
		m_params.spatialLengthY, m_params.spatialLengthZ, QStringLiteral("YZ"));
}

void ChartPanel::renderFieldPlot(const DATA_STRUCTS::Frame& frame) {
	m_fieldSeries->clear();
	const int numGrid = m_params.numGrid;
	const double dx = m_params.spatialLength / numGrid;
	const int fieldSize = static_cast<int>(frame.electricField.size());
	const int plotPoints = m_is3D ? numGrid : numGrid + 1;
	for (int j = 0; j < plotPoints && j < fieldSize; ++j) {
		m_fieldSeries->append(j * dx, frame.electricField[static_cast<size_t>(j)]);
	}
	m_fieldChart->setTitle(m_is3D
		? QStringLiteral("Electric field Ex (mid-plane slice)")
		: QStringLiteral("Electric field E(x)"));
}

void ChartPanel::renderEnergyPlot(const nlohmann::json& result) {
	const auto ke = result.at("ke").get<std::vector<std::vector<double>>>();
	const auto ese = result.at("ese").get<std::vector<double>>();
	const double dt = m_params.timeStepSize;

	for (size_t species = 0; species < ke.size(); ++species) {
		auto* series = new QLineSeries();
		series->setName(QString("KE species %1").arg(static_cast<int>(species)));
		series->setColor(SPECIES_COLORS[species % 4]);
		for (size_t i = 0; i < ke[species].size(); ++i) {
			series->append(static_cast<double>(i) * dt, ke[species][i]);
		}
		m_energyChart->addSeries(series);
		m_keSeries.push_back(series);
	}

	m_eseSeries = new QLineSeries();
	m_eseSeries->setName(QStringLiteral("Field energy"));
	m_eseSeries->setColor(QColor("#7c3aed"));
	for (size_t i = 0; i < ese.size(); ++i) {
		m_eseSeries->append(static_cast<double>(i) * dt, ese[i]);
	}
	m_energyChart->addSeries(m_eseSeries);

	m_totalSeries = new QLineSeries();
	m_totalSeries->setName(QStringLiteral("Total energy"));
	m_totalSeries->setColor(QColor("#0f766e"));
	for (size_t i = 0; i < ese.size(); ++i) {
		double totalKe = 0.0;
		for (const auto& speciesKe : ke) {
			if (i < speciesKe.size()) {
				totalKe += speciesKe[i];
			}
		}
		m_totalSeries->append(static_cast<double>(i) * dt, totalKe + ese[i]);
	}
	m_energyChart->addSeries(m_totalSeries);
	m_energyChart->createDefaultAxes();
}

double ChartPanel::fourierModeAmplitude(const std::vector<double>& field, int mode, double domainLength) const {
	const int n = static_cast<int>(field.size()) - 1;
	if (n <= 0) {
		return 0.0;
	}
	const double k = (2.0 * std::numbers::pi * mode) / domainLength;
	double re = 0.0;
	double im = 0.0;
	for (int j = 0; j < n; ++j) {
		const double x = (static_cast<double>(j) / n) * domainLength;
		const double phase = k * x;
		re += field[static_cast<size_t>(j)] * std::cos(phase);
		im += field[static_cast<size_t>(j)] * std::sin(phase);
	}
	re /= n;
	im /= n;
	return std::sqrt(re * re + im * im);
}

void ChartPanel::renderModePlot(const nlohmann::json& result) {
	(void)result;
	m_modeSeries->clear();
	m_modeAmplitudes.clear();
	m_modeTimes.clear();

	const int mode = m_params.spatialPerturbationMode;
	const double dt = m_params.timeStepSize;
	const double domainLength = m_params.spatialLength;

	bool useLog = false;
	for (const auto& frame : m_frames) {
		const double amplitude = fourierModeAmplitude(frame.electricField, mode, domainLength);
		m_modeAmplitudes.push_back(amplitude);
		m_modeTimes.push_back(frame.frameNumber * dt);
		if (amplitude > 0.0) {
			useLog = true;
		}
	}

	for (size_t i = 0; i < m_modeAmplitudes.size(); ++i) {
		m_modeSeries->append(m_modeTimes[i], m_modeAmplitudes[i]);
	}

	for (auto* axis : m_modeChart->axes()) {
		m_modeChart->removeAxis(axis);
		delete axis;
	}

	auto* axisX = new QValueAxis();
	axisX->setTitleText(QStringLiteral("time"));
	m_modeChart->addAxis(axisX, Qt::AlignBottom);
	m_modeSeries->attachAxis(axisX);

	if (useLog && !m_modeAmplitudes.empty()) {
		auto* axisY = new QLogValueAxis();
		axisY->setTitleText(QStringLiteral("|E_k|"));
		const double maxVal = *std::max_element(m_modeAmplitudes.begin(), m_modeAmplitudes.end());
		double minPositive = maxVal;
		for (const double value : m_modeAmplitudes) {
			if (value > 0.0) {
				minPositive = std::min(minPositive, value);
			}
		}
		axisY->setRange(std::max(minPositive * 0.5, 1e-12), maxVal * 1.1);
		m_modeChart->addAxis(axisY, Qt::AlignLeft);
		m_modeSeries->attachAxis(axisY);
	} else {
		auto* axisY = new QValueAxis();
		axisY->setTitleText(QStringLiteral("|E_k|"));
		m_modeChart->addAxis(axisY, Qt::AlignLeft);
		m_modeSeries->attachAxis(axisY);
	}

	m_modeChart->setTitle(QString("Fourier amplitude |E_k| (mode %1)").arg(mode));
}
