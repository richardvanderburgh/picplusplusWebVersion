#include "ChartPanel.h"

#include <QChart>
#include <QChartView>
#include <QGridLayout>
#include <QLabel>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QScatterSeries>
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

	m_phaseChart = new QChart();
	m_phaseChart->setTitle(QStringLiteral("Phase space"));
	m_phaseXAxis = new QValueAxis();
	m_phaseYAxis = new QValueAxis();
	m_phaseXAxis->setTitleText(QStringLiteral("position x"));
	m_phaseYAxis->setTitleText(QStringLiteral("velocity v"));
	m_phaseChart->addAxis(m_phaseXAxis, Qt::AlignBottom);
	m_phaseChart->addAxis(m_phaseYAxis, Qt::AlignLeft);
	m_phaseView = makeChartView(m_phaseChart);
	m_tabs->addTab(m_phaseView, QStringLiteral("Phase space"));

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

void ChartPanel::clearResults() {
	m_result = nlohmann::json();
	m_frames.clear();
	m_modeAmplitudes.clear();
	m_modeTimes.clear();

	for (auto* series : m_phaseSeries) {
		m_phaseChart->removeSeries(series);
		delete series;
	}
	m_phaseSeries.clear();

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
}

void ChartPanel::renderResults(const nlohmann::json& result) {
	clearResults();
	m_result = result;

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

	renderEnergyPlot(result);
	renderModePlot(result);

	if (!m_frames.empty()) {
		showFrame(0);
	}
}

void ChartPanel::showFrame(int frameIndex) {
	if (m_frames.empty()) {
		return;
	}
	const int clamped = std::clamp(frameIndex, 0, static_cast<int>(m_frames.size()) - 1);
	renderPhasePlot(m_frames[static_cast<size_t>(clamped)]);
	renderFieldPlot(m_frames[static_cast<size_t>(clamped)]);
}

void ChartPanel::renderPhasePlot(const DATA_STRUCTS::Frame& frame) {
	for (auto* series : m_phaseSeries) {
		m_phaseChart->removeSeries(series);
		delete series;
	}
	m_phaseSeries.clear();

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

void ChartPanel::renderFieldPlot(const DATA_STRUCTS::Frame& frame) {
	m_fieldSeries->clear();
	const int numGrid = m_params.numGrid;
	const double dx = m_params.spatialLength / numGrid;
	for (int j = 0; j <= numGrid && j < static_cast<int>(frame.electricField.size()); ++j) {
		m_fieldSeries->append(j * dx, frame.electricField[static_cast<size_t>(j)]);
	}
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
			series->append(i * dt, ke[species][i]);
		}
		m_energyChart->addSeries(series);
		m_keSeries.push_back(series);
	}

	m_eseSeries = new QLineSeries();
	m_eseSeries->setName(QStringLiteral("Field energy"));
	m_eseSeries->setColor(QColor("#7c3aed"));
	for (size_t i = 0; i < ese.size(); ++i) {
		m_eseSeries->append(i * dt, ese[i]);
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
		m_totalSeries->append(i * dt, totalKe + ese[i]);
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
