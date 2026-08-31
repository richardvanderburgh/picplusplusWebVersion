#include "ChartPanel.h"

#include "ParticleView3D.h"

#include <QChart>
#include <QChartView>
#include <QGridLayout>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QPen>
#include <QScatterSeries>
#include <QSplitter>
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

QColor trailColor(const QColor& base) {
	return QColor(base.red(), base.green(), base.blue(), 90);
}

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

void setupFieldChart(QChart* chart, QLineSeries* series, QValueAxis*& axisX, QValueAxis*& axisY) {
	chart->setAnimationOptions(QChart::NoAnimation);
	axisX = new QValueAxis();
	axisY = new QValueAxis();
	axisX->setTitleText(QStringLiteral("position x"));
	axisY->setTitleText(QStringLiteral("E(x)"));
	chart->addSeries(series);
	chart->addAxis(axisX, Qt::AlignBottom);
	chart->addAxis(axisY, Qt::AlignLeft);
	series->attachAxis(axisX);
	series->attachAxis(axisY);
}

QLineSeries* makeTimeCursor(QChart* chart, QAbstractAxis* axisX, QAbstractAxis* axisY, const QColor& color) {
	auto* cursor = new QLineSeries();
	cursor->setName(QStringLiteral("current time"));
	cursor->setPen(QPen(color, 2, Qt::DashLine));
	chart->addSeries(cursor);
	cursor->attachAxis(axisX);
	cursor->attachAxis(axisY);
	return cursor;
}

double markerSizeForCount(int particleCount) {
	if (particleCount > 4000) {
		return 2.5;
	}
	if (particleCount > 1500) {
		return 3.5;
	}
	if (particleCount > 500) {
		return 5.0;
	}
	return 6.5;
}

template <typename XAccessor, typename YAccessor>
void fillScatterSeries(QScatterSeries* series, QScatterSeries* trailSeries, int species,
	const DATA_STRUCTS::Frame& frame, const DATA_STRUCTS::Frame* previousFrame,
	XAccessor xAccessor, YAccessor yAccessor) {
	QVector<QPointF> points;
	points.reserve(static_cast<int>(frame.particles.size()));
	for (const auto& particle : frame.particles) {
		if (particle.species == species) {
			points.append(QPointF(xAccessor(particle), yAccessor(particle)));
		}
	}
	series->replace(points);

	QVector<QPointF> trailPoints;
	if (previousFrame != nullptr) {
		trailPoints.reserve(static_cast<int>(previousFrame->particles.size()));
		for (const auto& particle : previousFrame->particles) {
			if (particle.species == species) {
				trailPoints.append(QPointF(xAccessor(particle), yAccessor(particle)));
			}
		}
	}
	trailSeries->replace(trailPoints);
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

	m_particleView3D = new ParticleView3D();
	m_view3DTabs = new QTabWidget();
	m_view3DTabs->addTab(m_particleView3D, QStringLiteral("3D particles"));

	auto* projectionTabs = new QTabWidget();
	m_projXYChart = new QChart();
	m_projXYChart->setTitle(QStringLiteral("Particle projection XY"));
	setupScatterChart(m_projXYChart, m_projXYAxisX, m_projXYAxisY, QStringLiteral("x"), QStringLiteral("y"));
	m_projXYView = makeChartView(m_projXYChart);
	projectionTabs->addTab(m_projXYView, QStringLiteral("XY"));

	m_projXZChart = new QChart();
	m_projXZChart->setTitle(QStringLiteral("Particle projection XZ"));
	setupScatterChart(m_projXZChart, m_projXZAxisX, m_projXZAxisY, QStringLiteral("x"), QStringLiteral("z"));
	m_projXZView = makeChartView(m_projXZChart);
	projectionTabs->addTab(m_projXZView, QStringLiteral("XZ"));

	m_projYZChart = new QChart();
	m_projYZChart->setTitle(QStringLiteral("Particle projection YZ"));
	setupScatterChart(m_projYZChart, m_projYZAxisX, m_projYZAxisY, QStringLiteral("y"), QStringLiteral("z"));
	m_projYZView = makeChartView(m_projYZChart);
	projectionTabs->addTab(m_projYZView, QStringLiteral("YZ"));

	m_view3DTabs->addTab(projectionTabs, QStringLiteral("2D slices"));
	m_phaseStack->addWidget(m_view3DTabs);

	m_fieldChart = new QChart();
	m_fieldChart->setTitle(QStringLiteral("Electric field E(x)"));
	m_fieldSeries = new QLineSeries();
	m_fieldSeries->setColor(QColor("#7c3aed"));
	m_fieldSeries->setName(QStringLiteral("E(x)"));
	setupFieldChart(m_fieldChart, m_fieldSeries, m_fieldAxisX, m_fieldAxisY);
	m_fieldView = makeChartView(m_fieldChart);

	auto* animationSplitter = new QSplitter(Qt::Vertical);
	animationSplitter->addWidget(m_phaseStack);
	animationSplitter->addWidget(m_fieldView);
	animationSplitter->setStretchFactor(0, 3);
	animationSplitter->setStretchFactor(1, 2);
	animationSplitter->setSizes({420, 280});
	m_tabs->addTab(animationSplitter, QStringLiteral("Animation"));

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
	m_is3D = false;
	m_currentFrameIndex = -1;

	if (m_particleView3D != nullptr) {
		m_particleView3D->clearView();
	}

	auto clearScatter = [](std::vector<QScatterSeries*>& current, std::vector<QScatterSeries*>& trail, QChart* chart) {
		for (auto* series : current) {
			chart->removeSeries(series);
			delete series;
		}
		current.clear();
		for (auto* series : trail) {
			chart->removeSeries(series);
			delete series;
		}
		trail.clear();
	};

	clearScatter(m_phaseSeries, m_phaseTrailSeries, m_phaseChart);
	clearScatter(m_projXYSeries, m_projXYTrailSeries, m_projXYChart);
	clearScatter(m_projXZSeries, m_projXZTrailSeries, m_projXZChart);
	clearScatter(m_projYZSeries, m_projYZTrailSeries, m_projYZChart);

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
	if (m_energyTimeCursor != nullptr) {
		m_energyChart->removeSeries(m_energyTimeCursor);
		delete m_energyTimeCursor;
		m_energyTimeCursor = nullptr;
	}

	m_modeSeries->clear();
	if (m_modeTimeCursor != nullptr) {
		m_modeChart->removeSeries(m_modeTimeCursor);
		delete m_modeTimeCursor;
		m_modeTimeCursor = nullptr;
	}

	for (auto* axis : m_energyChart->axes()) {
		m_energyChart->removeAxis(axis);
		delete axis;
	}
	m_energyAxisX = nullptr;
	m_energyAxisY = nullptr;

	for (auto* axis : m_modeChart->axes()) {
		m_modeChart->removeAxis(axis);
		delete axis;
	}
	m_modeAxisX = nullptr;
	m_modeAxisY = nullptr;
}

void ChartPanel::setSimulationContext(const FormParams& params) {
	m_params = params;
	m_is3D = params.dimension == 3;
}

double ChartPanel::frameTime(int frameIndex) const {
	if (m_frames.empty() || frameIndex < 0 || frameIndex >= static_cast<int>(m_frames.size())) {
		return 0.0;
	}
	return m_frames[static_cast<size_t>(frameIndex)].frameNumber * m_params.timeStepSize;
}

void ChartPanel::computeFieldBounds() {
	m_fieldYMin = std::numeric_limits<double>::infinity();
	m_fieldYMax = -std::numeric_limits<double>::infinity();
	for (const auto& frame : m_frames) {
		for (const double value : frame.electricField) {
			m_fieldYMin = std::min(m_fieldYMin, value);
			m_fieldYMax = std::max(m_fieldYMax, value);
		}
	}
	if (!std::isfinite(m_fieldYMin) || !std::isfinite(m_fieldYMax)) {
		m_fieldYMin = -1.0;
		m_fieldYMax = 1.0;
	} else if (m_fieldYMin == m_fieldYMax) {
		m_fieldYMin -= 1.0;
		m_fieldYMax += 1.0;
	} else {
		const double pad = std::max((m_fieldYMax - m_fieldYMin) * 0.08, 1e-6);
		m_fieldYMin -= pad;
		m_fieldYMax += pad;
	}
}

void ChartPanel::prepareAnimationSeries() {
	const double markerSize = markerSizeForCount(m_params.numParticles);

	auto createSpeciesSeries = [&](QChart* chart, QValueAxis* axisX, QValueAxis* axisY,
		std::vector<QScatterSeries*>& current, std::vector<QScatterSeries*>& trail) {
		for (int species = 0; species < m_params.numSpecies; ++species) {
			auto* series = new QScatterSeries();
			series->setName(QString("Species %1").arg(species));
			series->setColor(SPECIES_COLORS[species % 4]);
			series->setMarkerSize(markerSize);
			series->setBorderColor(series->color());
			chart->addSeries(series);
			series->attachAxis(axisX);
			series->attachAxis(axisY);
			current.push_back(series);

			auto* trailSeries = new QScatterSeries();
			trailSeries->setName(QString("Species %1 (prev)").arg(species));
			trailSeries->setColor(trailColor(SPECIES_COLORS[species % 4]));
			trailSeries->setMarkerSize(markerSize * 0.75);
			trailSeries->setBorderColor(trailSeries->color());
			chart->addSeries(trailSeries);
			trailSeries->attachAxis(axisX);
			trailSeries->attachAxis(axisY);
			trail.push_back(trailSeries);
		}
	};

	createSpeciesSeries(m_phaseChart, m_phaseXAxis, m_phaseYAxis, m_phaseSeries, m_phaseTrailSeries);
	createSpeciesSeries(m_projXYChart, m_projXYAxisX, m_projXYAxisY, m_projXYSeries, m_projXYTrailSeries);
	createSpeciesSeries(m_projXZChart, m_projXZAxisX, m_projXZAxisY, m_projXZSeries, m_projXZTrailSeries);
	createSpeciesSeries(m_projYZChart, m_projYZAxisX, m_projYZAxisY, m_projYZSeries, m_projYZTrailSeries);

	if (m_is3D && m_particleView3D != nullptr) {
		m_particleView3D->configureDomain(
			m_params.spatialLength,
			m_params.spatialLengthY,
			m_params.spatialLengthZ,
			m_params.numSpecies,
			m_params.numParticles);
	}

	m_phaseXAxis->setRange(0.0, m_params.spatialLength);
	m_phaseYAxis->setRange(m_phaseVelocityMin, m_phaseVelocityMax);
	m_projXYAxisX->setRange(0.0, m_params.spatialLength);
	m_projXYAxisY->setRange(0.0, m_params.spatialLengthY);
	m_projXZAxisX->setRange(0.0, m_params.spatialLength);
	m_projXZAxisY->setRange(0.0, m_params.spatialLengthZ);
	m_projYZAxisX->setRange(0.0, m_params.spatialLengthY);
	m_projYZAxisY->setRange(0.0, m_params.spatialLengthZ);
	m_fieldAxisX->setRange(0.0, m_params.spatialLength);
	m_fieldAxisY->setRange(m_fieldYMin, m_fieldYMax);
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

	computeFieldBounds();
	m_phaseStack->setCurrentIndex(m_is3D ? 1 : 0);
	if (m_is3D && m_view3DTabs != nullptr) {
		m_view3DTabs->setCurrentIndex(0);
	}
	m_tabs->setTabVisible(m_tabs->indexOf(m_modeView), !m_is3D);

	renderEnergyPlot(result);
	if (!m_is3D) {
		renderModePlot(result);
	}

	prepareAnimationSeries();

	if (!m_frames.empty()) {
		m_tabs->setCurrentIndex(0);
		showFrame(0);
	}
}

void ChartPanel::showFrame(int frameIndex) {
	if (m_frames.empty()) {
		return;
	}
	const int clamped = std::clamp(frameIndex, 0, static_cast<int>(m_frames.size()) - 1);
	const auto& frame = m_frames[static_cast<size_t>(clamped)];
	const DATA_STRUCTS::Frame* previousFrame = nullptr;
	if (clamped > 0) {
		previousFrame = &m_frames[static_cast<size_t>(clamped - 1)];
	}

	if (m_is3D) {
		updateParticleView3D(frame, previousFrame);
		updateProjectionPlots(frame, previousFrame);
	} else {
		updatePhasePlot(frame, previousFrame);
	}
	updateFieldPlot(frame);

	const double time = frameTime(clamped);
	updateTimeCursors(time);
	updateChartTitles(time);
	m_currentFrameIndex = clamped;
}

void ChartPanel::updateParticleView3D(const DATA_STRUCTS::Frame& frame, const DATA_STRUCTS::Frame* previousFrame) {
	if (m_particleView3D != nullptr) {
		m_particleView3D->setFrame(frame, previousFrame);
	}
}

void ChartPanel::updatePhasePlot(const DATA_STRUCTS::Frame& frame, const DATA_STRUCTS::Frame* previousFrame) {
	for (int species = 0; species < m_params.numSpecies; ++species) {
		fillScatterSeries(
			m_phaseSeries[static_cast<size_t>(species)],
			m_phaseTrailSeries[static_cast<size_t>(species)],
			species, frame, previousFrame,
			[](const DATA_STRUCTS::Particle& p) { return p.position; },
			[](const DATA_STRUCTS::Particle& p) { return p.velocity; });
	}
}

void ChartPanel::updateProjectionPlots(const DATA_STRUCTS::Frame& frame, const DATA_STRUCTS::Frame* previousFrame) {
	auto updateProjection = [&](QChart* chart, std::vector<QScatterSeries*>& seriesList,
		std::vector<QScatterSeries*>& trailList, auto xAccessor, auto yAccessor) {
		(void)chart;
		for (int species = 0; species < m_params.numSpecies; ++species) {
			fillScatterSeries(
				seriesList[static_cast<size_t>(species)],
				trailList[static_cast<size_t>(species)],
				species, frame, previousFrame, xAccessor, yAccessor);
		}
	};

	updateProjection(
		m_projXYChart, m_projXYSeries, m_projXYTrailSeries,
		[](const DATA_STRUCTS::Particle& p) { return p.position; },
		[](const DATA_STRUCTS::Particle& p) { return p.positionY; });

	updateProjection(
		m_projXZChart, m_projXZSeries, m_projXZTrailSeries,
		[](const DATA_STRUCTS::Particle& p) { return p.position; },
		[](const DATA_STRUCTS::Particle& p) { return p.positionZ; });

	updateProjection(
		m_projYZChart, m_projYZSeries, m_projYZTrailSeries,
		[](const DATA_STRUCTS::Particle& p) { return p.positionY; },
		[](const DATA_STRUCTS::Particle& p) { return p.positionZ; });
}

void ChartPanel::updateFieldPlot(const DATA_STRUCTS::Frame& frame) {
	QVector<QPointF> points;
	const int numGrid = m_params.numGrid;
	const double dx = m_params.spatialLength / numGrid;
	const int fieldSize = static_cast<int>(frame.electricField.size());
	const int plotPoints = m_is3D ? numGrid : numGrid + 1;
	points.reserve(plotPoints);
	for (int j = 0; j < plotPoints && j < fieldSize; ++j) {
		points.append(QPointF(j * dx, frame.electricField[static_cast<size_t>(j)]));
	}
	m_fieldSeries->replace(points);
}

void ChartPanel::updateChartTitles(double time) {
	const QString timeSuffix = QStringLiteral(" · t = %1").arg(time, 0, 'f', 2);
	if (m_is3D) {
		m_view3DTabs->setTabText(0, QStringLiteral("3D particles") + timeSuffix);
		m_projXYChart->setTitle(QStringLiteral("Particle projection XY") + timeSuffix);
		m_projXZChart->setTitle(QStringLiteral("Particle projection XZ") + timeSuffix);
		m_projYZChart->setTitle(QStringLiteral("Particle projection YZ") + timeSuffix);
	} else {
		m_phaseChart->setTitle(QStringLiteral("Phase space (x–v)") + timeSuffix);
	}
	m_fieldChart->setTitle((m_is3D
		? QStringLiteral("Electric field Ex (mid-plane slice)")
		: QStringLiteral("Electric field E(x)")) + timeSuffix);
}

void ChartPanel::updateTimeCursors(double time) {
	if (m_energyTimeCursor != nullptr && m_energyAxisY != nullptr) {
		m_energyTimeCursor->clear();
		m_energyTimeCursor->append(time, 0.0);
		m_energyTimeCursor->append(time, m_energyYMax);
	}

	if (m_modeTimeCursor != nullptr && m_modeAxisY != nullptr && !m_is3D) {
		double yMin = 0.0;
		double yMax = 1.0;
		if (auto* valueAxis = qobject_cast<QValueAxis*>(m_modeAxisY)) {
			yMin = valueAxis->min();
			yMax = valueAxis->max();
		} else if (auto* logAxis = qobject_cast<QLogValueAxis*>(m_modeAxisY)) {
			yMin = logAxis->min();
			yMax = logAxis->max();
		}
		m_modeTimeCursor->clear();
		m_modeTimeCursor->append(time, yMin);
		m_modeTimeCursor->append(time, yMax);
	}
}

void ChartPanel::renderEnergyPlot(const nlohmann::json& result) {
	const auto ke = result.at("ke").get<std::vector<std::vector<double>>>();
	const auto ese = result.at("ese").get<std::vector<double>>();
	const double dt = m_params.timeStepSize;

	m_energyYMax = 0.0;
	for (size_t species = 0; species < ke.size(); ++species) {
		auto* series = new QLineSeries();
		series->setName(QString("KE species %1").arg(static_cast<int>(species)));
		series->setColor(SPECIES_COLORS[species % 4]);
		for (size_t i = 0; i < ke[species].size(); ++i) {
			const double value = ke[species][i];
			series->append(static_cast<double>(i) * dt, value);
			m_energyYMax = std::max(m_energyYMax, value);
		}
		m_energyChart->addSeries(series);
		m_keSeries.push_back(series);
	}

	m_eseSeries = new QLineSeries();
	m_eseSeries->setName(QStringLiteral("Field energy"));
	m_eseSeries->setColor(QColor("#7c3aed"));
	for (size_t i = 0; i < ese.size(); ++i) {
		m_eseSeries->append(static_cast<double>(i) * dt, ese[i]);
		m_energyYMax = std::max(m_energyYMax, ese[i]);
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
		const double total = totalKe + ese[i];
		m_totalSeries->append(static_cast<double>(i) * dt, total);
		m_energyYMax = std::max(m_energyYMax, total);
	}
	m_energyChart->addSeries(m_totalSeries);
	m_energyChart->createDefaultAxes();

	for (auto* axis : m_energyChart->axes()) {
		if (axis->alignment() == Qt::AlignBottom) {
			m_energyAxisX = qobject_cast<QValueAxis*>(axis);
		} else if (axis->alignment() == Qt::AlignLeft) {
			m_energyAxisY = qobject_cast<QValueAxis*>(axis);
		}
	}
	if (m_energyYMax <= 0.0) {
		m_energyYMax = 1.0;
	} else {
		m_energyYMax *= 1.05;
	}

	if (m_energyAxisX != nullptr && m_energyAxisY != nullptr) {
		m_energyTimeCursor = makeTimeCursor(m_energyChart, m_energyAxisX, m_energyAxisY, QColor("#ef4444"));
	}
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
	m_modeAxisX = axisX;

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
		m_modeAxisY = axisY;
	} else {
		auto* axisY = new QValueAxis();
		axisY->setTitleText(QStringLiteral("|E_k|"));
		m_modeChart->addAxis(axisY, Qt::AlignLeft);
		m_modeSeries->attachAxis(axisY);
		m_modeAxisY = axisY;
	}

	m_modeTimeCursor = makeTimeCursor(m_modeChart, m_modeAxisX, m_modeAxisY, QColor("#ef4444"));
	m_modeChart->setTitle(QString("Fourier amplitude |E_k| (mode %1)").arg(mode));
}
