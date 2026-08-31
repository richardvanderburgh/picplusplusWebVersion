#include "MainWindow.h"

#include "ChartPanel.h"
#include "SimulationWorker.h"

#include <QCheckBox>
#include <QKeySequence>
#include <QShortcut>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSlider>
#include <QSpinBox>
#include <QSplitter>
#include <QThread>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>

namespace {

QDoubleSpinBox* makeDoubleSpin(double value, double step, int decimals = 6) {
	auto* spin = new QDoubleSpinBox();
	spin->setDecimals(decimals);
	spin->setSingleStep(step);
	spin->setRange(-1e12, 1e12);
	spin->setValue(value);
	return spin;
}

QSpinBox* makeIntSpin(int value, int min, int max) {
	auto* spin = new QSpinBox();
	spin->setRange(min, max);
	spin->setValue(value);
	return spin;
}

double meanSlice(const std::vector<double>& values, size_t begin, size_t end) {
	if (values.empty()) {
		return 0.0;
	}
	begin = std::min(begin, values.size());
	end = std::min(end, values.size());
	if (begin >= end) {
		return 0.0;
	}
	double sum = 0.0;
	for (size_t i = begin; i < end; ++i) {
		sum += values[i];
	}
	return sum / static_cast<double>(end - begin);
}

} // namespace

MainWindow::MainWindow(const std::string& repoRoot, QWidget* parent)
	: QMainWindow(parent), m_repoRoot(repoRoot) {
	setWindowTitle(QStringLiteral("PIC++ Plasma Simulator"));
	resize(1280, 820);

	m_workerThread = new QThread(this);
	m_worker = new SimulationWorker();
	m_worker->moveToThread(m_workerThread);
	connect(m_workerThread, &QThread::finished, m_worker, &QObject::deleteLater);
	connect(this, &MainWindow::destroyed, m_workerThread, &QThread::quit);
	connect(m_worker, &SimulationWorker::progressUpdated, this, &MainWindow::onProgressUpdated);
	connect(m_worker, &SimulationWorker::simulationFinished, this, &MainWindow::onSimulationFinished);
	connect(m_worker, &SimulationWorker::simulationFailed, this, &MainWindow::onSimulationFailed);
	m_workerThread->start();

	m_animationTimer = new QTimer(this);
	connect(m_animationTimer, &QTimer::timeout, this, &MainWindow::onAnimationTick);

	auto* splitter = new QSplitter(Qt::Horizontal, this);
	setCentralWidget(splitter);

	auto* sidebar = new QWidget();
	auto* sidebarLayout = new QVBoxLayout(sidebar);

	auto* header = new QLabel(
		QStringLiteral("<h2>PIC++ Plasma Simulator</h2>"
		               "<p>1D/3D electrostatic particle-in-cell — pick a demo, run, explore.</p>"));
	header->setWordWrap(true);
	sidebarLayout->addWidget(header);

	m_demoSelect = new QComboBox();
	m_demoSelect->addItem(QStringLiteral("— choose a preset —"), QString());
	connect(m_demoSelect, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onDemoChanged);
	sidebarLayout->addWidget(new QLabel(QStringLiteral("Demo case")));
	sidebarLayout->addWidget(m_demoSelect);

	m_demoDescription = new QLabel();
	m_demoDescription->setWordWrap(true);
	m_demoDescription->setStyleSheet(QStringLiteral("color: #64748b;"));
	sidebarLayout->addWidget(m_demoDescription);

	auto* runRow = new QHBoxLayout();
	m_runButton = new QPushButton(QStringLiteral("Run simulation"));
	m_resetButton = new QPushButton(QStringLiteral("Reset"));
	m_resetButton->setEnabled(false);
	connect(m_runButton, &QPushButton::clicked, this, &MainWindow::onRunClicked);
	connect(m_resetButton, &QPushButton::clicked, this, &MainWindow::onResetClicked);
	runRow->addWidget(m_runButton, 1);
	runRow->addWidget(m_resetButton);
	sidebarLayout->addLayout(runRow);

	m_customizedBadge = new QLabel(QStringLiteral("Parameters edited"));
	m_customizedBadge->setStyleSheet(QStringLiteral("color: #b45309; font-weight: 600;"));
	m_customizedBadge->hide();
	sidebarLayout->addWidget(m_customizedBadge);

	m_progressBar = new QProgressBar();
	m_progressBar->setRange(0, 100);
	m_progressBar->hide();
	sidebarLayout->addWidget(m_progressBar);

	m_statusLabel = new QLabel();
	m_statusLabel->setWordWrap(true);
	sidebarLayout->addWidget(m_statusLabel);

	auto* paramsGroup = new QGroupBox(QStringLiteral("Simulation parameters"));
	m_paramsLayout = new QFormLayout(paramsGroup);

	m_spatialLength = makeDoubleSpin(6.28318530717958, 0.1);
	m_numParticles = makeIntSpin(500, 1, 1000000);
	m_timeSteps = makeIntSpin(500, 1, 1000000);
	m_timeStepSize = makeDoubleSpin(0.2, 0.01);
	m_numGrid = makeIntSpin(32, 2, 4096);
	m_numGridY = makeIntSpin(8, 2, 4096);
	m_numGridZ = makeIntSpin(8, 2, 4096);
	m_spatialLengthY = makeDoubleSpin(6.28318530717958, 0.1);
	m_spatialLengthZ = makeDoubleSpin(6.28318530717958, 0.1);
	m_spatialPerturbationMode = makeIntSpin(1, 0, 64);
	m_driftVelocity = makeDoubleSpin(1.0, 0.1);
	m_numSpecies = makeIntSpin(2, 1, 2);
	m_spatialPerturbationAmplitude = makeDoubleSpin(0.001, 0.0001, 8);
	m_thermalVelocity = makeDoubleSpin(0.0, 0.01);
	m_plasmaFrequency = makeDoubleSpin(1.0, 0.1);
	m_chargeMassRatio = makeDoubleSpin(-1.0, 0.1);
	m_framePeriod = makeIntSpin(5, 1, 10000);

	m_spatialPerturbationWaveform = new QComboBox();
	m_spatialPerturbationWaveform->addItem(QStringLiteral("cos — density perturbation"), QStringLiteral("cos"));
	m_spatialPerturbationWaveform->addItem(QStringLiteral("sin — standing wave (Landau)"), QStringLiteral("sin"));

	m_dimension = new QComboBox();
	m_dimension->addItem(QStringLiteral("1D"), 1);
	m_dimension->addItem(QStringLiteral("3D"), 3);
	connect(m_dimension, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onDimensionChanged);

	m_paramsLayout->addRow(QStringLiteral("Dimension"), m_dimension);
	m_paramsLayout->addRow(QStringLiteral("Domain length (x)"), m_spatialLength);
	m_length3DRowY = new QWidget();
	auto* lengthYLayout = new QHBoxLayout(m_length3DRowY);
	lengthYLayout->setContentsMargins(0, 0, 0, 0);
	lengthYLayout->addWidget(m_spatialLengthY);
	m_paramsLayout->addRow(QStringLiteral("Domain length (y)"), m_length3DRowY);
	m_length3DRowZ = new QWidget();
	auto* lengthZLayout = new QHBoxLayout(m_length3DRowZ);
	lengthZLayout->setContentsMargins(0, 0, 0, 0);
	lengthZLayout->addWidget(m_spatialLengthZ);
	m_paramsLayout->addRow(QStringLiteral("Domain length (z)"), m_length3DRowZ);
	m_paramsLayout->addRow(QStringLiteral("Particles / species"), m_numParticles);
	m_paramsLayout->addRow(QStringLiteral("Time steps"), m_timeSteps);
	m_paramsLayout->addRow(QStringLiteral("Δt"), m_timeStepSize);
	m_paramsLayout->addRow(QStringLiteral("Grid cells (x)"), m_numGrid);
	m_grid3DRowY = new QWidget();
	auto* gridYLayout = new QHBoxLayout(m_grid3DRowY);
	gridYLayout->setContentsMargins(0, 0, 0, 0);
	gridYLayout->addWidget(m_numGridY);
	m_paramsLayout->addRow(QStringLiteral("Grid cells (y)"), m_grid3DRowY);
	m_grid3DRowZ = new QWidget();
	auto* gridZLayout = new QHBoxLayout(m_grid3DRowZ);
	gridZLayout->setContentsMargins(0, 0, 0, 0);
	gridZLayout->addWidget(m_numGridZ);
	m_paramsLayout->addRow(QStringLiteral("Grid cells (z)"), m_grid3DRowZ);
	m_paramsLayout->addRow(QStringLiteral("Frame period"), m_framePeriod);
	m_paramsLayout->addRow(QStringLiteral("Drift v"), m_driftVelocity);
	m_paramsLayout->addRow(QStringLiteral("Species"), m_numSpecies);
	m_paramsLayout->addRow(QStringLiteral("Mode m"), m_spatialPerturbationMode);
	m_paramsLayout->addRow(QStringLiteral("Pert. amplitude"), m_spatialPerturbationAmplitude);
	m_paramsLayout->addRow(QStringLiteral("v_th"), m_thermalVelocity);
	m_paramsLayout->addRow(QStringLiteral("ω_p"), m_plasmaFrequency);
	m_paramsLayout->addRow(QStringLiteral("q/m"), m_chargeMassRatio);
	m_paramsLayout->addRow(QStringLiteral("Perturbation waveform"), m_spatialPerturbationWaveform);

	const auto markEdited = [this]() { markParamsCustomized(); };
	for (auto* widget : std::initializer_list<QDoubleSpinBox*>{
		m_spatialLength, m_spatialLengthY, m_spatialLengthZ, m_timeStepSize, m_driftVelocity,
		m_spatialPerturbationAmplitude, m_thermalVelocity, m_plasmaFrequency, m_chargeMassRatio}) {
		connect(widget, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, markEdited);
	}
	for (auto* widget : std::initializer_list<QSpinBox*>{
		m_numParticles, m_timeSteps, m_numGrid, m_numGridY, m_numGridZ, m_spatialPerturbationMode,
		m_numSpecies, m_framePeriod}) {
		connect(widget, QOverload<int>::of(&QSpinBox::valueChanged), this, markEdited);
	}
	connect(m_spatialPerturbationWaveform, QOverload<int>::of(&QComboBox::currentIndexChanged), this, markEdited);
	connect(m_dimension, QOverload<int>::of(&QComboBox::currentIndexChanged), this, markEdited);

	auto* paramsScroll = new QScrollArea();
	paramsScroll->setWidget(paramsGroup);
	paramsScroll->setWidgetResizable(true);
	paramsScroll->setMaximumHeight(360);
	sidebarLayout->addWidget(paramsScroll);

	auto* summaryGroup = new QGroupBox(QStringLiteral("Analysis summary"));
	auto* summaryLayout = new QFormLayout(summaryGroup);
	m_statEseRatio = new QLabel(QStringLiteral("—"));
	m_statEnergyDrift = new QLabel(QStringLiteral("—"));
	m_statEkRatio = new QLabel(QStringLiteral("—"));
	m_statEsePeak = new QLabel(QStringLiteral("—"));
	summaryLayout->addRow(QStringLiteral("Field energy late / early"), m_statEseRatio);
	summaryLayout->addRow(QStringLiteral("Total energy drift"), m_statEnergyDrift);
	summaryLayout->addRow(QStringLiteral("|E_k| end / start"), m_statEkRatio);
	summaryLayout->addRow(QStringLiteral("Peak field / initial"), m_statEsePeak);
	summaryGroup->hide();
	sidebarLayout->addWidget(summaryGroup);
	sidebarLayout->addStretch();

	auto* resultsPane = new QWidget();
	auto* resultsLayout = new QVBoxLayout(resultsPane);

	auto* controls = new QHBoxLayout();
	auto* firstButton = new QPushButton(QStringLiteral("⏮"));
	auto* prevButton = new QPushButton(QStringLiteral("◀"));
	m_playButton = new QPushButton(QStringLiteral("Play"));
	m_pauseButton = new QPushButton(QStringLiteral("Pause"));
	auto* nextButton = new QPushButton(QStringLiteral("▶"));
	auto* lastButton = new QPushButton(QStringLiteral("⏭"));
	m_frameSlider = new QSlider(Qt::Horizontal);
	m_frameSlider->setEnabled(false);
	m_frameLabel = new QLabel(QStringLiteral("—"));
	m_loopPlayback = new QCheckBox(QStringLiteral("Loop"));
	m_loopPlayback->setChecked(true);
	m_animSpeedControl = makeIntSpin(33, 10, 500);
	m_animSpeedControl->setSuffix(QStringLiteral(" ms"));
	m_animSpeedControl->setToolTip(QStringLiteral("Milliseconds between interpolated animation steps"));

	connect(firstButton, &QPushButton::clicked, this, &MainWindow::onFirstFrame);
	connect(prevButton, &QPushButton::clicked, this, &MainWindow::onPrevFrame);
	connect(m_playButton, &QPushButton::clicked, this, &MainWindow::onPlayClicked);
	connect(m_pauseButton, &QPushButton::clicked, this, &MainWindow::onPauseClicked);
	connect(nextButton, &QPushButton::clicked, this, &MainWindow::onNextFrame);
	connect(lastButton, &QPushButton::clicked, this, &MainWindow::onLastFrame);
	connect(m_frameSlider, &QSlider::valueChanged, this, &MainWindow::onFrameSliderChanged);
	connect(m_frameSlider, &QSlider::sliderPressed, this, &MainWindow::onFrameSliderPressed);
	connect(m_frameSlider, &QSlider::sliderReleased, this, &MainWindow::onFrameSliderReleased);
	connect(m_animSpeedControl, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::onAnimSpeedChanged);

	controls->addWidget(new QLabel(QStringLiteral("Explore:")));
	controls->addWidget(firstButton);
	controls->addWidget(prevButton);
	controls->addWidget(m_playButton);
	controls->addWidget(m_pauseButton);
	controls->addWidget(nextButton);
	controls->addWidget(lastButton);
	controls->addWidget(m_frameSlider, 1);
	controls->addWidget(m_loopPlayback);
	controls->addWidget(m_animSpeedControl);
	controls->addWidget(m_frameLabel);
	resultsLayout->addLayout(controls);

	auto* shortcutPlayPause = new QShortcut(QKeySequence(Qt::Key_Space), this);
	connect(shortcutPlayPause, &QShortcut::activated, this, [this]() {
		if (m_animationTimer->isActive()) {
			onPauseClicked();
		} else {
			onPlayClicked();
		}
	});
	auto* shortcutPrev = new QShortcut(QKeySequence(Qt::Key_Left), this);
	connect(shortcutPrev, &QShortcut::activated, this, &MainWindow::onPrevFrame);
	auto* shortcutNext = new QShortcut(QKeySequence(Qt::Key_Right), this);
	connect(shortcutNext, &QShortcut::activated, this, &MainWindow::onNextFrame);

	m_chartPanel = new ChartPanel();
	resultsLayout->addWidget(m_chartPanel, 1);

	splitter->addWidget(sidebar);
	splitter->addWidget(resultsPane);
	splitter->setStretchFactor(0, 0);
	splitter->setStretchFactor(1, 1);
	splitter->setSizes({360, 900});

	populateDemoSelector();

	std::string error;
	if (const auto params = SimulationConfig::formParamsFromDemo(m_repoRoot, "vortexTwoStream", error)) {
		applyFormParams(*params);
		m_savedDemoParams = *params;
		for (int i = 0; i < m_demoSelect->count(); ++i) {
			if (m_demoSelect->itemData(i).toString() == QStringLiteral("vortexTwoStream")) {
				m_demoSelect->setCurrentIndex(i);
				break;
			}
		}
	} else if (!error.empty()) {
		m_statusLabel->setText(QString::fromStdString(error));
	}
	updateDimensionControls();
}

void MainWindow::populateDemoSelector() {
	std::string error;
	m_demos = SimulationConfig::loadDemoManifest(m_repoRoot, error);
	if (!error.empty()) {
		m_statusLabel->setText(QString::fromStdString(error));
		return;
	}

	QString currentCategory;
	for (const auto& demo : m_demos) {
		const QString category = QString::fromStdString(demo.category);
		if (category != currentCategory) {
			currentCategory = category;
		}
		const QString label = category.isEmpty()
			? QString::fromStdString(demo.title)
			: QString("%1 — %2").arg(category, QString::fromStdString(demo.title));
		m_demoSelect->addItem(label, QString::fromStdString(demo.id));
	}
}

void MainWindow::updateDimensionControls() {
	const bool is3D = m_dimension->currentData().toInt() == 3;
	const auto setRowVisible = [this](QWidget* field, bool visible) {
		field->setVisible(visible);
		if (m_paramsLayout != nullptr) {
			if (auto* label = m_paramsLayout->labelForField(field)) {
				label->setVisible(visible);
			}
		}
	};
	setRowVisible(m_length3DRowY, is3D);
	setRowVisible(m_length3DRowZ, is3D);
	setRowVisible(m_grid3DRowY, is3D);
	setRowVisible(m_grid3DRowZ, is3D);
}

void MainWindow::onDimensionChanged(int /*index*/) {
	updateDimensionControls();
}

void MainWindow::applyFormParams(const FormParams& params) {
	const int dimIndex = m_dimension->findData(params.dimension);
	if (dimIndex >= 0) {
		m_dimension->blockSignals(true);
		m_dimension->setCurrentIndex(dimIndex);
		m_dimension->blockSignals(false);
	}
	m_spatialLengthY->setValue(params.spatialLengthY);
	m_spatialLengthZ->setValue(params.spatialLengthZ);
	m_spatialLength->setValue(params.spatialLength);
	m_numParticles->setValue(params.numParticles);
	m_timeSteps->setValue(params.timeSteps);
	m_timeStepSize->setValue(params.timeStepSize);
	m_numGrid->setValue(params.numGrid);
	m_numGridY->setValue(params.numGridY);
	m_numGridZ->setValue(params.numGridZ);
	m_spatialPerturbationMode->setValue(params.spatialPerturbationMode);
	m_driftVelocity->setValue(params.driftVelocity);
	m_numSpecies->setValue(params.numSpecies);
	m_spatialPerturbationAmplitude->setValue(params.spatialPerturbationAmplitude);
	m_thermalVelocity->setValue(params.thermalVelocity);
	m_plasmaFrequency->setValue(params.plasmaFrequency);
	m_chargeMassRatio->setValue(params.chargeMassRatio);
	m_framePeriod->setValue(params.framePeriod);

	const int waveformIndex = m_spatialPerturbationWaveform->findData(
		QString::fromStdString(params.spatialPerturbationWaveform));
	if (waveformIndex >= 0) {
		m_spatialPerturbationWaveform->setCurrentIndex(waveformIndex);
	}
	updateDimensionControls();
}

FormParams MainWindow::collectFormParams() const {
	FormParams params;
	params.dimension = m_dimension->currentData().toInt();
	params.spatialLength = m_spatialLength->value();
	params.spatialLengthY = m_spatialLengthY->value();
	params.spatialLengthZ = m_spatialLengthZ->value();
	params.numParticles = m_numParticles->value();
	params.timeSteps = m_timeSteps->value();
	params.timeStepSize = m_timeStepSize->value();
	params.numGrid = m_numGrid->value();
	params.numGridY = m_numGridY->value();
	params.numGridZ = m_numGridZ->value();
	params.spatialPerturbationMode = m_spatialPerturbationMode->value();
	params.driftVelocity = m_driftVelocity->value();
	params.numSpecies = m_numSpecies->value();
	params.spatialPerturbationAmplitude = m_spatialPerturbationAmplitude->value();
	params.thermalVelocity = m_thermalVelocity->value();
	params.plasmaFrequency = m_plasmaFrequency->value();
	params.chargeMassRatio = m_chargeMassRatio->value();
	params.framePeriod = m_framePeriod->value();
	params.spatialPerturbationWaveform = m_spatialPerturbationWaveform->currentData().toString().toStdString();
	return params;
}

void MainWindow::onDemoChanged(int index) {
	const QString demoId = m_demoSelect->itemData(index).toString();
	if (demoId.isEmpty()) {
		m_activeDemoId.clear();
		m_demoDescription->clear();
		return;
	}

	m_activeDemoId = demoId;
	for (const auto& demo : m_demos) {
		if (demo.id == demoId.toStdString()) {
			m_demoDescription->setText(QString::fromStdString(demo.description));
			break;
		}
	}

	std::string error;
	if (const auto params = SimulationConfig::formParamsFromDemo(m_repoRoot, demoId.toStdString(), error)) {
		applyFormParams(*params);
		m_savedDemoParams = *params;
		setParamsCustomized(false);
	} else {
		m_statusLabel->setText(QString::fromStdString(error));
	}
}

void MainWindow::onRunClicked() {
	const FormParams params = collectFormParams();
	if (const auto validationError = SimulationConfig::validateFormParams(params)) {
		QMessageBox::warning(this, QStringLiteral("Invalid parameters"), QString::fromStdString(*validationError));
		return;
	}

	m_animationTimer->stop();
	m_chartPanel->setSimulationContext(params);
	setRunning(true);
	m_statusLabel->clear();

	QMetaObject::invokeMethod(
		m_worker,
		"runSimulation",
		Qt::QueuedConnection,
		Q_ARG(FormParams, params));
}

void MainWindow::onResetClicked() {
	applyFormParams(m_savedDemoParams);
	setParamsCustomized(false);
}

void MainWindow::onProgressUpdated(int percent, const QString& message) {
	m_progressBar->setValue(percent);
	m_statusLabel->setText(message);
}

void MainWindow::onSimulationFinished(const nlohmann::json& result) {
	setRunning(false);
	m_statusLabel->setText(QStringLiteral("Simulation complete — playing animation."));
	m_chartPanel->renderResults(result);
	updateSummary(result);

	if (result.contains("phaseFrames")) {
		m_frameCount = static_cast<int>(result["phaseFrames"].size());
	} else {
		m_frameCount = 0;
	}

	m_currentFrameIndex = 0;
	m_playbackPosition = 0.0;
	m_frameSlider->setEnabled(m_frameCount > 0);
	m_frameSlider->setMaximum(std::max(0, m_frameCount - 1));
	m_updatingSlider = true;
	m_frameSlider->setValue(0);
	m_updatingSlider = false;
	seekPlayback(0.0, false);

	if (auto* summaryGroup = qobject_cast<QGroupBox*>(m_statEseRatio->parent()->parent())) {
		summaryGroup->show();
	}

	if (m_frameCount > 1) {
		onPlayClicked();
	}
}

void MainWindow::onSimulationFailed(const QString& error) {
	setRunning(false);
	m_statusLabel->setText(error);
	QMessageBox::critical(this, QStringLiteral("Simulation failed"), error);
}

void MainWindow::setRunning(bool running) {
	m_runButton->setDisabled(running);
	m_demoSelect->setDisabled(running);
	m_progressBar->setVisible(running);
	if (running) {
		m_progressBar->setValue(0);
	}
}

void MainWindow::updateSummary(const nlohmann::json& result) {
	const auto ese = result.at("ese").get<std::vector<double>>();
	const auto ke = result.at("ke").get<std::vector<std::vector<double>>>();
	const size_t n = ese.size();
	if (n == 0) {
		return;
	}

	const double earlyEse = meanSlice(ese, static_cast<size_t>(std::floor(static_cast<double>(n) * 0.1)), static_cast<size_t>(std::floor(static_cast<double>(n) * 0.3)));
	const double lateEse = meanSlice(ese, static_cast<size_t>(std::floor(static_cast<double>(n) * 0.8)), n);
	const double eseRatio = earlyEse > 0.0 ? lateEse / earlyEse : 0.0;

	double total0 = ese.front();
	double totalN = ese.back();
	for (const auto& speciesKe : ke) {
		if (!speciesKe.empty()) {
			total0 += speciesKe.front();
			totalN += speciesKe.back();
		}
	}
	const double energyDrift = total0 > 0.0 ? std::abs(totalN - total0) / total0 : 0.0;
	const double peakEse = *std::max_element(ese.begin(), ese.end());
	const double peakRatio = ese.front() > 0.0 ? peakEse / ese.front() : 0.0;

	double ekRatio = 0.0;
	const FormParams params = collectFormParams();
	if (result.contains("phaseFrames")) {
		const auto frames = result["phaseFrames"].get<std::vector<DATA_STRUCTS::Frame>>();
		if (frames.size() >= 2) {
			auto amplitudeAt = [&](size_t index) {
				const auto& field = frames[index].electricField;
				const int gridN = static_cast<int>(field.size()) - 1;
				if (gridN <= 0) {
					return 0.0;
				}
				const double k = (2.0 * 3.14159265358979323846 * params.spatialPerturbationMode) / params.spatialLength;
				double re = 0.0;
				double im = 0.0;
				for (int j = 0; j < gridN; ++j) {
					const double x = (static_cast<double>(j) / gridN) * params.spatialLength;
					const double phase = k * x;
					re += field[static_cast<size_t>(j)] * std::cos(phase);
					im += field[static_cast<size_t>(j)] * std::sin(phase);
				}
				re /= gridN;
				im /= gridN;
				return std::sqrt(re * re + im * im);
			};
			const double start = amplitudeAt(0);
			const double end = amplitudeAt(frames.size() - 1);
			if (start > 0.0) {
				ekRatio = end / start;
			}
		}
	}

	m_statEseRatio->setText(eseRatio > 0.0 ? QString::number(eseRatio, 'f', 2) : QStringLiteral("—"));
	m_statEnergyDrift->setText(QString("%1%").arg(energyDrift * 100.0, 0, 'f', 1));
	m_statEsePeak->setText(peakRatio > 0.0 ? QString("%1×").arg(peakRatio, 0, 'f', 2) : QStringLiteral("—"));
	m_statEkRatio->setText(params.dimension == 3
		? QStringLiteral("N/A")
		: (ekRatio > 0.0 ? QString::number(ekRatio, 'f', 2) : QStringLiteral("—")));
}

void MainWindow::markParamsCustomized() {
	if (m_savedDemoParams.numParticles == 0 && m_activeDemoId.isEmpty()) {
		return;
	}
	const FormParams current = collectFormParams();
	const FormParams saved = m_savedDemoParams;
	const bool same =
		current.dimension == saved.dimension &&
		current.spatialLength == saved.spatialLength &&
		current.spatialLengthY == saved.spatialLengthY &&
		current.spatialLengthZ == saved.spatialLengthZ &&
		current.numParticles == saved.numParticles &&
		current.timeSteps == saved.timeSteps &&
		current.timeStepSize == saved.timeStepSize &&
		current.numGrid == saved.numGrid &&
		current.numGridY == saved.numGridY &&
		current.numGridZ == saved.numGridZ &&
		current.spatialPerturbationMode == saved.spatialPerturbationMode &&
		current.driftVelocity == saved.driftVelocity &&
		current.numSpecies == saved.numSpecies &&
		current.spatialPerturbationAmplitude == saved.spatialPerturbationAmplitude &&
		current.thermalVelocity == saved.thermalVelocity &&
		current.plasmaFrequency == saved.plasmaFrequency &&
		current.chargeMassRatio == saved.chargeMassRatio &&
		current.framePeriod == saved.framePeriod &&
		current.spatialPerturbationWaveform == saved.spatialPerturbationWaveform;
	setParamsCustomized(!same);
}

void MainWindow::setParamsCustomized(bool customized) {
	m_paramsCustomized = customized;
	m_customizedBadge->setVisible(customized);
	m_resetButton->setEnabled(customized);
}

void MainWindow::updateFrameLabel(double continuousIndex) {
	const double simTime = m_chartPanel->continuousTime(continuousIndex);
	const int keyframe = static_cast<int>(std::floor(continuousIndex + 1e-9)) + 1;
	m_frameLabel->setText(QString("Frame %1 / %2 · t ≈ %3")
		.arg(std::min(keyframe, m_frameCount))
		.arg(m_frameCount)
		.arg(simTime, 0, 'f', 2));
}

void MainWindow::seekPlayback(double continuousIndex, bool updateSlider) {
	if (m_frameCount <= 0) {
		return;
	}
	const double maxIndex = static_cast<double>(m_frameCount - 1);
	m_playbackPosition = std::clamp(continuousIndex, 0.0, maxIndex);
	m_currentFrameIndex = static_cast<int>(std::floor(m_playbackPosition + 1e-9));
	m_chartPanel->showPlaybackPosition(m_playbackPosition);
	updateFrameLabel(m_playbackPosition);

	if (updateSlider) {
		m_updatingSlider = true;
		m_frameSlider->setValue(m_currentFrameIndex);
		m_updatingSlider = false;
	}
}

void MainWindow::onPlayClicked() {
	if (m_frameCount <= 0) {
		return;
	}
	if (m_playbackPosition >= static_cast<double>(m_frameCount - 1) - 1e-9) {
		seekPlayback(0.0, true);
	}
	m_animationTimer->start(m_animSpeedControl->value());
}

void MainWindow::onPauseClicked() {
	m_animationTimer->stop();
}

void MainWindow::onAnimSpeedChanged(int value) {
	if (m_animationTimer->isActive()) {
		m_animationTimer->setInterval(value);
	}
}

void MainWindow::onFirstFrame() {
	if (m_frameCount <= 0) {
		return;
	}
	m_animationTimer->stop();
	seekPlayback(0.0, true);
}

void MainWindow::onLastFrame() {
	if (m_frameCount <= 0) {
		return;
	}
	m_animationTimer->stop();
	seekPlayback(static_cast<double>(m_frameCount - 1), true);
}

void MainWindow::onPrevFrame() {
	if (m_frameCount <= 0) {
		return;
	}
	m_animationTimer->stop();
	double next = m_playbackPosition - 1.0;
	if (next < 0.0) {
		next = m_loopPlayback->isChecked() ? static_cast<double>(m_frameCount - 1) : 0.0;
	}
	seekPlayback(next, true);
}

void MainWindow::onNextFrame() {
	if (m_frameCount <= 0) {
		return;
	}
	double next = m_playbackPosition + 1.0;
	if (next > static_cast<double>(m_frameCount - 1)) {
		next = m_loopPlayback->isChecked() ? 0.0 : static_cast<double>(m_frameCount - 1);
	}
	seekPlayback(next, true);
}

void MainWindow::onFrameSliderPressed() {
	m_wasPlayingBeforeScrub = m_animationTimer->isActive();
	m_animationTimer->stop();
}

void MainWindow::onFrameSliderReleased() {
	if (m_wasPlayingBeforeScrub) {
		m_animationTimer->start(m_animSpeedControl->value());
	}
}

void MainWindow::onFrameSliderChanged(int value) {
	if (m_updatingSlider) {
		return;
	}
	seekPlayback(static_cast<double>(value), false);
}

void MainWindow::onAnimationTick() {
	if (m_frameCount <= 0) {
		return;
	}
	const double maxIndex = static_cast<double>(m_frameCount - 1);
	const double step = 1.0 / static_cast<double>(kInterpSteps);
	double next = m_playbackPosition + step;
	if (next > maxIndex + 1e-9) {
		if (m_loopPlayback->isChecked()) {
			next = 0.0;
		} else {
			seekPlayback(maxIndex, true);
			m_animationTimer->stop();
			m_statusLabel->setText(QStringLiteral("Animation finished."));
			return;
		}
	}
	seekPlayback(next, true);
}
