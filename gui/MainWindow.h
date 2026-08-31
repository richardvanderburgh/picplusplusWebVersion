#pragma once

#include <QMainWindow>
#include <QTimer>

#include <nlohmann/json.hpp>

#include "ChartPanel.h"
#include "SimulationConfig.h"

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QFormLayout;
class QGroupBox;
class QLabel;
class QProgressBar;
class QPushButton;
class QSlider;
class QSpinBox;
class SimulationWorker;

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	explicit MainWindow(const std::string& repoRoot, QWidget* parent = nullptr);

private slots:
	void onDemoChanged(int index);
	void onDimensionChanged(int index);
	void onRunClicked();
	void onResetClicked();
	void onProgressUpdated(int percent, const QString& message);
	void onSimulationFinished(const nlohmann::json& result);
	void onSimulationFailed(const QString& error);
	void onPlayClicked();
	void onPauseClicked();
	void onPrevFrame();
	void onNextFrame();
	void onFrameSliderChanged(int value);
	void onFirstFrame();
	void onLastFrame();
	void onFrameSliderPressed();
	void onFrameSliderReleased();
	void onAnimationTick();
	void onAnimSpeedChanged(int value);

	void onCopyConfig();
	void onExportEnergyCsv();
	void onExportModeCsv();
	void onExportFramesCsv();
	void onExportChartPng();
	void onLessonChanged(int index);
	void onLessonPrev();
	void onLessonNext();
	void onSweepClicked();

private:
	void populateDemoSelector();
	void populateLessonSelector();
	void setupMenus();
	void applyFormParams(const FormParams& params);
	FormParams collectFormParams() const;
	nlohmann::json currentRunConfig() const;
	void startSimulation(const nlohmann::json& config, const FormParams& contextParams);
	void setRunning(bool running);
	void updateSummary(const nlohmann::json& result);
	void markParamsCustomized();
	void setParamsCustomized(bool customized);
	void updateDimensionControls();
	void updateFrameLabel(double continuousIndex);
	void seekPlayback(double continuousIndex, bool updateSlider);
	void applyLessonStep();
	void continueSweepIfNeeded(const nlohmann::json& result);
	QString exportPath(const QString& suffix) const;

	static constexpr int kInterpSteps = 4;

	std::string m_repoRoot;
	std::vector<DemoEntry> m_demos;
	std::vector<LessonEntry> m_lessons;
	FormParams m_savedDemoParams;
	QString m_activeDemoId;
	bool m_paramsCustomized = false;
	nlohmann::json m_lastConfig;
	nlohmann::json m_lastResult;

	QComboBox* m_demoSelect = nullptr;
	QLabel* m_demoDescription = nullptr;
	QPushButton* m_runButton = nullptr;
	QPushButton* m_resetButton = nullptr;
	QProgressBar* m_progressBar = nullptr;
	QLabel* m_statusLabel = nullptr;
	QLabel* m_customizedBadge = nullptr;

	QComboBox* m_lessonSelect = nullptr;
	QLabel* m_lessonPrompt = nullptr;
	QPushButton* m_lessonPrevButton = nullptr;
	QPushButton* m_lessonNextButton = nullptr;
	int m_lessonIndex = -1;
	int m_lessonStep = 0;

	QComboBox* m_sweepParam = nullptr;
	QDoubleSpinBox* m_sweepStart = nullptr;
	QDoubleSpinBox* m_sweepEnd = nullptr;
	QSpinBox* m_sweepCount = nullptr;
	QPushButton* m_sweepButton = nullptr;
	bool m_sweepActive = false;
	int m_sweepIndex = 0;
	std::vector<FormParams> m_sweepJobs;
	std::vector<SweepSeriesData> m_sweepResults;

	QDoubleSpinBox* m_spatialLength = nullptr;
	QSpinBox* m_numParticles = nullptr;
	QSpinBox* m_timeSteps = nullptr;
	QDoubleSpinBox* m_timeStepSize = nullptr;
	QSpinBox* m_numGrid = nullptr;
	QSpinBox* m_spatialPerturbationMode = nullptr;
	QDoubleSpinBox* m_driftVelocity = nullptr;
	QSpinBox* m_numSpecies = nullptr;
	QDoubleSpinBox* m_spatialPerturbationAmplitude = nullptr;
	QDoubleSpinBox* m_thermalVelocity = nullptr;
	QDoubleSpinBox* m_plasmaFrequency = nullptr;
	QDoubleSpinBox* m_chargeMassRatio = nullptr;
	QComboBox* m_spatialPerturbationWaveform = nullptr;
	QComboBox* m_dimension = nullptr;
	QSpinBox* m_framePeriod = nullptr;
	QSpinBox* m_numGridY = nullptr;
	QSpinBox* m_numGridZ = nullptr;
	QDoubleSpinBox* m_spatialLengthY = nullptr;
	QDoubleSpinBox* m_spatialLengthZ = nullptr;
	QFormLayout* m_paramsLayout = nullptr;
	QWidget* m_grid3DRowY = nullptr;
	QWidget* m_grid3DRowZ = nullptr;
	QWidget* m_length3DRowY = nullptr;
	QWidget* m_length3DRowZ = nullptr;

	QLabel* m_statEseRatio = nullptr;
	QLabel* m_statEnergyDrift = nullptr;
	QLabel* m_statEkRatio = nullptr;
	QLabel* m_statEsePeak = nullptr;
	QLabel* m_statGrowth = nullptr;
	QLabel* m_statTheory = nullptr;
	QGroupBox* m_summaryGroup = nullptr;

	QSlider* m_frameSlider = nullptr;
	QLabel* m_frameLabel = nullptr;
	QPushButton* m_playButton = nullptr;
	QPushButton* m_pauseButton = nullptr;
	QCheckBox* m_loopPlayback = nullptr;
	QSpinBox* m_animSpeedControl = nullptr;

	ChartPanel* m_chartPanel = nullptr;

	QThread* m_workerThread = nullptr;
	SimulationWorker* m_worker = nullptr;

	QTimer* m_animationTimer = nullptr;
	bool m_wasPlayingBeforeScrub = false;
	bool m_updatingSlider = false;
	double m_playbackPosition = 0.0;
	int m_currentFrameIndex = 0;
	int m_frameCount = 0;
};
