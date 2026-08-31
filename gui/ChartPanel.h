#pragma once

#include <QChartView>
#include <QStackedWidget>
#include <QTabWidget>
#include <QWidget>

#include <nlohmann/json.hpp>

#include <DataStructs.h>
#include "SimulationConfig.h"

QT_BEGIN_NAMESPACE
class QChart;
class QLineSeries;
class QScatterSeries;
class QValueAxis;
class QLogValueAxis;
QT_END_NAMESPACE

class ChartPanel : public QWidget {
	Q_OBJECT

public:
	explicit ChartPanel(QWidget* parent = nullptr);

	void clearResults();
	void setSimulationContext(const FormParams& params);
	void renderResults(const nlohmann::json& result);
	void showFrame(int frameIndex);

private:
	void setupCharts();
	void renderEnergyPlot(const nlohmann::json& result);
	void renderModePlot(const nlohmann::json& result);
	void renderPhasePlot(const DATA_STRUCTS::Frame& frame);
	void renderFieldPlot(const DATA_STRUCTS::Frame& frame);
	void renderProjectionPlots(const DATA_STRUCTS::Frame& frame);
	void clearScatterSeries(std::vector<QScatterSeries*>& seriesList, QChart* chart);
	double fourierModeAmplitude(const std::vector<double>& field, int mode, double domainLength) const;

	FormParams m_params;
	bool m_is3D = false;
	nlohmann::json m_result;
	std::vector<DATA_STRUCTS::Frame> m_frames;
	std::vector<double> m_modeAmplitudes;
	std::vector<double> m_modeTimes;
	double m_phaseVelocityMin = -1.0;
	double m_phaseVelocityMax = 1.0;

	QTabWidget* m_tabs = nullptr;
	QStackedWidget* m_phaseStack = nullptr;

	QChartView* m_phaseView = nullptr;
	QChartView* m_projXYView = nullptr;
	QChartView* m_projXZView = nullptr;
	QChartView* m_projYZView = nullptr;
	QChartView* m_fieldView = nullptr;
	QChartView* m_energyView = nullptr;
	QChartView* m_modeView = nullptr;

	QChart* m_phaseChart = nullptr;
	QChart* m_projXYChart = nullptr;
	QChart* m_projXZChart = nullptr;
	QChart* m_projYZChart = nullptr;
	QChart* m_fieldChart = nullptr;
	QChart* m_energyChart = nullptr;
	QChart* m_modeChart = nullptr;

	std::vector<QScatterSeries*> m_phaseSeries;
	std::vector<QScatterSeries*> m_projXYSeries;
	std::vector<QScatterSeries*> m_projXZSeries;
	std::vector<QScatterSeries*> m_projYZSeries;
	QLineSeries* m_fieldSeries = nullptr;
	std::vector<QLineSeries*> m_keSeries;
	QLineSeries* m_eseSeries = nullptr;
	QLineSeries* m_totalSeries = nullptr;
	QLineSeries* m_modeSeries = nullptr;

	QValueAxis* m_phaseXAxis = nullptr;
	QValueAxis* m_phaseYAxis = nullptr;
};
