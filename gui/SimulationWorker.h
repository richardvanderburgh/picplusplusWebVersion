#pragma once

#include <QObject>
#include <QString>

#include <nlohmann/json.hpp>

#include "SimulationConfig.h"

class SimulationWorker : public QObject {
	Q_OBJECT

public:
	explicit SimulationWorker(QObject* parent = nullptr);

public slots:
	void runSimulation(const FormParams& params);

signals:
	void progressUpdated(int percent, const QString& message);
	void simulationFinished(const nlohmann::json& result);
	void simulationFailed(const QString& error);
};
