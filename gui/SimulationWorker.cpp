#include "SimulationWorker.h"

#include <PICPlusPlus.h>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <streambuf>
#include <string>

#include "JsonConfigLoader.hpp"

namespace {

class ProgressStreambuf final : public std::streambuf {
public:
	explicit ProgressStreambuf(std::function<void(int)> onProgress)
		: m_onProgress(std::move(onProgress)) {}

protected:
	int overflow(int character) override {
		if (character == '\n') {
			flushLine();
		} else if (character != EOF) {
			m_buffer.push_back(static_cast<char>(character));
		}
		return character;
	}

	int sync() override {
		flushLine();
		return 0;
	}

private:
	void flushLine() {
		if (m_buffer.rfind("PICPP_PROGRESS ", 0) == 0 && m_onProgress) {
			try {
				m_onProgress(std::stoi(m_buffer.substr(15)));
			} catch (...) {
			}
		}
		m_buffer.clear();
	}

	std::string m_buffer;
	std::function<void(int)> m_onProgress;
};

} // namespace

SimulationWorker::SimulationWorker(QObject* parent)
	: QObject(parent) {}

void SimulationWorker::runSimulation(const FormParams& params) {
	if (const auto validationError = SimulationConfig::validateFormParams(params)) {
		emit simulationFailed(QString::fromStdString(*validationError));
		return;
	}

	const nlohmann::json config = SimulationConfig::buildConfig(params);
	DATA_STRUCTS::InputVariables inputVariables = loadJSONFile(config);

	if (const auto validationError = validateSimulationParams(inputVariables.simulationParams)) {
		emit simulationFailed(QString::fromStdString(*validationError));
		return;
	}

	qputenv("PICPP_PROGRESS", "1");

	ProgressStreambuf progressBuf([this](int percent) {
		emit progressUpdated(percent, QString("Time loop %1%").arg(percent));
	});
	const auto previousBuf = std::cerr.rdbuf(&progressBuf);

	emit progressUpdated(0, QStringLiteral("Starting simulation…"));

	try {
		PIC_PLUS_PLUS::PICPlusPlus picPlusPlus(inputVariables);
		const auto jsonResult = picPlusPlus.initialize();

		std::cerr.rdbuf(previousBuf);

		if (!jsonResult.has_value()) {
			emit simulationFailed(QStringLiteral("Simulation failed to produce results."));
			return;
		}

		emit progressUpdated(100, QStringLiteral("Done"));
		emit simulationFinished(*jsonResult);
	} catch (const std::exception& ex) {
		std::cerr.rdbuf(previousBuf);
		emit simulationFailed(QString::fromStdString(ex.what()));
	} catch (...) {
		std::cerr.rdbuf(previousBuf);
		emit simulationFailed(QStringLiteral("Unknown simulation error."));
	}
}
