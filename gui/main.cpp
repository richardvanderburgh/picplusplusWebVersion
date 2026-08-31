#include "MainWindow.h"

#include <QApplication>

#include <filesystem>

#include "SimulationConfig.h"

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	qRegisterMetaType<FormParams>("FormParams");
	QApplication::setApplicationName(QStringLiteral("PIC++ GUI"));
	QApplication::setOrganizationName(QStringLiteral("PIC++"));

	std::string repoRoot = std::filesystem::current_path().string();
	if (argc >= 2) {
		repoRoot = argv[1];
	}

	MainWindow window(repoRoot);
	window.show();
	return app.exec();
}
