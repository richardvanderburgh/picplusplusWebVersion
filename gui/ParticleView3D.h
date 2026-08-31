#pragma once

#include <QWidget>

#include <DataStructs.h>

#include <vector>

class Q3DScatter;
class QCheckBox;
class QScatter3DSeries;
class QTimer;

class ParticleView3D : public QWidget {
	Q_OBJECT

public:
	explicit ParticleView3D(QWidget* parent = nullptr);

	void configureDomain(double lengthX, double lengthY, double lengthZ, int numSpecies, int particlesPerSpecies);
	void clearView();
	void setFrame(const DATA_STRUCTS::Frame& frame, const DATA_STRUCTS::Frame* previousFrame);

public slots:
	void resetCamera();

private slots:
	void onTrailsToggled(bool enabled);
	void onAutoRotateToggled(bool enabled);
	void onAutoRotateTick();

private:
	void clearSeries();
	void fillSeries(QScatter3DSeries* series, int species, const DATA_STRUCTS::Frame& frame);
	float itemSizeForDomain() const;

	Q3DScatter* m_scatter = nullptr;
	QWidget* m_container = nullptr;
	QCheckBox* m_showTrails = nullptr;
	QCheckBox* m_autoRotate = nullptr;
	QTimer* m_rotateTimer = nullptr;

	std::vector<QScatter3DSeries*> m_series;
	std::vector<QScatter3DSeries*> m_trailSeries;

	double m_lengthX = 1.0;
	double m_lengthY = 1.0;
	double m_lengthZ = 1.0;
	int m_numSpecies = 0;
	int m_particlesPerSpecies = 0;
};
