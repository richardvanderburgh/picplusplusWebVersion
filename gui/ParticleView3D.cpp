#include "ParticleView3D.h"

#include <QCheckBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QTimer>
#include <QVBoxLayout>

#include <cmath>
#include <algorithm>

#include <QtDataVisualization/Q3DCamera>
#include <QtDataVisualization/Q3DScatter>
#include <QtDataVisualization/Q3DTheme>
#include <QtDataVisualization/QScatter3DSeries>
#include <QtDataVisualization/QScatterDataItem>

namespace {

const QColor SPECIES_COLORS[] = {
	QColor("#2563eb"),
	QColor("#dc2626"),
	QColor("#16a34a"),
	QColor("#ea580c"),
};

QColor trailColor(const QColor& base) {
	return QColor(base.red(), base.green(), base.blue(), 140);
}

} // namespace

ParticleView3D::ParticleView3D(QWidget* parent)
	: QWidget(parent) {
	m_scatter = new Q3DScatter();
	m_scatter->setShadowQuality(QAbstract3DGraph::ShadowQualityNone);
	m_scatter->setSelectionMode(QAbstract3DGraph::SelectionNone);
	m_scatter->activeTheme()->setType(Q3DTheme::ThemeQt);
	m_scatter->activeTheme()->setBackgroundEnabled(true);
	m_scatter->activeTheme()->setGridEnabled(true);
	m_scatter->activeTheme()->setLabelBorderEnabled(true);
	m_scatter->activeTheme()->setColorStyle(Q3DTheme::ColorStyleUniform);

	auto* resetButton = new QPushButton(QStringLiteral("Reset view"));
	m_showTrails = new QCheckBox(QStringLiteral("Show trails"));
	m_showTrails->setChecked(true);
	m_autoRotate = new QCheckBox(QStringLiteral("Auto-rotate"));
	m_autoRotate->setChecked(true);
	m_rotateTimer = new QTimer(this);
	m_rotateTimer->setInterval(33);
	connect(resetButton, &QPushButton::clicked, this, &ParticleView3D::resetCamera);
	connect(m_showTrails, &QCheckBox::toggled, this, &ParticleView3D::onTrailsToggled);
	connect(m_autoRotate, &QCheckBox::toggled, this, &ParticleView3D::onAutoRotateToggled);
	connect(m_rotateTimer, &QTimer::timeout, this, &ParticleView3D::onAutoRotateTick);

	auto* hint = new QLabel(QStringLiteral("Drag to rotate · scroll to zoom · right-drag to pan"));
	hint->setStyleSheet(QStringLiteral("color: #64748b;"));

	auto* toolbar = new QHBoxLayout();
	toolbar->addWidget(resetButton);
	toolbar->addWidget(m_showTrails);
	toolbar->addWidget(m_autoRotate);
	toolbar->addStretch(1);
	toolbar->addWidget(hint);

	m_container = QWidget::createWindowContainer(m_scatter, this);
	m_container->setMinimumSize(420, 320);
	m_container->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	m_container->setFocusPolicy(Qt::StrongFocus);

	auto* layout = new QVBoxLayout(this);
	layout->setContentsMargins(0, 0, 0, 0);
	layout->addLayout(toolbar);
	layout->addWidget(m_container, 1);
}

void ParticleView3D::clearSeries() {
	for (auto* series : m_series) {
		m_scatter->removeSeries(series);
		delete series;
	}
	m_series.clear();

	for (auto* series : m_trailSeries) {
		m_scatter->removeSeries(series);
		delete series;
	}
	m_trailSeries.clear();
}

void ParticleView3D::clearView() {
	clearSeries();
	m_numSpecies = 0;
	m_particlesPerSpecies = 0;
}

float ParticleView3D::itemSizeForDomain() const {
	const double volume = std::max(m_lengthX * m_lengthY * m_lengthZ, 1e-6);
	const double referenceVolume = 6.28318530717958 * 6.28318530717958 * 6.28318530717958;
	const float scale = static_cast<float>(std::cbrt(referenceVolume / volume));
	float size = 0.12f * scale;
	if (m_particlesPerSpecies > 2000) {
		size *= 0.65f;
	} else if (m_particlesPerSpecies > 800) {
		size *= 0.8f;
	}
	return std::clamp(size, 0.03f, 0.2f);
}

void ParticleView3D::configureDomain(double lengthX, double lengthY, double lengthZ, int numSpecies, int particlesPerSpecies) {
	clearSeries();

	m_lengthX = lengthX;
	m_lengthY = lengthY;
	m_lengthZ = lengthZ;
	m_numSpecies = numSpecies;
	m_particlesPerSpecies = particlesPerSpecies;

	m_scatter->axisX()->setTitle(QStringLiteral("x"));
	m_scatter->axisY()->setTitle(QStringLiteral("y"));
	m_scatter->axisZ()->setTitle(QStringLiteral("z"));
	m_scatter->axisX()->setRange(0.0, lengthX);
	m_scatter->axisY()->setRange(0.0, lengthY);
	m_scatter->axisZ()->setRange(0.0, lengthZ);
	m_scatter->axisX()->setSegmentCount(4);
	m_scatter->axisY()->setSegmentCount(4);
	m_scatter->axisZ()->setSegmentCount(4);
	m_scatter->axisX()->setLabelFormat(QStringLiteral("%.1f"));
	m_scatter->axisY()->setLabelFormat(QStringLiteral("%.1f"));
	m_scatter->axisZ()->setLabelFormat(QStringLiteral("%.1f"));

	const float itemSize = itemSizeForDomain();
	for (int species = 0; species < numSpecies; ++species) {
		auto* series = new QScatter3DSeries();
		series->setItemSize(itemSize);
		series->setMesh(QAbstract3DSeries::MeshSphere);
		series->setBaseColor(SPECIES_COLORS[species % 4]);
		series->setName(QStringLiteral("Species %1").arg(species));
		m_scatter->addSeries(series);
		m_series.push_back(series);

		auto* trail = new QScatter3DSeries();
		trail->setItemSize(itemSize * 0.55f);
		trail->setMesh(QAbstract3DSeries::MeshPoint);
		trail->setBaseColor(trailColor(SPECIES_COLORS[species % 4]));
		trail->setName(QStringLiteral("Species %1 (prev)").arg(species));
		trail->setVisible(m_showTrails->isChecked());
		m_scatter->addSeries(trail);
		m_trailSeries.push_back(trail);
	}

	resetCamera();
	if (m_autoRotate->isChecked()) {
		m_rotateTimer->start();
	}
}

void ParticleView3D::fillSeries(QScatter3DSeries* series, int species, const DATA_STRUCTS::Frame& frame) {
	auto* items = new QScatterDataArray();
	items->reserve(static_cast<int>(frame.particles.size()));
	for (const auto& particle : frame.particles) {
		if (particle.species == species) {
			items->append(QScatterDataItem(QVector3D(
				static_cast<float>(particle.position),
				static_cast<float>(particle.positionY),
				static_cast<float>(particle.positionZ))));
		}
	}
	series->dataProxy()->resetArray(items);
}

void ParticleView3D::setFrame(const DATA_STRUCTS::Frame& frame, const DATA_STRUCTS::Frame* previousFrame) {
	if (m_series.empty()) {
		return;
	}

	for (int species = 0; species < m_numSpecies; ++species) {
		fillSeries(m_series[static_cast<size_t>(species)], species, frame);
		if (m_showTrails->isChecked() && previousFrame != nullptr) {
			fillSeries(m_trailSeries[static_cast<size_t>(species)], species, *previousFrame);
		} else if (!m_trailSeries.empty()) {
			m_trailSeries[static_cast<size_t>(species)]->dataProxy()->resetArray(new QScatterDataArray());
		}
	}
}

void ParticleView3D::resetCamera() {
	if (m_scatter == nullptr) {
		return;
	}
	auto* camera = m_scatter->scene()->activeCamera();
	camera->setCameraPreset(Q3DCamera::CameraPresetIsometricLeft);
	camera->setZoomLevel(100.0f);
}

void ParticleView3D::onTrailsToggled(bool enabled) {
	for (auto* series : m_trailSeries) {
		series->setVisible(enabled);
	}
}

void ParticleView3D::onAutoRotateToggled(bool enabled) {
	if (enabled) {
		m_rotateTimer->start();
	} else {
		m_rotateTimer->stop();
	}
}

void ParticleView3D::onAutoRotateTick() {
	if (m_scatter == nullptr) {
		return;
	}
	auto* camera = m_scatter->scene()->activeCamera();
	float xRot = camera->xRotation() + 0.6f;
	if (xRot > 180.0f) {
		xRot -= 360.0f;
	}
	camera->setXRotation(xRot);
}
