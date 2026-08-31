#pragma once

#include <cmath>

// Boris rotation in physical velocity units.
// tVec = (q/m) * B * (dt/2); returns updated (vx, vy, vz).
inline void borisRotate(
	double& vx,
	double& vy,
	double& vz,
	double tx,
	double ty,
	double tz) {

	const double vPrimeX = vx + (vy * tz - vz * ty);
	const double vPrimeY = vy + (vz * tx - vx * tz);
	const double vPrimeZ = vz + (vx * ty - vy * tx);

	const double tSquared = tx * tx + ty * ty + tz * tz;
	const double sx = 2.0 * tx / (1.0 + tSquared);
	const double sy = 2.0 * ty / (1.0 + tSquared);
	const double sz = 2.0 * tz / (1.0 + tSquared);

	vx += vPrimeY * sz - vPrimeZ * sy;
	vy += vPrimeZ * sx - vPrimeX * sz;
	vz += vPrimeX * sy - vPrimeY * sx;
}

// Full Boris push with electric half-kicks.
// ae{X,Y,Z} convert sampled E → code-unit Δv (same as electrostatic accel).
// Velocities enter/leave in code units (cells/timestep): v_code = v_phys * dt/d{ℓ}.
inline void borisPushCodeUnits(
	double& vxCode,
	double& vyCode,
	double& vzCode,
	double exSample,
	double eySample,
	double ezSample,
	double aeX,
	double aeY,
	double aeZ,
	double qm,
	double dt,
	double dxdt,
	double dydt,
	double dzdt,
	double bx,
	double by,
	double bz) {

	// Electric half-kick in code units.
	vxCode += 0.5 * aeX * exSample;
	vyCode += 0.5 * aeY * eySample;
	vzCode += 0.5 * aeZ * ezSample;

	// Magnetic rotation in physical velocity units.
	double vx = vxCode * dxdt;
	double vy = vyCode * dydt;
	double vz = vzCode * dzdt;

	const double halfDt = 0.5 * dt;
	const double tx = qm * bx * halfDt;
	const double ty = qm * by * halfDt;
	const double tz = qm * bz * halfDt;
	borisRotate(vx, vy, vz, tx, ty, tz);

	vxCode = vx / dxdt;
	vyCode = vy / dydt;
	vzCode = vz / dzdt;

	// Second electric half-kick.
	vxCode += 0.5 * aeX * exSample;
	vyCode += 0.5 * aeY * eySample;
	vzCode += 0.5 * aeZ * ezSample;
}
