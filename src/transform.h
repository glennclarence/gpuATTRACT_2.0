/*
 * CoordTrafo.h
 *
 *  Created on: Aug 11, 2016
 *      Author: uwe
 */

#ifndef SRC_TRANSFORM_H_
#define SRC_TRANSFORM_H_

#include "Vec3.h"
#include "RotMat.h"
#include "MatrixFunctions.h"

namespace as {

template<typename REAL>
void rotate_translate(
		REAL const* x,
		REAL const* y,
		REAL const* z,
		Vec3<REAL> const& displacement,
		Vec3<REAL> const& ang,
		unsigned const& numAtoms,
		REAL* xTr,
		REAL* yTr,
		REAL* zTr)
{
	const RotMat<REAL> rotMat = euler2rotmat(ang.x, ang.y, ang.z);
	for (unsigned i = 0; i < numAtoms; ++i) {
		Vec3<REAL> posAtom(x[i], y[i], z[i]);

		posAtom = rotMat*posAtom;
		posAtom += displacement;

		xTr[i] = posAtom.x;
		yTr[i] = posAtom.y;
		zTr[i] = posAtom.z;
	}
}

} // namespace as




#endif /* SRC_TRANSFORM_H_ */