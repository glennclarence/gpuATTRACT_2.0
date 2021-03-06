/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef META_H_
#define META_H_

#include <Eigen/Core>
#include "Types_6D.h"
#include "Types_6D_Modes.h"

#define OBJGRAD(dof, energy)	\
	do { 						\
		/*std::cout << "\t" << "request=" << Vector2extDOF(dof) << std::endl;*/ \
		state = dof;			\
		ca();					\
		energy = objective;		\
		/*std::cout << "\t" << "result=" << ObjGrad2extEnGrad(energy) << std::endl;*/						\
	} while(0)


#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x,y) ((x < y) ? y : x)
#endif

namespace as {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Scalar = Eigen::VectorXd::Scalar;

struct ObjGrad {
	double obj; // function value
	Vector grad; // gradients
};


template<typename FIRST, typename SECOND>
class TypesConverter {
	static FIRST toFirst(SECOND const&);
	static SECOND toSecond(FIRST const&);
};


template<typename REAL>
class TypesConverter<DOF_6D<REAL>, Vector> {
public:
	static DOF_6D<REAL> toFirst(Vector const& vec) {
		DOF_6D<REAL> dof;
		dof.ang.x = vec(0);
		dof.ang.y = vec(1);
		dof.ang.z = vec(2);
		dof.pos.x = vec(3);
		dof.pos.y = vec(4);
		dof.pos.z = vec(5);
		return dof;
	}

	static Vector toSecond(DOF_6D<REAL> const& dof) {
		Vector vec(6);
		vec  << dof.ang.x, dof.ang.y, dof.ang.z,
				dof.pos.x, dof.pos.y , dof.pos.z;
		return vec;
	}
};

template<typename REAL>
class TypesConverter<Vector, DOF_6D<REAL>> {
public:
	static Vector toFirst(DOF_6D<REAL> const& dof) {
		return TypesConverter<DOF_6D<REAL>, Vector>::toSecond(dof);
	}

	static DOF_6D<REAL> toSecond(Vector const& vec) {
		return TypesConverter<DOF_6D<REAL>, Vector>::toFirst(vec);
	}
};

template<typename REAL>
class TypesConverter<Result_6D<REAL>, ObjGrad> {
public:
	static Result_6D<REAL> toFirst(ObjGrad const& objGrad) {
		Result_6D<REAL> enGrad;
		enGrad.E = objGrad.obj;
		enGrad.ang.x = objGrad.grad(0);
		enGrad.ang.y = objGrad.grad(1);
		enGrad.ang.z = objGrad.grad(2);
		enGrad.pos.x = objGrad.grad(3);
		enGrad.pos.y = objGrad.grad(4);
		enGrad.pos.z = objGrad.grad(5);
		return enGrad;
	}

	static ObjGrad toSecond (Result_6D<REAL> const& enGrad) {
		ObjGrad objGrad;
		objGrad.obj = enGrad.E;
		objGrad.grad = Vector(6);
		objGrad.grad  << -enGrad.ang.x, -enGrad.ang.y,  -enGrad.ang.z,
						 -enGrad.pos.x, -enGrad.pos.y , -enGrad.pos.z;
		return objGrad;
	}
};

template<typename REAL>
class TypesConverter<ObjGrad, Result_6D<REAL>> {
public:

	static ObjGrad toFirst (Result_6D<REAL> const& enGrad) {
		return TypesConverter<Result_6D<REAL>, ObjGrad>::toSecond(enGrad);
	}

	static Result_6D<REAL> toSecond(ObjGrad const& objGrad) {
		return TypesConverter<Result_6D<REAL>, ObjGrad>::toFirst(objGrad);
	}
};


// Typesconverter for modes

template<typename REAL>
class TypesConverter<DOF_6D_Modes<REAL>, Vector> {
public:
	static DOF_6D_Modes<REAL> toFirst(Vector const& vec) {
		DOF_6D_Modes<REAL> dof;
		dof._6D.ang.x = vec(0);
		dof._6D.ang.y = vec(1);
		dof._6D.ang.z = vec(2);
		dof._6D.pos.x = vec(3);
		dof._6D.pos.y = vec(4);
		dof._6D.pos.z = vec(5);
		for(int mode=0;mode< Common_Modes::numModesRec; mode++){
			dof.modesRec[mode]=vec(6 + mode);
		}
		for(int mode=0;mode< Common_Modes::numModesLig; mode++){
			dof.modesLig[mode]=vec(6 + mode + Common_Modes::numModesRec);
		}
		return dof;
	}

	static Vector toSecond(DOF_6D_Modes<REAL> const& dof) {
		Vector vec(DOF_MAX_NUMBER);
		vec(0) = dof._6D.ang.x;
		vec(1) = dof._6D.ang.y;
		vec(2) = dof._6D.ang.z;
		vec(3) = dof._6D.pos.x;
		vec(4) = dof._6D.pos.y;
		vec(5) = dof._6D.pos.z;
		for(int mode=0;mode< Common_Modes::numModesRec; mode++){
			vec(6 + mode)  = dof.modesRec[mode];
		}
		for(int mode=0;mode< Common_Modes::numModesLig; mode++){
			vec(6 + Common_Modes::numModesRec + mode)  = dof.modesLig[mode];
		}
		return vec;
	}
};

template<typename REAL>
class TypesConverter<Vector, DOF_6D_Modes<REAL>> {
public:
	static Vector toFirst(DOF_6D_Modes<REAL> const& dof) {
		return TypesConverter<DOF_6D_Modes<REAL>, Vector>::toSecond(dof);
	}

	static DOF_6D_Modes<REAL> toSecond(Vector const& vec) {
		return TypesConverter<DOF_6D_Modes<REAL>, Vector>::toFirst(vec);
	}
};

template<typename REAL>
class TypesConverter<Result_6D_Modes<REAL>, ObjGrad> {
public:
	static Result_6D_Modes<REAL> toFirst(ObjGrad const& objGrad) {
		Result_6D_Modes<REAL> enGrad;
		enGrad._6D.E = objGrad.obj;
		enGrad._6D.ang.x = objGrad.grad(0);
		enGrad._6D.ang.y = objGrad.grad(1);
		enGrad._6D.ang.z = objGrad.grad(2);
		enGrad._6D.pos.x = objGrad.grad(3);
		enGrad._6D.pos.y = objGrad.grad(4);
		enGrad._6D.pos.z = objGrad.grad(5);
		for(int mode=0;mode< Common_Modes::numModesRec; mode++){
			enGrad.modesRec[mode] = objGrad.grad(6 + mode);
		}
		for(int mode=0;mode< Common_Modes::numModesLig; mode++){
			enGrad.modesLig[mode] = objGrad.grad(6 + mode + Common_Modes::numModesRec);
		}
		return enGrad;
	}

	static ObjGrad toSecond (Result_6D_Modes<REAL> const& enGrad) {
		ObjGrad objGrad;
		objGrad.obj = enGrad._6D.E;
		objGrad.grad = Vector( DOF_MAX_NUMBER );
		objGrad.grad(0) = -enGrad._6D.ang.x;
		objGrad.grad(1) = -enGrad._6D.ang.y;
		objGrad.grad(2) = -enGrad._6D.ang.z;
		objGrad.grad(3) = -enGrad._6D.pos.x;
		objGrad.grad(4) = -enGrad._6D.pos.y;
		objGrad.grad(5) = -enGrad._6D.pos.z;

		for(int mode=0;mode< Common_Modes::numModesRec; mode++){
			objGrad.grad(6 + mode)  = -enGrad.modesRec[mode];
		}
		for(int mode=0;mode< Common_Modes::numModesLig; mode++){
			objGrad.grad(6 + Common_Modes::numModesRec +  mode)  = -enGrad.modesLig[mode];
		}
		return objGrad;
	}
};

template<typename REAL>
class TypesConverter<ObjGrad, Result_6D_Modes<REAL>> {
public:

	static ObjGrad toFirst (Result_6D_Modes<REAL> const& enGrad) {
		return TypesConverter<Result_6D_Modes<REAL>, ObjGrad>::toSecond(enGrad);
	}

	static Result_6D_Modes<REAL> toSecond(ObjGrad const& objGrad) {
		return TypesConverter<Result_6D_Modes<REAL>, ObjGrad>::toFirst(objGrad);
	}
};

} //namespace

#endif /* META_H_ */
