#ifndef PARAMTABLE_H_
#define PARAMTABLE_H_

#include <stdexcept>
#include <type_traits>

#include "DataItem.h"

namespace as {

template<typename REAL>
class ParamTable : public DataItem {
	// Check if REAL is of floating-point type
	using real_t = typename std::enable_if<std::is_floating_point<REAL>::value, REAL>::type;
public:

	struct attractFFParams_t {
		real_t rc; // A = rc
		real_t ac; // B = ac
		int ipon;
		real_t rmin2;
		real_t emin;
	};

	/* Potential shape type of the ATTRACT force field */
	enum PotShape {
		_12_6,
		_8_6,
		undefined
	};

	using type_t = attractFFParams_t;

	/* Constructor */
	ParamTable() : _paramTable(nullptr), _numTypes(0),
			_shape(undefined), _swiOn(0), _swiOff(0) {}



	/* Destructor */
	~ParamTable() {
		delete[] _paramTable;
	}

	/***************
	* G E T T E R
	***************/
	unsigned numTypes() const noexcept {
		return _numTypes;
	}

	/* read only access */
	const type_t* table() const noexcept {
		return _paramTable;
	}

	PotShape potShape() const noexcept {
		return _shape;
	}

	real_t swiOn() const noexcept {
		return _swiOn;
	}

	real_t swiOff() const noexcept {
		return _swiOff;
	}

	void setNumTypes(unsigned numTypes) noexcept {
		_numTypes = numTypes;
	}

	void setPotShape(PotShape shape) noexcept {
		_shape = shape;
	}

	void setSwiOn(real_t swiOn) noexcept {
		_swiOn = swiOn;
	}

	void setSwiOff(real_t swiOff) noexcept {
		_swiOff = swiOff;
	}

	/****************************
	 * public member functions
	 ****************************/
	inline const type_t& getParams(const int& typeA, const int& typeB) const noexcept {
		return _paramTable[_numTypes*typeA + typeB];
	}

	/*
	 * Read and write access.
	 * Should be used for initialization
	 */
	type_t* getOrCreateTable() {
		if (_paramTable == nullptr) {
			if (_numTypes == 0) {
				throw std::runtime_error("Error: getOrCreateTypePtr(): the number of types must be set before");
			}
			_paramTable = new ParamTable::type_t[_numTypes * _numTypes];
		}
		return _paramTable;
	}

private:

	type_t* _paramTable;
	unsigned _numTypes; /** number of particle/atom types */

	PotShape _shape; /** potential shape 12:6 or 8:6 or undefined */

	real_t _swiOn;	/** switching potential parameter. Not used at the moment */
	real_t _swiOff;
};

}


#endif /* PARAMTABLE_H_ */
