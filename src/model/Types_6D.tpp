/*
 * Types_6D.tpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_6D_TPP_
#define SRC_TYPES_6D_TPP_

#include <ostream>
#include <iomanip>
#include <iostream>
#include "Types_6D.h"

namespace as {

void Common::printDofHeader()
{
    std::cout << "#pivot 1 "<< " "<< pivotRec.x << " "<< pivotRec.y << " "<< pivotRec.z << std::endl;
    std::cout << "#pivot 2 "<< " "<< pivotLig.x << " "<< pivotLig.y << " "<< pivotLig.z << std::endl;
    std::cout << "#centered receptor: "<< centeredRec << std::endl;
    std::cout << "#centered ligand: "<<  centeredLig << std::endl;
}

void Common::setPivotRec(Vec3<double> pivot, bool centered)
{
    pivotRec = pivot;
    centeredRec = centered;
}

void Common::setPivotLig(Vec3<double> pivot, bool centered)
{
    pivotLig = pivot;
    centeredLig = centered;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_6D<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;
    outStream 	<< setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << std::endl;
    outStream   << setw(w) << dof.pos.x + Common::pivotRec.x - Common::pivotLig.x 
                << setw(w) << dof.pos.y + Common::pivotRec.y - Common::pivotLig.y
                << setw(w) << dof.pos.z + Common::pivotRec.z - Common::pivotLig.z;
	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result_6D<REAL> const& enGrad) {
	using namespace std;
	int precisionSetting = s.precision( );
	ios::fmtflags flagSettings = s.flags();

	s.setf(ios::scientific);
	s.precision(8);
	s << " Energy: " << enGrad.E << endl;
	s.unsetf(ios::scientific);

	s.setf(ios::fixed);
	s.precision(3);
	s << setw(12) << enGrad.E << endl;
	s.unsetf(ios::fixed);

	s.setf(ios::scientific);
	s.precision(8);
	int width = 20;
    s   << setw(width) << 0 << setw(width) << 0 << setw(width) << 0 << setw(width) << 0 << setw(width) << 0 << setw(width) << 0 << setw(width) << std::endl;
    s	<< setw(width) << enGrad.ang.x  << setw(width) << enGrad.ang.y  << setw(width) << enGrad.ang.z
        << setw(width) << enGrad.pos.x  << setw(width) << enGrad.pos.y  << setw(width) << enGrad.pos.z;
	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}

template<typename dof_t, typename result_t>
void printResults(dof_t dof, result_t result, int index)
{
    std::cout << "#"<<index<< std::endl;
    std::cout << "## Energy: "<<result.E<< std::endl;
    std::cout << dof << std::endl;
}

} // namespace as



#endif /* SRC_TYPES_6D_TPP_ */
