/*
 * Types_6D.tpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#ifndef SRC_TYPES_6D_MODES_TPP_
#define SRC_TYPES_6D_MODES_TPP_

#include <ostream>
#include <iomanip>
#include <sstream>
#include "Types_6D_Modes.h"
#include <iostream>
namespace as {


template<typename REAL>
std::ostream& operator<<(std::ostream& outStream, DOF_6D_Modes<REAL> const& dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;


	outStream 	<< setw(w) << 0 << setw(w) << 0 << setw(w) << 0
				<< setw(w) << 0 << setw(w) << 0 << setw(w) << 0;
				for(int mode=0;mode<Common_Modes::numModesRec;mode++){
					outStream<< setw(w) << dof.modesRec[mode];
				}

	outStream   << std::endl << setw(w) << dof._6D.ang.x << setw(w) << dof._6D.ang.y << setw(w) << dof._6D.ang.z
			                 << setw(w) << dof._6D.pos.x + Common_Modes::pivotRec.x - Common_Modes::pivotLig.x 
                             << setw(w) << dof._6D.pos.y + Common_Modes::pivotRec.y - Common_Modes::pivotLig.y
                             << setw(w) << dof._6D.pos.z + Common_Modes::pivotRec.z - Common_Modes::pivotLig.z;
				for(int mode=0;mode<Common_Modes::numModesLig;mode++){
					outStream<< setw(w) << dof.modesLig[mode];
				}


	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template<typename REAL>
std::ostream& operator<<(std::ostream& s, Result_6D_Modes<REAL> const& enGrad) {
	using namespace std;
	int precisionSetting = s.precision( );
	ios::fmtflags flagSettings = s.flags();

	s.setf(ios::scientific);
	s.precision(8);
	s << " Energy: " << enGrad._6D.E << endl;
	s.unsetf(ios::scientific);

	s.setf(ios::fixed);
	s.precision(3);
	s << setw(12) << enGrad._6D.E << endl;
	s.unsetf(ios::fixed);

	s.setf(ios::scientific);
	s.precision(8);
	int width = 20;


	s 	<< setw(width) << 0  << setw(width) << 0  << setw(width) << 0
		<< setw(width) << 0  << setw(width) << 0  << setw(width) << 0;
		for(int mode=0;mode<Common_Modes::numModesRec;mode++){
			s<< setw(width) << enGrad.modesRec[mode];
		}

	s 	<< std::endl 
        << setw(width) << enGrad._6D.ang.x  << setw(width) << enGrad._6D.ang.y  << setw(width) << enGrad._6D.ang.z
		<< setw(width) << enGrad._6D.pos.x  << setw(width) << enGrad._6D.pos.y  << setw(width) << enGrad._6D.pos.z;
		for(int mode=0;mode<Common_Modes::numModesLig;mode++){
			s<< setw(width) << enGrad.modesLig[mode];
		}

	s.unsetf(ios::scientific);

	s.precision(precisionSetting);
	s.flags(flagSettings);

	return s;
}

template<typename dof_t, typename result_t>
void printResults(dof_t dof, result_t result, int index)
{
    std::cout << "#"<<index<< std::endl;
    std::cout << "## Energy: "<<result._6D.E<< std::endl;
    std::cout << dof << std::endl;
}

} // namespace as



#endif /* SRC_TYPES_6D_TPP_ */
