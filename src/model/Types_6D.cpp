/*
 * DOF_6D.cpp
 *
 *  Created on: Aug 24, 2016
 *      Author: uwe
 */

#include "Types_6D.tpp"

using namespace std;

namespace as {

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D<float> const& dof);

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D<double> const& dof);


template
std::ostream& operator<< (std::ostream& s, Result_6D<float> const& enGrad);

template
std::ostream& operator<< (std::ostream& s, Result_6D<double> const& enGrad);

Vec3<double> Common::pivotRec = Vec3<double>(0,0,0);
Vec3<double> Common::pivotLig = Vec3<double>(0,0,0);
bool Common::centeredRec = false;
bool Common::centeredLig = false;



}  // namespace as


