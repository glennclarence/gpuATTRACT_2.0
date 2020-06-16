#include "Types_6D_Modes.tpp"

using namespace std;

namespace as {

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D_Modes<float> const& dof);

template
std::ostream& operator <<(std::ostream& outStream, DOF_6D_Modes<double> const& dof);


template
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<float> const& enGrad);

template
std::ostream& operator<< (std::ostream& s, Result_6D_Modes<double> const& enGrad);

template
void printResults<DOF_6D_Modes<float>, Result_6D_Modes<float>>(DOF_6D_Modes<float> dof, Result_6D_Modes<float> result, int index);

template
void printResults<DOF_6D_Modes<double>, Result_6D_Modes<double>>(DOF_6D_Modes<double> dof, Result_6D_Modes<double> result, int index);

unsigned int Common_Modes::numModesRec = 0;
unsigned int Common_Modes::numModesLig = 0;

}  // namespace as

