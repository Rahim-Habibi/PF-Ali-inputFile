
#pragma once

#include "libmesh/libmesh_config.h"
#include "AuxKernel.h"
#include "libmesh/fparser.hh"
#include "RankTwoTensor.h"

class TigerFailureCriteria : public AuxKernel
{
public:
  static InputParameters validParams();

  TigerFailureCriteria(const InputParameters & parameters);


protected:
  //void handleMessage();
virtual Real computeValue() override;


  const bool _use_displaced_mesh;
  /// cohesion of the material in Pa
 Real  _cohesion;
 /// a function to calculate normal stress
 Real  sigma_nn ();
 /// friction coefficient of material in degree
 Real  _phi;
 Real  _mu_s;
 /// Normal vector acting on 2D elements specially on fault
 const MaterialProperty<RealVectorValue> & _normal_vector;
 /// normal traction
 const MaterialProperty<RealVectorValue>  &  _Normal_Traction;
 const MaterialProperty<RankTwoTensor> * _tensor;
 /// normal stress
const MaterialProperty<Real>  &  _sig_n;
   /// shear stress
  const MaterialProperty<Real>  &  _sig_t;
private:
// enum to select failure criterion
const enum class CriterionType { Mohr_Coulomb, Slip_Tedency, Fracture_Potential } _criterion_type;
};
