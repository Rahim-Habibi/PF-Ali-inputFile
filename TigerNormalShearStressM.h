#pragma once
#include "Material.h"
#include "libmesh/quadrature_gauss.h"
#include "RankTwoTensor.h"

class TigerNormalShearStressM : public Material
{
public:
  static InputParameters validParams();

TigerNormalShearStressM(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// For 1D elements, the value computed will be perpendicular to this vector
  const RealVectorValue _1D_perp;
  /// Value used for 3D elements
  const RealVectorValue _3D_default;
  /// a function to calculate the components of the normal vector
  RealVectorValue normal_fun ();
  /// normal vector acting on 2D elements
  MaterialProperty <RealVectorValue> & _normal_vector;
const bool _use_displaced_mesh;
 std::string _base_name;
 /// total stress in the form of tensor
const MaterialProperty<RankTwoTensor> & _TenMech_total_stress;
/// a function to calculate the components of the normal traction
 RealVectorValue  normal_traction();
 /// a vector to store the components of the normal traction
MaterialProperty <RealVectorValue> & _Normal_Traction;
/// normal stress perpendicular on the element
MaterialProperty <Real> & _sig_n;
/// shear stress parallel to the element
MaterialProperty<Real>  &  _sig_t;
private:
  // enum to select whether shear stress calculated or not?
  bool _DoUWantShear;
};
