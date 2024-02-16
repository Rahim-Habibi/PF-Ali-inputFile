

#include "TigerNormalShearStressM.h"

registerMooseObject("TigerApp", TigerNormalShearStressM);

InputParameters
TigerNormalShearStressM::validParams()
{
  InputParameters params = Material::validParams();

  params.addParam<RealVectorValue>("3D_default",
                                   RealVectorValue(0, 0, 1),
                                   "The value that will be produced for 3D elements, since such "
                                   "elements do not have a 'normal direction'");
  params.addParam<RealVectorValue>(
      "1D_perp",
      RealVectorValue(0, 0, 1),
      "The normal for all 1D elements will be perpendicular to this vector");

      params.addParam<std::string>("base_name", "the identical base name provided "
            "in TensorMechanics Action");
      params.addRequiredParam<MaterialPropertyName>("total_stress",
                                                    "The rank two material tensor name");
       params.addRequiredParam<bool>("DoUWantShear",
            "The first calculates just normal stress but the second calculates normal and shear required for TigerTerminator [JustNormal, NormalAndShear]");
      params.set<bool>("use_displaced_mesh") = true;
  params.addClassDescription(
      "Material to compute Normal and shear stress on the 2D element.  This is mostly designed for 2D "
      "elements living in 3D space, however, the 1D-element and 3D-element cases are set  "
      "zero intentionally.  The Variable for this Material must be an elemental Variable");

  return params;
}

TigerNormalShearStressM::TigerNormalShearStressM(const InputParameters & parameters)
  : Material(parameters),
    _1D_perp(getParam<RealVectorValue>("1D_perp")),
    _3D_default(getParam<RealVectorValue>("3D_default")),
    _use_displaced_mesh(getParam<bool>("use_displaced_mesh")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Normal_Traction(declareProperty<RealVectorValue>("Normal_Traction")),
    _TenMech_total_stress(getMaterialProperty<RankTwoTensor>(_base_name + "total_stress")),
    _DoUWantShear(getParam<bool>("DoUWantShear")),
    _sig_n(declareProperty<Real>("NORMAL_STRESS")),
    _sig_t(declareProperty<Real>("SHEAR_STRESS")),
    _normal_vector(declareProperty<RealVectorValue>("normal_vector"))
{
  if (isNodal())
    paramError("variable", "The variable must be an elemental variable");
  if (_1D_perp.norm() == 0.0)
    paramError("1D_perp", "Must not be the zero vector");
  if (_3D_default.norm() == 0.0)
    paramError("3D_default", "Must not be the zero vector");
}

void
TigerNormalShearStressM::computeQpProperties()
{
  _normal_vector [_qp] = normal_fun(); //component
  _Normal_Traction [_qp] = (_TenMech_total_stress[_qp]) * _normal_vector[_qp];
  /// a function to calculate power of normal traction
  Real traction2 = pow (_Normal_Traction[_qp](0) ,2) +
               pow (_Normal_Traction[_qp](1) ,2) +
               pow (_Normal_Traction[_qp](2) ,2);
  /// a function to calculate normal stress in this scope
  Real  sigma_n =  (_Normal_Traction[_qp](0) * _normal_vector[_qp](0)) +
                        (_Normal_Traction[_qp](1) * _normal_vector[_qp](1)) +
                        (_Normal_Traction[_qp](2) * _normal_vector[_qp](2));

   if (!_DoUWantShear)

      _sig_n [_qp] = (_Normal_Traction[_qp](0) * _normal_vector[_qp](0)) +
      (_Normal_Traction[_qp](1) * _normal_vector[_qp](1)) +
      (_Normal_Traction[_qp](2) * _normal_vector[_qp](2));

    if (_DoUWantShear)
    //This line calculates normal stress only. If you wish to have both normal and shear values output, you should keep it activated,
    //however if you only intreseted in shear stress and failure criterion, you can deactivate this line.
      _sig_n [_qp] =  (_Normal_Traction[_qp](0) * _normal_vector[_qp](0)) +
             (_Normal_Traction[_qp](1) * _normal_vector[_qp](1)) +
             (_Normal_Traction[_qp](2) * _normal_vector[_qp](2));
     /// calculating shear stress using normal stress and traction components
       _sig_t [_qp] =  abs ( sqrt (traction2 - sigma_n * sigma_n));

}

  RealVectorValue         // if you define it as Real, it will print only one component.
  TigerNormalShearStressM::normal_fun()
  {
  RealVectorValue _n;
  const auto num_nodes = _current_elem->n_nodes();
  switch (_current_elem->dim())
  {
    case 1:
    {
      for (unsigned i = 0; i < num_nodes - 1; ++i)
      {
        RealVectorValue v = _current_elem->point((i + 1) % num_nodes) - _current_elem->point(i);
        _n += v.cross(_1D_perp);
      }
      break;
    }
    case 2:
    {
      for (unsigned i = 0; i < num_nodes - 2; ++i)
      {
        RealVectorValue v1 = _current_elem->point((i + 1) % num_nodes) - _current_elem->point(i);
        RealVectorValue v2 =
            _current_elem->point((i + 2) % num_nodes) - _current_elem->point((i + 1) % num_nodes);
        _n += v1.cross(v2);
      }
      break;
    }
    default:
      _n = _3D_default;
  }
   return _n.unit(); //component
}
