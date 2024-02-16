

#include "libmesh/libmesh_config.h"
#include "TigerFailureCriteria.h"
#include "RankTwoScalarTools.h"
#include "MooseEnum.h"
#include "Executioner.h"
#include "FEProblem.h"
#include <cfloat>
#include "Function.h"
# define PI          3.141592653589793238462643383279502884

registerMooseObject("TigerApp", TigerFailureCriteria);

InputParameters
TigerFailureCriteria::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Requests termination of the current solve based on the evaluation of"
                             " a parsed logical expression of the Postprocessor value(s).");

      params.addParam<std::string>("base_name", "the identical base name provided "
            "in TensorMechanics Action");
      params.addRequiredParam<MaterialPropertyName>("total_stress",
                                                    "The rank two material tensor name");
      params.addParam<Real>("cohesion",0.0,
                                                  "The value of cohesion should be given by user");
      params.addParam<Real>("phi", 35.0,
                                                  "The value of internal friction coefficient should be given by user in degree");
      params.addParam<Real>("mu_s", 0.6,
                                                  "The value of friction coefficient should be given by user");
      params.addParam<MooseEnum>("criterion_type",
                                                                       MooseEnum("Mohr_Coulomb Slip_Tedency", "Mohr_Coulomb", "Fracture_Potential"),
                                                                       "Criterion to use for the threshold");
      params.set<bool>("use_displaced_mesh") = true;
  return params;
}

TigerFailureCriteria::TigerFailureCriteria(const InputParameters & parameters)
  : AuxKernel(parameters),
    _use_displaced_mesh(getParam<bool>("use_displaced_mesh")),
    _cohesion(getParam<Real>("cohesion")),
    _phi(getParam<Real>("phi")),
    _mu_s (getParam<Real>("mu_s")),
    _criterion_type(getParam<MooseEnum>("criterion_type").getEnum<CriterionType>()),
    _normal_vector(getMaterialProperty<RealVectorValue>("normal_vector")),
    _Normal_Traction(getMaterialProperty<RealVectorValue>("Normal_Traction")),
    _sig_n(getMaterialProperty<Real>("NORMAL_STRESS")),
    _sig_t(getMaterialProperty<Real>("SHEAR_STRESS")),
    _tensor(nullptr)
{
  if (hasMaterialProperty<RankTwoTensor>("rank_two_tensor"))
    _tensor = &getMaterialProperty<RankTwoTensor>("rank_two_tensor");
  else
    mooseError("Error in RankTwoBasedFailureCriteriaNOSPD! Required rank two tensor is not "
               "available for current peridynamics model!");

  if (isNodal())
    paramError("variable", "This AuxKernel only supports Elemental fields");
}


Real
TigerFailureCriteria::computeValue()
{
Real ratio = 1.0;
switch (_criterion_type)
{
  case CriterionType::Mohr_Coulomb:
{
  Real _tau_strength = _cohesion - (_sig_n [_qp] * tan (_phi*PI/180));


      ratio =   abs ( _sig_t[_qp]/_tau_strength);

  break;
  }
  case CriterionType::Slip_Tedency:
{
   ratio =   abs (( _sig_t[_qp]/_sig_n[_qp])/_mu_s);

break;
    }
    case CriterionType::Fracture_Potential:
    {
      RankTwoTensor avg_tensor = 0.5 * ((*_tensor)[0] + (*_tensor)[1]);
      Real _max_shear = 2.0 * RankTwoScalarTools::maxShear(avg_tensor);         // calculates 'sigma 1 - sigma 3' in tems of mohr circle
      Real _distance = _sig_n[_qp] * (tan(_phi * PI / 180)/ cos(_phi * PI /180));
      ratio = _max_shear / (_max_shear + _distance);

      break;
    }
  }
  return ratio;
}
