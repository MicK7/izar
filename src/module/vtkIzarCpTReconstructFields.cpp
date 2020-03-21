#include "vtkIzarCpTReconstructFields.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointSet.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkSetGet.h"
#include <iostream>

vtkStandardNewMacro(vtkIzarCpTReconstructFields)

vtkIzarCpTReconstructFields::vtkIzarCpTReconstructFields():
  ReplaceVars(1), AdvancedMode(0), SuthConst(110.4),
  SuthMuRef(1.717e-5), SuthTRef(273.0),
  ComputePs(1), ComputePtabs(1), ComputeTtabs(1),
  ComputePtrel(1), ComputeTtrel(1), ComputeS(1), ComputeV(1),
  ComputeW(1), ComputeM(1), ComputeMrel(1), ComputePhi(1),
  ComputeAlpha(1), ComputeBeta(1), ComputeMu(1), ComputeGamma(1),
  ComputeHtabs(0), ComputeHtrel(0)
{
  IZAR_WARNING
}

vtkIzarCpTReconstructFields::~vtkIzarCpTReconstructFields()
{

}

void vtkIzarCpTReconstructFields::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "vtkIzarCpTReconstructFields";
}

int vtkIzarCpTReconstructFields::RequestData(vtkInformation* request,
                                                   vtkInformationVector** inVector, vtkInformationVector* outVector)
{
  vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
  if(!inInfo)
    {
    return 0;
    }
  vtkPointSet* inData = vtkPointSet::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if(!inData)
    {
    return 0;
    }

  vtkInformation* outInfo = outVector->GetInformationObject(0);
  vtkPointSet* outData = vtkPointSet::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if(!outData)
    {
    return 0;
    }

  outData->ShallowCopy(inData);

  if(outData->GetNumberOfPoints() == 0) { return 1; } // Do nothing

  vtkDataArray* coords = outData->GetPoints()->GetData();
  vtkDataArray* rho = outData->GetPointData()->GetArray(RHO);
  vtkDataArray* rhov = outData->GetPointData()->GetArray(RHOV);
  vtkDataArray* rhoe = outData->GetPointData()->GetArray(RHOE);

  if(!coords) { vtkErrorMacro("Cannot find coordinates array"); return 0; }
  if(coords->GetNumberOfComponents() != 3) { vtkErrorMacro("Cannot find vector coordinates"); return 0; }
  if((!rho) || (rho->GetNumberOfComponents() != 1)) { vtkErrorMacro("Cannot find scalar Density field '" << RHO << "'"); return 0; }
  if((!rhov) || (rhov->GetNumberOfComponents() != 3)) { vtkErrorMacro("Cannot find vector Momentum field '" << RHOV "'"); return 0; }
  if((!rhoe) || (rhoe->GetNumberOfComponents() != 1)) { vtkErrorMacro("Cannot find scalar EnergyStagnationDensity field '" << RHOE << "'"); return 0; }

  int N = coords->GetNumberOfTuples();
  if((rho->GetNumberOfTuples() != N) || (rhov->GetNumberOfTuples() != N) || (rhoe->GetNumberOfTuples() != N))
    {
    vtkErrorMacro("Number of tuples mismatch");
    return 0;
    }

  // This filter could be extended to get rid of this assumption
  if((rho->GetDataType() != rhov->GetDataType()) ||
     (rho->GetDataType() != rhoe->GetDataType()))
    {
    vtkErrorMacro("The Density, Momentum, and EnergyStagnationDensity must have the same data type (eg. double)");
    return 0;
    }

  // Get the field data
  double omega;
  double rgas;
  double rfuel;
  double rwater;
  double far = 0.0;
  double war = 0.0;
  double rmix;

  if(!Izar::GetDoubleFieldData(outData, OMEGA, omega))
    {
    vtkErrorMacro("Cannot find omega in the field data: assuming omega = 0");
    omega = 0.;
    }
  if(!Izar::GetDoubleFieldData(outData, RGAS, rgas))
    {
    vtkErrorMacro("Cannot find rgas in the field data: assuming rgas = 287.05 J/kg/K");
    rgas = 287.05;
    }
  if(!Izar::GetDoubleFieldData(outData, RFUEL, rfuel))
    {
    vtkErrorMacro("Cannot find rfuel in the field data: assuming rfuel = 291.356 J/kg/K");
    rfuel = 291.356;
    }
  if(!Izar::GetDoubleFieldData(outData, RWATER, rwater))
    {
    vtkErrorMacro("Cannot find rwater in the field data: assuming rwater = 461.51 J/kg/K");
    rwater = 461.51;
    }
  if(!Izar::GetDoubleFieldData(outData, FAR, far))
    {
    vtkErrorMacro("Cannot find FAR (fuel air ratio) in the field data: assuming far = 0.0");
    far = 0.0;
    }
  if(!Izar::GetDoubleFieldData(outData, WAR, rwater))
    {
    vtkErrorMacro("Cannot find WAR (water vapor air ratio) in the field data: assuming war = 0.0 (dry air)");
    war = 0.0;
    }
  // Get Polynomials for CP_{air,kero,water}
  double * cpAirPoly;
  double * cpFuelPoly;
  double * cpWaterPoly;
  double * cpMixPoly;

  int sizeAirPoly;
  int sizeFuelPoly;
  int sizeWaterPoly;
  double defaultCpAir[8] = {1005.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if(!Izar::GetDoubleFieldDataArray(inData, CPAIR, &cpAirPoly, sizeAirPoly))
    {
    vtkErrorMacro("Cannot find CPAIR (polynomial for air Cp) in the field data:"
                  " assuming cpAirPoly = [1005. ,0.0, ...]");
    cpAirPoly = defaultCpAir;
    sizeAirPoly = 8;
    }
  if(!Izar::GetDoubleFieldDataArray(outData, CPFUEL, &cpFuelPoly, sizeFuelPoly))
    {
    vtkErrorMacro("Cannot find CPFUEL (polynomial for fuel Cp) in the field data:"
                  " assuming cpFuelPoly = [1005. ,0.0, ...]");
    cpFuelPoly = defaultCpAir;
    sizeFuelPoly = 8;
    }
  if(!Izar::GetDoubleFieldDataArray(outData, CPWATER, &cpWaterPoly, sizeWaterPoly))
    {
    vtkErrorMacro("Cannot find CPWATER (polynomial for water vapor Cp) in the field data:"
                  " assuming cpWaterPoly = [1005. ,0.0, ...]");
    cpWaterPoly = defaultCpAir;
    sizeWaterPoly = 8;
    }

  if ((sizeAirPoly != sizeFuelPoly) ||
      (sizeAirPoly != sizeWaterPoly))
    {
    vtkErrorMacro("Polynomial for air, fuel and water should be of same order.");
    }

  // -- Compute Mixing Gas Constant :
  rmix = (rgas+far*rfuel+war*rwater)/(1.0+far+war);
  // -- Save it as Field Data
  vtkDoubleArray* arRmix = vtkDoubleArray::New();
  arRmix->SetNumberOfComponents(1);
  arRmix->SetNumberOfTuples(1);
  double* data = reinterpret_cast<double*>(arRmix->GetVoidPointer(0));
  data[0] = rmix;
  arRmix->SetName(RMIX);
  outData->GetFieldData()->AddArray(arRmix);
  arRmix->Delete();
  // -- Compute CpPolyMix as Field Data
  vtkDoubleArray* arCpMix = vtkDoubleArray::New();
  arCpMix->SetNumberOfComponents(1);
  arCpMix->SetNumberOfTuples(sizeAirPoly);
  cpMixPoly = reinterpret_cast<double*>(arCpMix->GetVoidPointer(0));
  for (int ii=0; ii<sizeAirPoly; ii++)
    {
    cpMixPoly[ii] = cpAirPoly[ii] + far*cpFuelPoly[ii] + war*cpWaterPoly[ii];
    }
  arCpMix->SetName(CPMIX);
  outData->GetFieldData()->AddArray(arCpMix);
  arCpMix->Delete();
  // ---

  // Generate the new arrays
  vtkDataArray* ts = NULL;
  vtkDataArray* cp = NULL;
  vtkDataArray* ps = NULL;
  vtkDataArray* ptabs = NULL;
  vtkDataArray* ttabs = NULL;
  vtkDataArray* ptrel = NULL;
  vtkDataArray* ttrel = NULL;
  vtkDataArray* s = NULL;
  vtkDataArray* v = NULL;
  vtkDataArray* w = NULL;
  vtkDataArray* M = NULL;
  vtkDataArray* Mrel = NULL;
  vtkDataArray* phi = NULL;
  vtkDataArray* alpha = NULL;
  vtkDataArray* beta = NULL;
  vtkDataArray* mu = NULL;
  vtkDataArray* gamma = NULL;
  vtkDataArray* htabs = NULL;
  vtkDataArray* htrel = NULL;

  // Mandatory:
  ts = rho->NewInstance(); ts->SetNumberOfComponents(1); ts->SetNumberOfTuples(N);
  cp = rho->NewInstance(); cp->SetNumberOfComponents(1); cp->SetNumberOfTuples(N);

  // Optionals
  if(ComputePs) { ps = rho->NewInstance(); ps->SetNumberOfComponents(1); ps->SetNumberOfTuples(N); }
  if(ComputePtabs) { ptabs = rho->NewInstance(); ptabs->SetNumberOfComponents(1); ptabs->SetNumberOfTuples(N); }
  if(ComputeTtabs) { ttabs = rho->NewInstance(); ttabs->SetNumberOfComponents(1); ttabs->SetNumberOfTuples(N); }
  if(ComputePtrel) { ptrel = rho->NewInstance(); ptrel->SetNumberOfComponents(1); ptrel->SetNumberOfTuples(N); }
  if(ComputeTtrel) { ttrel = rho->NewInstance(); ttrel->SetNumberOfComponents(1); ttrel->SetNumberOfTuples(N); }
  if(ComputeS) { s = rho->NewInstance(); s->SetNumberOfComponents(1); s->SetNumberOfTuples(N); }
  if(ComputeV) { v = rho->NewInstance(); v->SetNumberOfComponents(3); v->SetNumberOfTuples(N); }
  if(ComputeW) { w = rho->NewInstance(); w->SetNumberOfComponents(3); w->SetNumberOfTuples(N); }
  if(ComputeM) { M = rho->NewInstance(); M->SetNumberOfComponents(1); M->SetNumberOfTuples(N); }
  if(ComputeMrel) { Mrel = rho->NewInstance(); Mrel->SetNumberOfComponents(1); Mrel->SetNumberOfTuples(N); }
  if(ComputePhi) { phi = rho->NewInstance(); phi->SetNumberOfComponents(1); phi->SetNumberOfTuples(N); }
  if(ComputeAlpha) { alpha = rho->NewInstance(); alpha->SetNumberOfComponents(1); alpha->SetNumberOfTuples(N); }
  if(ComputeBeta) { beta = rho->NewInstance(); beta->SetNumberOfComponents(1); beta->SetNumberOfTuples(N); }
  if(ComputeMu) { mu = rho->NewInstance(); mu->SetNumberOfComponents(1); mu->SetNumberOfTuples(N); }
  if(ComputeGamma) { gamma = rho->NewInstance(); gamma->SetNumberOfComponents(1); gamma->SetNumberOfTuples(N); }
  if(ComputeHtabs) { htabs = rho->NewInstance(); htabs->SetNumberOfComponents(1); htabs->SetNumberOfTuples(N); }
  if(ComputeHtrel) { htrel = rho->NewInstance(); htrel->SetNumberOfComponents(1); htrel->SetNumberOfTuples(N); }

  switch(vtkTemplate2PackMacro(coords->GetDataType(), rho->GetDataType()))
    {
    vtkTemplate2MacroCase2(VTK_DOUBLE, double, VTK_DOUBLE, double,
                                              (this->CallSMPCpTReconstructOp<VTK_T1, VTK_T2>(
                                               coords, rho, rhov, rhoe,
                                               ps, ts, ptabs, ttabs, ptrel, ttrel, s, v, w, M, Mrel, phi, alpha, beta, mu,
                                               htabs, htrel, cp, gamma,
                                               ComputePs, ComputePtabs, ComputeTtabs, ComputePtrel,
                                               ComputeTtrel, ComputeS, ComputeV, ComputeW, ComputeM, ComputeMrel, ComputePhi,
                                               ComputeAlpha, ComputeBeta, ComputeMu,
                                               ComputeGamma, ComputeHtrel, ComputeHtabs,
                                               cpMixPoly, sizeAirPoly,
                                               rmix, omega, N)));
    vtkTemplate2MacroCase2(VTK_FLOAT, float, VTK_FLOAT, float,
                                              (this->CallSMPCpTReconstructOp<VTK_T1, VTK_T2>(
                                               coords, rho, rhov, rhoe,
                                               ps, ts, ptabs, ttabs, ptrel, ttrel, s, v, w, M, Mrel, phi, alpha, beta, mu,
                                               htabs, htrel, cp, gamma,
                                               ComputePs, ComputePtabs, ComputeTtabs, ComputePtrel,
                                               ComputeTtrel, ComputeS, ComputeV, ComputeW, ComputeM, ComputeMrel, ComputePhi,
                                               ComputeAlpha, ComputeBeta, ComputeMu,
                                               ComputeGamma, ComputeHtrel, ComputeHtabs,
                                               cpMixPoly, sizeAirPoly,
                                               rmix, omega, N)));
    default:
      vtkErrorMacro("Should not happen");
      return 0;
      break;
    }

  if(ts) { this->ReplaceVarByArrayIfNeeded(outData, ts, TS); ts->Delete(); }
  if(cp) { this->ReplaceVarByArrayIfNeeded(outData, cp, CP); cp->Delete(); }

  if(ps) { this->ReplaceVarByArrayIfNeeded(outData, ps, PS); ps->Delete(); }
  if(ptabs) { this->ReplaceVarByArrayIfNeeded(outData, ptabs, PTABS); ptabs->Delete(); }
  if(ttabs) { this->ReplaceVarByArrayIfNeeded(outData, ttabs, TTABS); ttabs->Delete(); }
  if(ptrel) { this->ReplaceVarByArrayIfNeeded(outData, ptrel, PTREL); ptrel->Delete(); }
  if(ttrel) { this->ReplaceVarByArrayIfNeeded(outData, ttrel, TTREL); ttrel->Delete(); }
  if(s) { this->ReplaceVarByArrayIfNeeded(outData, s, ENTROPY); s->Delete(); }
  if(v) { this->ReplaceVarByArrayIfNeeded(outData, v, VELOCITYABS); v->Delete(); }
  if(w) { this->ReplaceVarByArrayIfNeeded(outData, w, VELOCITYREL); w->Delete(); }
  if(M) { this->ReplaceVarByArrayIfNeeded(outData, M, MACHABS); M->Delete(); }
  if(Mrel) { this->ReplaceVarByArrayIfNeeded(outData, Mrel, MACHREL); Mrel->Delete(); }
  if(phi) { this->ReplaceVarByArrayIfNeeded(outData, phi, PHI); phi->Delete(); }
  if(alpha) { this->ReplaceVarByArrayIfNeeded(outData, alpha, ALPHA); alpha->Delete(); }
  if(beta) { this->ReplaceVarByArrayIfNeeded(outData, beta, BETA); beta->Delete(); }
  if(mu) { this->ReplaceVarByArrayIfNeeded(outData, mu, MU); mu->Delete(); }

  if(gamma) { this->ReplaceVarByArrayIfNeeded(outData, gamma, GAMMA); gamma->Delete(); }
  if(htabs) { this->ReplaceVarByArrayIfNeeded(outData, htabs, HTABS); htabs->Delete(); }
  if(htrel) { this->ReplaceVarByArrayIfNeeded(outData, htrel, HTREL); htrel->Delete(); }

  return 1;
}

void vtkIzarCpTReconstructFields::ReplaceVarByArrayIfNeeded(vtkDataSet* data, vtkDataArray* ar, const char* name)
{
  if(data->GetPointData()->HasArray(name))
    {
    if(!this->ReplaceVars) { return; }
    data->GetPointData()->RemoveArray(name);
    }
  ar->SetName(name);
  data->GetPointData()->AddArray(ar);
}
