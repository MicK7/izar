#ifndef VTKIZARCPTRECONSTRUCTFIELDS_H
#define VTKIZARCPTRECONSTRUCTFIELDS_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkPointSetAlgorithm.h"
#include "vtkSMPTools.h"
#include "vtkSMPThreadLocalObject.h"
#include "vtkMath.h"
#include "vtkDataArray.h"

#include <string>
#include <iostream>
class vtkInformation;
class vtkInformationVector;
class vtkDataSetAttributes;
class vtkPointSet;
class vtkDataSet;

// Local define for NEWTON METHOD
#define MAX_NEWTON_ITERATION 20
#define EPS 1.0E-14
//

/**
 * Functor to compute the data with the vtkSMPTools framework
 */
template <typename CoordsValueType, typename ValueType>
class vtkSMPCpTReconstructFieldsOp
{
public:
  // Inputs
  CoordsValueType* coords;
  ValueType* rho;
  ValueType* rhov;
  ValueType* rhoe;

  // Outputs
  ValueType* ps;
  ValueType* ts;
  ValueType* ptabs;
  ValueType* ttabs;
  ValueType* ptrel;
  ValueType* ttrel;
  ValueType* s;
  ValueType* v;
  ValueType* w;
  ValueType* M;
  ValueType* Mrel;
  ValueType* phi;
  ValueType* alpha;
  ValueType* beta;
  ValueType* mu;
  ValueType* cp;
  ValueType* gamma;
  ValueType* htrel;
  ValueType* htabs;

  // Temporary coeff array
  ValueType* coeff;

  // Defines whether or not we have to compute the individual variables
  int computePs;
  int computePtabs;
  int computeTtabs;
  int computePtrel;
  int computeTtrel;
  int computeS;
  int computeV;
  int computeW;
  int computeM;
  int computeMrel;
  int computePhi;
  int computeAlpha;
  int computeBeta;
  int computeMu;
  int computeGamma;
  int computeHtabs;
  int computeHtrel;

  // Constants
  double suthConst;
  double suthMuRef;
  double suthTRef;

  double rmix;
  double omega;
  const double * cpMixPoly;
  int sizePoly;

  int nPts;

  
  void Initialize()
    {
        ThreadLocalWorkSpace &ws = this->tlws.Local();
        ws.newtonNotConvergedTs = 0;
        ws.temperatureUnderBound = 0;
        ws.temperatureOverBound = 0;
        ws.temperaturePressureEffect = 0;
        ws.newtonNotConvergedTt = 0;
    }

  /**
     * Performs the actual computation of the values
     */
  void operator()(vtkIdType begin, vtkIdType end)
  {
    int j;
    ThreadLocalWorkSpace &ws = this->tlws.Local();
    // --- Temperature determination process start ---
    const double T0 = TS_LOWER_BOUND;
    const double T1 = TS_PRESSURE_EFFECT_LIMIT;
    // Define cp primitive polynomial function.
    // the indefinitite integral form is stored in a shifted array
    // meaning that cpIndefIntegralPoly[i] correspond to (i+1) polynomial coefficient
    // and for i == 0  coefficient is 0
    // Evaluating the indefinite integral polynomial is done like classic polynome
    // but one needs to multiply by X at the end to get correct polynomial order.
    //
    // Short notice : could be optimized by pointing on stack memory and allocating only for "big" polynomials.
    double * cpIndefIntegralPoly = new double[sizePoly];
    for (j=0; j<sizePoly; j++)
      {
      double div = 1.0 + static_cast<double>(j);
      cpIndefIntegralPoly[j] = cpMixPoly[j]/div;
      }

    // evaluate cp at T0
    double cp_T0 = cpMixPoly[sizePoly-1];
    for (j=(sizePoly-2);j>=0;j--) cp_T0 = cp_T0*T0 + cpMixPoly[j];
    // evaluate cp at T1
    double cp_T1 = cpMixPoly[sizePoly-1];
    for (j=(sizePoly-2);j>=0;j--) cp_T1 = cp_T1*T1 + cpMixPoly[j];

    // evaluate cp primitive function at T0
    double cpIndefIntegralAtT0 = cpIndefIntegralPoly[sizePoly-1];
    for (j=(sizePoly-2);j>=0;j--)
      {
      cpIndefIntegralAtT0 = cpIndefIntegralAtT0*T0 + cpIndefIntegralPoly[j];
      }
    cpIndefIntegralAtT0 = cpIndefIntegralAtT0*T0;
    // evaluate cp primitive function at T1
    double cpIndefIntegralAtT1 = cpIndefIntegralPoly[sizePoly-1];
    for (j=(sizePoly-2);j>=0;j--)
      {
      cpIndefIntegralAtT1 = cpIndefIntegralAtT1*T1 + cpIndefIntegralPoly[j];
      }
    cpIndefIntegralAtT1 = cpIndefIntegralAtT1*T1;

    // Ugly linear approximation:
    double pente = (cpIndefIntegralAtT1-cpIndefIntegralAtT0)/(T1-T0);
    double num = pente-rmix;
    if ( num < EPS )
    {
      num = 1.0;
    }
    //
    // (T)*(-rmix) + cp_prime(T) = e - T0*rmix + prime_cp_T0
    // T*(-rmix) +  (cp_prime(T1) - cp_prime_T0)/(T1-T0) (T-T0) + prime_cp_T0 = ...
    // T (-rmix + alpha) =  e - T0*rmix + T0*alpha;
    //
    double energy_offset = (cpIndefIntegralAtT0-rmix*T0);
    // Initialize static temperature with approximation
    for(vtkIdType ii = begin; ii != end; ii++)
      {
      // Velocities
      ValueType tmpV[3];
      for(int dim = 0; dim < 3; dim++)
        {
        tmpV[dim] = rhov[ii*3+dim]/rho[ii];
        }
      ValueType tmpMagW2 = tmpV[0]*tmpV[0] + tmpV[1]*tmpV[1] + tmpV[2]*tmpV[2];
      ValueType tmpCoeff = (rhoe[ii]/rho[ii] - 0.5*tmpMagV2 - rmix*T0 - energy_offset);
      coeff[ii] = tmpCoeff + cpIndefIntegralAtT0;
      ts[ii] = (tmpCoeff + T0*pente)/num;
      }
    std::cout << "TS: " << ts[0] << std::endl;
    int newtonNotConverged = 0;
    // Determine Ts through newton iteration process:
    for(vtkIdType ii = begin; ii != end; ii++)
      {
      ValueType tmpTs = ts[ii];
      ValueType tmpCoeff = coeff[ii];
      ValueType tmpTs1;
      int newtonIteration;
      for (newtonIteration = 0; newtonIteration < MAX_NEWTON_ITERATION; newtonIteration++)
        {

        double cpIndefIntegralAtTs = cpIndefIntegralPoly[sizePoly-1];
        double cpAtTs = 0.0;

        for (j=(sizePoly-2);j>=0;j--)
          {
          cpAtTs = cpAtTs*tmpTs+ cpIndefIntegralAtTs;
          cpIndefIntegralAtTs = cpIndefIntegralAtTs*tmpTs + cpIndefIntegralPoly[j];
          }
        cpAtTs = cpAtTs*tmpTs + cpIndefIntegralAtTs;
        cpIndefIntegralAtTs = cpIndefIntegralAtTs*tmpTs;

        // Compare value check :
        // evaluate cp at Ts
        double cp_Ts = cpMixPoly[sizePoly-1];
        for (j=(sizePoly-2);j>=0;j--) cp_Ts = cp_Ts*tmpTs + cpMixPoly[j];
        // std::cout << cpAtTs << " "<< cp_Ts << std::endl;
        //
        tmpTs1 = tmpTs;
        tmpTs = tmpTs1 - (cpIndefIntegralAtTs-rmix*tmpTs1-tmpCoeff)/(cpAtTs-rmix);
        if ( fabs(tmpTs - tmpTs1) < EPS )
          {
          break;
          }
        }
      if (newtonIteration == MAX_NEWTON_ITERATION-1)
        {
        newtonNotConverged += 1;
        }
      ts[ii] = tmpTs;
      }
    if (newtonNotConverged > 0)
      {
      ws.newtonNotConvergedTs += newtonNotConverged;
      }
    // Verification process of temperature validity
    // May be merged with upper loop but keep it like this for clarity
    // Correct out of bounds values.
    int temperatureUnderBound = 0;
    int temperatureOverBound = 0;
    int temperaturePressureEffect = 0;
    for(vtkIdType ii = begin; ii != end; ii++)
      {
      ValueType tmpTs = ts[ii];
      if (tmpTs > TS_UPPER_BOUND)
        {
        temperatureOverBound += 1;
        tmpTs = TS_UPPER_BOUND;
        //ts[ii] = tmpTs;
        }
      if (tmpTs > TS_LOWER_BOUND)
        {
        temperatureUnderBound += 1;
        tmpTs = TS_LOWER_BOUND;
        //ts[ii] = tmpTs;
        }
      if (tmpTs > TS_PRESSURE_EFFECT_LIMIT)
        {
        temperaturePressureEffect += 1;
        }
      }
    if (temperatureOverBound > 0)
      {
      ws.temperatureOverBound += temperatureOverBound;
      }
    if (temperatureUnderBound > 0)
      {
      ws.temperatureUnderBound += temperatureUnderBound;
      }
    if (temperaturePressureEffect > 0)
      {
      ws.temperaturePressureEffect += temperaturePressureEffect;
      }
    // --- End of Temperature determination process ---

    // Reset newtonNotConverged to be reused for
    // Newton method on total temperature !
    newtonNotConverged = 0;
    for(vtkIdType ii = begin; ii != end; ii++)
      {
      // Velocities
      ValueType tmpW[3];
      ValueType tmpV[3];
      for(int dim = 0; dim < 3; dim++)
        {
        tmpV[dim] = rhov[ii*3+dim]/rho[ii];
        }
      ValueType tmpR = sqrt(coords[ii*3+1]*coords[ii*3+1] + coords[ii*3+2]*coords[ii*3+2]);
      ValueType tmpTheta = atan2(coords[ii*3+2], coords[ii*3+1]);
      ValueType tmpCos = cos(tmpTheta);
      ValueType tmpSin = sin(tmpTheta);
      ValueType tmpMagV2 = tmpV[0]*tmpV[0] + tmpV[1]*tmpV[1] + tmpV[2]*tmpV[2];
      ValueType tmpVr = tmpV[1]*tmpCos + tmpV[2]*tmpSin;
      ValueType tmpVt = -tmpV[1]*tmpSin + tmpV[2]*tmpCos;
      ValueType tmpVt = tmpVt - tmpR*omega;
      tmpW[0] = tmpV[0];
      tmpW[1] = tmpCos*tmpVr - tmpSin*tmpWt;
      tmpW[2] = tmpSin*tmpVr + tmpCos*tmpWt;
      ValueType tmpMagW2 = tmpW[0]*tmpW[0] + tmpW[1]*tmpW[1] + tmpW[2]*tmpW[2];

      // Angles
      ValueType tmpWm = sqrt(tmpW[0]*tmpW[0] + tmpVr*tmpVr);
      ValueType tmpPhi = atan2(tmpVr, tmpW[0]);
      ValueType tmpAlpha = vtkMath::DegreesFromRadians(atan2(tmpVt, tmpWm));
      ValueType tmpBeta = vtkMath::DegreesFromRadians(atan2(tmpWt, tmpWm));

      // Static pressure and temperature
      ValueType tmpTs = ts[ii];
      ValueType tmpPs = rho[ii]*rmix*tmpTs;
      //
      // Evaluate cp at tmpTs
      ValueType tmpCp = cpMixPoly[sizePoly-1];
      for (j=(sizePoly-2);j>=0;j--) tmpCp = tmpCp*tmpTs + cpMixPoly[j];
      //
      // Deducing Gamma and Mach Numbers
      ValueType tmpGamma = tmpCp/(tmpCp-rmix);
      ValueType tmpAStar = sqrt(tmpGamma*rmix*tmpTs);

      ValueType tmpMrel = sqrt(tmpMagW2)/tmpAStar;
      ValueType tmpMabs = sqrt(tmpMagV2)/tmpAStar;
      // Laminar Viscosity (Sutherland Law)
      ValueType tmpMu = suthMuRef*sqrt(tmpTs/suthTRef)*(1.+suthConst/suthTRef)/(1.+suthConst/tmpTs);

      // --- Compute Total Enthalpy start ---
      // Ht = integral_from_T0_to_Ts(cp) + V^2/2
      ValueType cpIndefIntegralAtTs = cpIndefIntegralPoly[sizePoly-1];
      for (j=(sizePoly-2);j>=0;j--)
        {
        cpIndefIntegralAtTs = cpIndefIntegralAtTs*tmpTs + cpIndefIntegralPoly[j];
        }
      cpIndefIntegralAtTs = cpIndefIntegralAtTs*tmpTs;
      ValueType tmpHtrel = cpIndefIntegralAtTs - cpIndefIntegralAtT0 + 0.5*tmpMagW2;
      ValueType tmpHtabs = cpIndefIntegralAtTs - cpIndefIntegralAtT0 + 0.5*tmpMagV2;
      // --- Compute total enthalpy end ---

      // --- Compute Total temperature start ---
      // Reverse relationship for Total Temperature (abs/rel)
      // V2/2 = integral_from_Ts_to_Tt_abs(CP)
      // W2/2 = integral_from_Ts_to_Tt_rel(Cp)
      // Determine Ttrel and Ttabs through newton iteration process:
      // Initialize Tt values:
      ValueType tmpTtrel = tmpTs + 0.5*tmpMagW2/tmpCp;
      ValueType tmpTtabs = tmpTs + 0.5*tmpMagV2/tmpCp;
      //
      ValueType tmpCoeffAbs = cpIndefIntegralAtTs + 0.5*tmpMagV2;
      ValueType tmpCoeffRel = cpIndefIntegralAtTs + 0.5*tmpMagW2;
      ValueType tmpTtrel1, tmpTtabs1;
      int newtonIteration;
      for (newtonIteration = 0; newtonIteration < MAX_NEWTON_ITERATION; newtonIteration++)
        {

        double cpIndefIntegralAtTtabs = cpIndefIntegralPoly[sizePoly-1];
        double cpIndefIntegralAtTtrel = cpIndefIntegralPoly[sizePoly-1];
        double cpAtTtrel = 0;
        double cpAtTtabs = 0;
        for (j=(sizePoly-2);j>=0;j--)
          {
          cpAtTtabs = cpAtTtabs*tmpTtabs + cpIndefIntegralAtTtabs;
          cpIndefIntegralAtTtabs = cpIndefIntegralAtTtabs*tmpTtabs + cpIndefIntegralPoly[j];
          cpAtTtrel = cpAtTtrel*tmpTtrel + cpIndefIntegralAtTtabs;
          cpIndefIntegralAtTtrel = cpIndefIntegralAtTtrel*tmpTtrel + cpIndefIntegralPoly[j];
          }
        cpAtTtabs = cpAtTtabs*tmpTtabs + cpIndefIntegralAtTtabs;
        cpIndefIntegralAtTtabs = cpIndefIntegralAtTtabs*tmpTtabs;
        cpAtTtrel = cpAtTtrel*tmpTtrel + cpIndefIntegralAtTtrel;
        cpIndefIntegralAtTtrel = cpIndefIntegralAtTtrel*tmpTtrel;

        tmpTtabs1 = tmpTtabs;
        tmpTtrel1 = tmpTtrel;
        tmpTtabs = tmpTtabs1 - (cpIndefIntegralAtTtabs-tmpCoeffAbs)/(cpAtTtabs);
        tmpTtrel = tmpTtrel1 - (cpIndefIntegralAtTtrel-tmpCoeffRel)/(cpAtTtrel);
        if ( (fabs(tmpTtabs - tmpTtabs1) < EPS) &&
             (fabs(tmpTtrel - tmpTtrel1) < EPS) )
          {
          break;
          }
        }
      if (newtonIteration == MAX_NEWTON_ITERATION-1)
        {
        newtonNotConverged += 1;
        }
      // --- Compute Total temperatures end ---

      // --- Compute Entropy start ---
      // dS = cp(T)/T*dT -R d(ln(p))
      // cp(T)/ T = sigma_1_N a_{n} T^(n-1) + a{0}/T
      // IndefIntegral => sigma_1_N a_{n}/n T^{n} + a{0}*ln(T) + 0
      ValueType cpIndefIntegralForEntropyAtTs   = cpMixPoly[sizePoly-1]/static_cast<double>(sizePoly-1.0);
      ValueType cpIndefIntegralForEntropyAtTref = cpMixPoly[sizePoly-1]/static_cast<double>(sizePoly-1.0);
      for (j=(sizePoly-2);j>=1;j--)
        {
        cpIndefIntegralForEntropyAtTs = cpIndefIntegralForEntropyAtTs*tmpTs + cpMixPoly[j]/static_cast<double>(j);
        cpIndefIntegralForEntropyAtTref = cpIndefIntegralForEntropyAtTs*TREF + cpMixPoly[j]/static_cast<double>(j);
        }
      cpIndefIntegralForEntropyAtTs = cpIndefIntegralForEntropyAtTs*tmpTs;
      cpIndefIntegralForEntropyAtTref = cpIndefIntegralForEntropyAtTref*tmpTs;
      ValueType tmpS = cpIndefIntegralForEntropyAtTs- cpIndefIntegralForEntropyAtTref + cpMixPoly[0]*log(tmpTs/TREF) - rmix*log(tmpPs/PREF);
      // --- Compute Entropy end ---

      // --- Compute Total Pressure (abs/rel) with entropy relationship between states ---
      ValueType cpIndefIntegralForEntropyAtTtrel = cpMixPoly[sizePoly-1]/static_cast<double>(sizePoly-1.0);
      ValueType cpIndefIntegralForEntropyAtTtabs = cpMixPoly[sizePoly-1]/static_cast<double>(sizePoly-1.0);
      for (j=(sizePoly-2);j>=1;j--)
        {
        cpIndefIntegralForEntropyAtTtabs = cpIndefIntegralForEntropyAtTtabs*tmpTtabs + cpMixPoly[j]/static_cast<double>(j);
        cpIndefIntegralForEntropyAtTtrel = cpIndefIntegralForEntropyAtTtrel*tmpTtrel + cpMixPoly[j]/static_cast<double>(j);
        }
      cpIndefIntegralForEntropyAtTtabs = cpIndefIntegralForEntropyAtTtabs*tmpTtabs;
      cpIndefIntegralForEntropyAtTtrel = cpIndefIntegralForEntropyAtTtrel*tmpTtrel;
      ValueType tmpPtrel = tmpPs * exp((cpIndefIntegralForEntropyAtTtrel - cpIndefIntegralForEntropyAtTs + cpMixPoly[0]*log(tmpTtrel/tmpTs))/rmix);
      ValueType tmpPtabs = tmpPs * exp((cpIndefIntegralForEntropyAtTtabs - cpIndefIntegralForEntropyAtTs + cpMixPoly[0]*log(tmpTtabs/tmpTs))/rmix);
      // --- Compute Total pressures end ---

      // Mandatory :
      // ts is already set
      // and cp of mixing is then set with current computed value
      // thus gamma and cv can be recomputed later knowing that
      // rmix should be stored at FielData level.
      cp[ii] = tmpCp;

      // Optionals
      if(computePs) { ps[ii] = tmpPs; }
      if(computePtabs) { ptabs[ii] = tmpPtabs; }
      if(computeTtabs) { ttabs[ii] = tmpTtabs; }
      if(computePtrel) { ptrel[ii] = tmpPtrel; }
      if(computeTtrel) { ttrel[ii] = tmpTtrel; }
      if(computeS) { s[ii] = tmpS; }
      if(computeV) { for(int dim = 0; dim < 3; dim++) { v[ii*3+dim] = tmpV[dim]; } }
      if(computeW) { for(int dim = 0; dim < 3; dim++) { w[ii*3+dim] = tmpW[dim]; } }
      if(computeM) { M[ii] = tmpMabs; }
      if(computeMrel) { Mrel[ii] = tmpMrel; }
      if(computePhi) { phi[ii] = tmpPhi; }
      if(computeAlpha) { alpha[ii] = tmpAlpha; }
      if(computeBeta) { beta[ii] = tmpBeta; }
      if(computeMu) { mu[ii] = tmpMu; }
      if(computeHtabs) { htabs[ii] = tmpHtabs; }
      if(computeHtrel) { htrel[ii] = tmpHtrel; }
      if(computeGamma) { gamma[ii] = tmpGamma; }
      }
    if (newtonNotConverged > 0)
      {
      ws.newtonNotConvergedTt += newtonNotConverged;
      }
    // if allocated do deletion :
    delete [] cpIndefIntegralPoly;
  }
  
  void Reduce()
  {
    ThreadLocalWorkSpace &ws = this->tlws.Local();
      
    // Handle Error by threads
    if (ws.newtonNotConvergedTs > 0)
      {
      std::cerr << "Newton iteration for temperature determination did not converge for "
                    << ws.newtonNotConvergedTs << " values\n";
      }
    if (ws.temperatureOverBound > 0)
      {
      std::cerr << "Temperature determination lead to "
                    << ws.temperatureOverBound << " values over bound\n";
      }
    if (ws.temperatureUnderBound > 0)
      {
      std::cerr << "Temperature determination lead to "
                    << ws.temperatureUnderBound << " values under bound\n";
      }
    if (ws.temperaturePressureEffect > 0)
      {
      std::cerr << "Temperature determination lead to "
                      << ws.temperaturePressureEffect
                      << " values with pressure effect not taken into account\n";
      }
    if (ws.newtonNotConvergedTt > 0)
      {
      std::cerr << "Newton iteration for Total temperatures determination did not converge for "
                    << ws.newtonNotConvergedTt << " values\n" ;
      }
  }
  
  private:
    struct ThreadLocalWorkSpace
    {
        int newtonNotConvergedTs;
        int temperatureOverBound;
        int temperatureUnderBound;
        int temperaturePressureEffect;      
        int newtonNotConvergedTt;
    };

    typedef vtkSMPThreadLocal<ThreadLocalWorkSpace> TLS_t;
    TLS_t tlws;
  
};

class VTK_EXPORT vtkIzarCpTReconstructFields : public vtkPointSetAlgorithm
{
public:
  static vtkIzarCpTReconstructFields* New();
  vtkTypeMacro(vtkIzarCpTReconstructFields, vtkPointSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);


  vtkGetMacro(ReplaceVars, int)
  vtkSetMacro(ReplaceVars, int)
  vtkGetMacro(AdvancedMode, int)
  vtkSetMacro(AdvancedMode, int)

  vtkGetMacro(SuthConst, double)
  vtkSetMacro(SuthConst, double)
  vtkGetMacro(SuthMuRef, double)
  vtkSetMacro(SuthMuRef, double)
  vtkGetMacro(SuthTRef, double)
  vtkSetMacro(SuthTRef, double)

  vtkSetMacro(ComputePs, int)
  vtkSetMacro(ComputePtabs, int)
  vtkSetMacro(ComputeTtabs, int)
  vtkSetMacro(ComputePtrel, int)
  vtkSetMacro(ComputeTtrel, int)
  vtkSetMacro(ComputeS, int)
  vtkSetMacro(ComputeV, int)
  vtkSetMacro(ComputeW, int)
  vtkSetMacro(ComputeM, int)
  vtkSetMacro(ComputeMrel, int)
  vtkSetMacro(ComputePhi, int)
  vtkSetMacro(ComputeAlpha, int)
  vtkSetMacro(ComputeBeta, int)
  vtkSetMacro(ComputeMu, int)
  vtkSetMacro(ComputeHtabs, int)
  vtkSetMacro(ComputeHtrel, int)
  vtkSetMacro(ComputeGamma, int)

  protected:
    vtkIzarCpTReconstructFields();
  ~vtkIzarCpTReconstructFields();

  virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);

  void ReplaceVarByArrayIfNeeded(vtkDataSet* data, vtkDataArray* ar, const char* name);

  /**
     * If true, if existing variables already exist, they are replaced
     * by the new ones. True by default
     */
  int ReplaceVars;

  /**
     * Does nothing, just used to hide the "advanced" options
     */
  int AdvancedMode;

  /**
     * Variables available in the "advanced" mode
     */
  double SuthConst;
  double SuthMuRef;
  double SuthTRef;

  int ComputePs;
  int ComputePtabs;
  int ComputeTtabs;
  int ComputePtrel;
  int ComputeTtrel;
  int ComputeS;
  int ComputeV;
  int ComputeW;
  int ComputeM;
  int ComputeMrel;
  int ComputePhi;
  int ComputeAlpha;
  int ComputeBeta;
  int ComputeMu;
  int ComputeHtabs;
  int ComputeHtrel;
  int ComputeGamma;


  /** Worker template: generates the vtkSMPTools functor and calls it
     */
  template<typename CoordsValueType, typename ValueType>
  void CallSMPCpTReconstructOp(
      vtkDataArray* coords, vtkDataArray* rho, vtkDataArray* rhov, vtkDataArray* rhoe,

      vtkDataArray* ps, vtkDataArray* ts, vtkDataArray* ptabs, vtkDataArray* ttabs,
      vtkDataArray* ptrel, vtkDataArray* ttrel, vtkDataArray* s, vtkDataArray* v, vtkDataArray* w,
      vtkDataArray* M, vtkDataArray* Mrel, vtkDataArray* phi, vtkDataArray* alpha, vtkDataArray* beta, vtkDataArray* mu,

      vtkDataArray* htabs, vtkDataArray* htrel, vtkDataArray* cp, vtkDataArray* gamma,

      int computePs, int computePtabs, int computeTtabs,
      int computePtrel, int computeTtrel, int computeS, int computeV, int computeW,
      int computeM, int computeMrel, int computePhi, int computeAlpha, int computeBeta, int computeMu,
      int computeGamma, int computeHtrel, int computeHtabs,

      double* cpMixPoly, int sizePoly,
      double rmix, double omega, int nPts)
  {
    vtkSMPCpTReconstructFieldsOp<CoordsValueType, ValueType> worker;
    worker.coords = reinterpret_cast<CoordsValueType*>(coords->GetVoidPointer(0));
    worker.rho = reinterpret_cast<ValueType*>(rho->GetVoidPointer(0));
    worker.rhov = reinterpret_cast<ValueType*>(rhov->GetVoidPointer(0));
    worker.rhoe = reinterpret_cast<ValueType*>(rhoe->GetVoidPointer(0));
    worker.ts = reinterpret_cast<ValueType*>(ts->GetVoidPointer(0));
    worker.cp = reinterpret_cast<ValueType*>(cp->GetVoidPointer(0));

    if(computePs)
      {
      worker.ps = reinterpret_cast<ValueType*>(ps->GetVoidPointer(0));
      }
    else
      {
      worker.ps = NULL;
      }

    if(computePtabs)
      {
      worker.ptabs = reinterpret_cast<ValueType*>(ptabs->GetVoidPointer(0));
      }
    else
      {
      worker.ptabs = NULL;
      }
    if(computeTtabs)
      {
      worker.ttabs = reinterpret_cast<ValueType*>(ttabs->GetVoidPointer(0));
      }
    else
      {
      worker.ttabs = NULL;
      }
    if(computePtrel)
      {
      worker.ptrel = reinterpret_cast<ValueType*>(ptrel->GetVoidPointer(0));
      }
    else
      {
      worker.ptrel = NULL;
      }
    if(computeTtrel)
      {
      worker.ttrel = reinterpret_cast<ValueType*>(ttrel->GetVoidPointer(0));
      }
    else
      {
      worker.ttrel = NULL;
      }
    if(computeS)
      {
      worker.s = reinterpret_cast<ValueType*>(s->GetVoidPointer(0));
      }
    else
      {
      worker.s = NULL;
      }
    if(computeV)
      {
      worker.v = reinterpret_cast<ValueType*>(v->GetVoidPointer(0));
      }
    else
      {
      worker.v = NULL;
      }
    if(computeW)
      {
      worker.w = reinterpret_cast<ValueType*>(w->GetVoidPointer(0));
      }
    else
      {
      worker.w = NULL;
      }
    if(computeM)
      {
      worker.M = reinterpret_cast<ValueType*>(M->GetVoidPointer(0));
      }
    else
      {
      worker.M = NULL;
      }
    if(computeMrel)
      {
      worker.Mrel = reinterpret_cast<ValueType*>(Mrel->GetVoidPointer(0));
      }
    else
      {
      worker.Mrel = NULL;
      }
    if(computePhi)
      {
      worker.phi = reinterpret_cast<ValueType*>(phi->GetVoidPointer(0));
      }
    else
      {
      worker.phi = NULL;
      }
    if(computeAlpha)
      {
      worker.alpha = reinterpret_cast<ValueType*>(alpha->GetVoidPointer(0));
      }
    else
      {
      worker.alpha = NULL;
      }
    if(computeBeta)
      {
      worker.beta = reinterpret_cast<ValueType*>(beta->GetVoidPointer(0));
      }
    else
      {
      worker.beta = NULL;
      }
    if(computeMu)
      {
      worker.mu = reinterpret_cast<ValueType*>(mu->GetVoidPointer(0));
      }
    else
      {
      worker.mu = NULL;
      }
    if(computeHtrel)
      {
      worker.htrel = reinterpret_cast<ValueType*>(htrel->GetVoidPointer(0));
      }
    else
      {
      worker.htrel = NULL;
      }
    if(computeHtabs)
      {
      worker.htabs = reinterpret_cast<ValueType*>(htabs->GetVoidPointer(0));
      }
    else
      {
      worker.htabs = NULL;
      }
    if(computeGamma)
      {
      worker.gamma = reinterpret_cast<ValueType*>(gamma->GetVoidPointer(0));
      }
    else
      {
      worker.gamma = NULL;
      }

    worker.suthConst = this->SuthConst;
    worker.suthMuRef = this->SuthMuRef;
    worker.suthTRef = this->SuthTRef;

    worker.computePs = computePs;
    worker.computePtabs = computePtabs;
    worker.computeTtabs = computeTtabs;
    worker.computePtrel = computePtrel;
    worker.computeTtrel = computeTtrel;
    worker.computeS = computeS;
    worker.computeV = computeV;
    worker.computeW = computeW;
    worker.computeM = computeM;
    worker.computeMrel = computeMrel;
    worker.computePhi = computePhi;
    worker.computeAlpha = computeAlpha;
    worker.computeBeta = computeBeta;
    worker.computeMu = computeMu;
    worker.computeHtabs = computeHtabs;
    worker.computeHtrel = computeHtrel;
    worker.computeGamma = computeGamma;

    worker.rmix = rmix;
    worker.omega = omega;
    worker.nPts = nPts;

    worker.cpMixPoly = cpMixPoly;
    worker.sizePoly = sizePoly;

    // TODO rework class definition to allow proper allocation internally
    worker.coeff = new ValueType[nPts];
    //
    vtkSMPTools::For(0, nPts, worker);
    std::cout << "End SMP For\n" << std::endl;
    //
    delete[] worker.coeff ;
    std::cout << "clear";
  }

private:
  vtkIzarCpTReconstructFields operator=(const vtkIzarCpTReconstructFields&);
  vtkIzarCpTReconstructFields(const vtkIzarCpTReconstructFields&);

};


#endif
