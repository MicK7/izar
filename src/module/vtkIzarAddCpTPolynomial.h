#ifndef VTKIZARADDCPTPOLYNOMIAL_H
#define VTKIZARADDCPTPOLYNOMIAL_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkDataSetAlgorithm.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;

/**
 *  -- WIP --
 * Allow the user to add specific Cp polynomial for
 *   Air, Fuel and Water.
 */
class VTK_EXPORT vtkIzarAddCpTPolynomial : public vtkDataSetAlgorithm
{
public:
    static vtkIzarAddCpTPolynomial* New();
    vtkTypeMacro(vtkIzarAddCpTPolynomial, vtkDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    void SetCpAirPoly(double coeffs[8])
    {
      for(int ii=0; ii<8; ii++)
      {
        this->polynomialCoeffs[0][ii] = coeffs[ii];
      }
      this->Modified();
    };

    void SetCpFuelPoly(double coeffs[8])
    {
      for(int ii=0; ii<8; ii++)
      {
        this->polynomialCoeffs[1][ii] = coeffs[ii];
      }
      this->Modified();
    };

    void SetCpWaterPoly(double coeffs[8])
    {
      for(int ii=0; ii<8; ii++)
      {
        this->polynomialCoeffs[2][ii] = coeffs[ii];
      }
      this->Modified();
    };
    
    /**
      * Add a double to the data array to add in the FieldData
      */
    void AddDoubleDataArray(const std::string& name, double* data, int size);

    /**
     * Clear the data names and values previously added
     */
    void clearData();
    
protected:
  vtkIzarAddCpTPolynomial();
  virtual ~vtkIzarAddCpTPolynomial();
  
  virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
  
  
  /**
   * Names of the polynomials to add in the FieldData
   */
  std::vector<std::string> polynomialDataNames;
  /**
   * Polynomial coefficient to add in the FieldData
   */
  std::vector<std::vector<double>> polynomialCoeffs;

private:
    vtkIzarAddCpTPolynomial operator=(const vtkIzarAddCpTPolynomial&);
    vtkIzarAddCpTPolynomial(const vtkIzarAddCpTPolynomial&);
};

#endif
