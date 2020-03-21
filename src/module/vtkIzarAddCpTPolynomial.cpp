
#include "vtkIzarAddCpTPolynomial.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDoubleArray.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"


vtkStandardNewMacro(vtkIzarAddCpTPolynomial)

vtkIzarAddCpTPolynomial::vtkIzarAddCpTPolynomial()
{
  this->polynomialDataNames.resize(3);
  this->polynomialCoeffs.resize(3);

  this->polynomialDataNames[0] = CPAIR;
  this->polynomialDataNames[1] = CPFUEL;
  this->polynomialDataNames[2] = CPWATER;

  this->polynomialCoeffs[0].resize(8);
  this->polynomialCoeffs[1].resize(8);
  this->polynomialCoeffs[2].resize(8);

  IZAR_WARNING
}

vtkIzarAddCpTPolynomial::~vtkIzarAddCpTPolynomial()
{
}

void vtkIzarAddCpTPolynomial::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "vtkIzarAddCpTPolynomial\n";
  os << indent << "List of FieldData Polynomials to add\n";
  vtkIndent nextIndent = indent.GetNextIndent();
  os << indent << this->polynomialDataNames.size()
               << " polynomial coefficients:\n";
  for(int ii = 0; ii < this->polynomialDataNames[ii].size(); ii++)
  {
    os << nextIndent << this->polynomialDataNames[ii] << " -> " << "\n";
    for(int jj =0; jj < this->polynomialCoeffs[ii].size(); jj++)
    {
    os << nextIndent << "    A[" << jj <<"] = "
                     << this->polynomialCoeffs[ii][jj] << "\n";
 
    }
  }
}

void vtkIzarAddCpTPolynomial::AddDoubleDataArray(const std::string& name, double* data, int size)
{
  this->polynomialDataNames.push_back(name);
  std::vector<double> tmp;
  tmp.assign(data, data+size);
  this->polynomialCoeffs.push_back(tmp);
}

void vtkIzarAddCpTPolynomial::clearData()
{
  this->polynomialDataNames.clear();
  this->polynomialCoeffs.clear();
}

int vtkIzarAddCpTPolynomial::RequestData(vtkInformation* request,
                                          vtkInformationVector** inVector,
                                          vtkInformationVector* outVector)
{
  vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
  if(!inInfo)
  {
    return 0;
  }
  vtkDataSet* inData = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if(!inData)
  {
    return 0;
  }
  
  vtkInformation* outInfo = outVector->GetInformationObject(0);
  if(!outInfo)
  {
    return 0;
  }
  vtkDataSet* outData = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if(!outData)
  {
    return 0;
  }
  outData->ShallowCopy(inData);

  for(int ii = 0; ii < this->polynomialDataNames.size(); ii++)
  {
    vtkDoubleArray* ar = vtkDoubleArray::New();
    ar->SetNumberOfComponents(1);
    ar->SetNumberOfTuples(this->polynomialCoeffs[ii].size());
    double* data = reinterpret_cast<double*>(ar->GetVoidPointer(0));
    for(int jj =0; jj < this->polynomialCoeffs[ii].size(); jj++)
    {
       data[jj] = this->polynomialCoeffs[ii][jj];
    }
    ar->SetName(this->polynomialDataNames[ii].c_str());
    outData->GetFieldData()->AddArray(ar);
    ar->Delete();
  }
  return 1;
}
