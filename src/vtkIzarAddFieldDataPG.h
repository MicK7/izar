#ifndef VTKIZARADDFIELDDATAPG_H
#define VTKIZARADDFIELDDATAPG_H

#include "IzarDefines.h"

#include "vtkIzarAddFieldData.h"

/**
 * Interface of vtkIzarAddFieldData to be used as a PV plugin
 * 
 * Allows the user to add specific FieldData, related with a 
 * turbomachinery computation with perfect gas :
 *   - Gamma (heat capacity ratio)
 *   - Rgas (perfect gas constant)
 *   - Omega (rotation speed of the dataset)
 *   - Zsector (number of duplications needed to cover the full 360° domain)
 */
class VTK_EXPORT vtkIzarAddFieldDataPG : public vtkIzarAddFieldData
{
public:
    static vtkIzarAddFieldDataPG* New();
    vtkTypeMacro(vtkIzarAddFieldDataPG, vtkIzarAddFieldData);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    void SetGamma(double val) { this->doubleDataValues[0] = val; this->Modified(); };
    double GetGamma() { return this->doubleDataValues[0]; };
    void SetRgas(double val) { this->doubleDataValues[1] = val; this->Modified(); };
    double GetRgas() { return this->doubleDataValues[1]; };
    void SetOmega(double val) { this->doubleDataValues[2] = val; this->Modified(); };
    double GetOmega() { return this->doubleDataValues[2]; };
    void SetZsector(int val) { this->intDataValues[0] = val; this->Modified(); };
    int GetZsector() { return this->intDataValues[0]; };
protected:
	vtkIzarAddFieldDataPG();
	virtual ~vtkIzarAddFieldDataPG();
	
private:
    vtkIzarAddFieldDataPG operator=(const vtkIzarAddFieldDataPG&);
    vtkIzarAddFieldDataPG(const vtkIzarAddFieldDataPG&);
};

#endif
