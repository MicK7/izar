#ifndef VTKIZARADDFIELDDATACPT_H
#define VTKIZARADDFIELDDATACPT_H

#include "IzarDefines.h"

#include "vtkIzarAddFieldData.h"

/**
 * Interface of vtkIzarAddFieldData to be used as a PV plugin
 * 
 * Allows the user to add specific FieldData, related with a 
 * turbomachinery computation with perfect gas :
 *   - Rgas (perfect gas constant for air)
 *   - Rfuel (perfect gas constant for fuel)
 *   - Rwater (perfect gas constant for water)
 *   - FAR
 *   - WAR
 *   - Omega (rotation speed of the dataset)
 *   - Zsector (number of duplications needed to cover the full 360° domain)
 */
class VTK_EXPORT vtkIzarAddFieldDataCpT : public vtkIzarAddFieldData
{
public:
    static vtkIzarAddFieldDataCpT* New();
    vtkTypeMacro(vtkIzarAddFieldDataCpT, vtkIzarAddFieldData);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    void SetRgas(double val) { this->doubleDataValues[0] = val; this->Modified(); };
    double GetRgas() { return this->doubleDataValues[0]; };

    void SetRfuel(double val) { this->doubleDataValues[1] = val; this->Modified(); };
    double GetRfuel() { return this->doubleDataValues[1]; };

    void SetRwater(double val) { this->doubleDataValues[2] = val; this->Modified(); };
    double GetRwater() { return this->doubleDataValues[2]; };

    void SetFAR(double val) { this->doubleDataValues[3] = val; this->Modified(); };
    double GetFAR() { return this->doubleDataValues[3]; };

    void SetWAR(double val) { this->doubleDataValues[4] = val; this->Modified(); };
    double GetWAR() { return this->doubleDataValues[4]; };

    void SetOmega(double val) { this->doubleDataValues[5] = val; this->Modified(); };
    double GetOmega() { return this->doubleDataValues[5]; };

    void SetZsector(int val) { this->intDataValues[0] = val; this->Modified(); };
    int GetZsector() { return this->intDataValues[0]; };

protected:
	vtkIzarAddFieldDataCpT();
	virtual ~vtkIzarAddFieldDataCpT();
	
private:
    vtkIzarAddFieldDataCpT operator=(const vtkIzarAddFieldDataCpT&);
    vtkIzarAddFieldDataCpT(const vtkIzarAddFieldDataCpT&);
};

#endif
