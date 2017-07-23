#ifndef VTKIZARXRCUTFROMCSV_H
#define VTKIZARXRCUTFROMCSV_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkIzarGenericXRCut.h"

#include <string>

class VTK_EXPORT vtkIzarXRCutFromCSV : public vtkIzarGenericXRCut
{
public:
    static vtkIzarXRCutFromCSV* New();
    vtkTypeMacro(vtkIzarXRCutFromCSV, vtkIzarGenericXRCut);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    /**
     * Standard Setter
     */
    void SetFileName(const char* name) { this->FileName = name; this->Modified(); }
    /**
     * Standard Getter
     */
    const char* GetFileName() { return this->FileName.c_str(); }
    
protected:
    vtkIzarXRCutFromCSV();
    ~vtkIzarXRCutFromCSV();
    
    virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
    
    std::string FileName;
    
    /**
     * Parses a CSV file containing 2 columns : X and R
     */
    bool FillXRFromCSVFile();
    
private:
    vtkIzarXRCutFromCSV operator=(const vtkIzarXRCutFromCSV&);
    vtkIzarXRCutFromCSV(const vtkIzarXRCutFromCSV&);
};

#endif
