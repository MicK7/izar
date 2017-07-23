#ifndef VTKIZARCONECUT_H
#define VTKIZARCONECUT_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkIzarGenericXRCut.h"

class VTK_EXPORT vtkIzarConeCut : public vtkIzarGenericXRCut
{
public:
    static vtkIzarConeCut* New();
    vtkTypeMacro(vtkIzarConeCut, vtkIzarGenericXRCut);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    void SetXR1(double x, double r) { this->XList[0] = x; this->RList[0] = r; this->Modified(); }
    void SetXR2(double x, double r) { this->XList[1] = x; this->RList[1] = r; this->Modified(); }
    void SetClipExtremities(int b) { this->ClipFirstPoint = b; this->ClipLastPoint = b; this->Modified(); }
    
protected:
    vtkIzarConeCut();
    ~vtkIzarConeCut();
    
private:
    vtkIzarConeCut operator=(const vtkIzarConeCut&);
    vtkIzarConeCut(const vtkIzarConeCut&);
};

#endif
