#ifndef VTKIZARADDFIELDDATAGENERIC_H
#define VTKIZARADDFIELDDATAGENERIC_H

#include "IzarDefines.h"

#include "vtkIzarAddFieldData.h"

#include <string>

/**
 * Interface of vtkIzarAddFieldData to be used as a PV plugin
 * 
 * Allows the user to add generic FieldData. The interface is a text
 * field. The text field can contain one line per FieldData to add.
 * Each line contains :
 *   - The FieldData name. If it contains spaces, it must be enclosed in
 *     quotes.
 *   - One or more spaces
 *   - The value. If it does not contain a decimal separator ("."), it
 *     is interpreted as an int. Otherwise, it is interpreted as a
 *     double.
 *  Blank lines are allowed. Tab characters are interpreted as spaces.
 *  If the formatting of one or more lines is uncorrect, an error
 *  message is printed and no FieldData is added.
 */
class VTK_EXPORT vtkIzarAddFieldDataGeneric : public vtkIzarAddFieldData
{
public:
    static vtkIzarAddFieldDataGeneric* New();
    vtkTypeMacro(vtkIzarAddFieldDataGeneric, vtkIzarAddFieldData);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    /**
     * Standard getter
     */
    virtual const char* GetElementList() { return this->ElementList.c_str(); }
    /**
     * Setter, contains the FieldData list parsing logic.
     */
    virtual void SetElementList(const char* str);
    
protected:
	vtkIzarAddFieldDataGeneric();
	virtual ~vtkIzarAddFieldDataGeneric();
	
	std::string ElementList;
	
private:
    vtkIzarAddFieldDataGeneric operator=(const vtkIzarAddFieldDataGeneric&);
    vtkIzarAddFieldDataGeneric(const vtkIzarAddFieldDataGeneric&);
};

#endif
