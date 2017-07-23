
#include "vtkIzarAddFieldDataGeneric.h"

#include "vtkObjectFactory.h"

#include <sstream>
#include <stdexcept>

vtkStandardNewMacro(vtkIzarAddFieldDataGeneric)

vtkIzarAddFieldDataGeneric::vtkIzarAddFieldDataGeneric() : vtkIzarAddFieldData()
{
	IZAR_WARNING
}

vtkIzarAddFieldDataGeneric::~vtkIzarAddFieldDataGeneric()
{
}

void vtkIzarAddFieldDataGeneric::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarAddFieldDataGeneric";
	this->Superclass::PrintSelf(os, indent.GetNextIndent());
}

void vtkIzarAddFieldDataGeneric::SetElementList(const char* str)
{
	// We don't use the operator= directly because we want to manage the
	// replacement of the \t by spaces when copying str to ElementSize
	int n = strlen(str);
	this->ElementList.resize(n);
	for(int ii = 0; ii < n; ii++)
	{
		if(str[ii] == '\t')
		{
			this->ElementList[ii] = ' ';
		}
		else
		{
			this->ElementList[ii] = str[ii];
		}
	}
	this->ElementList = str;
	
	this->clearData();
	int iGlobal = 0;
	int iLine = 0;
	while (iGlobal < n)
	{
		// Extract current line
		int lineBegin = iGlobal;
		int lineEnd = iGlobal;
		while((lineEnd < n) && (this->ElementList[lineEnd] != '\n'))
		{
			lineEnd++;
		}
		std::string curLine = this->ElementList.substr(lineBegin, lineEnd-lineBegin);
		vtkDebugMacro("Read line " << curLine)
		// Test if line is blank
		bool isBlank = true;
		for(int ii = 0 ; ii < curLine.size() ; ii++)
		{
			if(curLine[ii] != ' ')
			{
				isBlank = false;
				break;
			}
		}
		// If it is not blank: process it
		if(!isBlank)
		{
			// First, get the field name
			int beginFieldName = curLine.find_first_not_of(' '); // Should not return npos because we tested that the line is not blank
			int endFieldName;
			if(curLine[beginFieldName] == '"')
			{
				// Quotes : take into account spaces
				beginFieldName++;
				endFieldName = curLine.find('"', beginFieldName);
				if(endFieldName == std::string::npos)
				{
					vtkErrorMacro("Error while parsing line " << iLine << ":\n" <<
						"\"" << curLine << "\"\n" << 
						"FieldData name started with a quote (\"), but does not contains the ending quote")
					this->clearData();
					return;
				}
			}
			else
			{
				// No quotes
				endFieldName = curLine.find(' ', beginFieldName);
				if(endFieldName == std::string::npos)
				{
					vtkErrorMacro("Error while parsing line " << iLine << ":\n" <<
						"\"" << curLine << "\"\n" <<
						"Cannot find the value")
					this->clearData();
					return;
				}
			}
			std::string fieldName = curLine.substr(beginFieldName, endFieldName-beginFieldName);
			vtkDebugMacro("FieldName read : " << fieldName)
			
			// Then, the value
			int beginValue = curLine.find_first_not_of(' ', endFieldName+1);
			if(beginValue == std::string::npos)
			{
				vtkErrorMacro("Error while parsing line " << iLine << ":\n" <<
					"\"" << curLine << "\"\n" <<
					"Cannot find the value")
				this->clearData();
				return;
			}
			int endValue = curLine.find(' ', beginValue);
			std::string valueStr;
			if(endValue != std::string::npos) // There are trailing spaces: check that there are no characters different from spaces after the value
			{
				int notBlankedChar = curLine.find_first_not_of(' ', endValue);
				if(notBlankedChar != std::string::npos)
				{
					vtkErrorMacro("Error while parsing line " << iLine << ":\n" <<
						"\"" << curLine << "\"\n" <<
						"Trailing characters after the value")
					this->clearData();
					return;
				}
				valueStr = curLine.substr(beginValue);
			}
			else
			{
				valueStr = curLine.substr(beginValue, endValue-beginValue);
			}
			// Convert the value to int or double
			if(valueStr.find('.') != std::string::npos)
			{
				// This is a double
				double val;
				try
				{
					size_t end;
					val = std::stod(valueStr, &end);
					if(end != valueStr.size())
					{
						throw std::invalid_argument("Invalid double");
					}
					
				}
				catch(const std::invalid_argument& ia)
				{
					vtkErrorMacro("Error while parsing line " << iLine << ":\n" <<
						"\"" << curLine << "\"\n" <<
						"Cannot convert \"" << valueStr << "\" to double\n" <<
						"(" << ia.what() << ")")
					this->clearData();
					return;
				}
				vtkDebugMacro("Double value read " << val)
				this->AddDoubleData(fieldName, val);
			}
			else
			{
				// This is an int
				int val;
				try
				{
					size_t end;
					val = std::stoi(valueStr, &end);
					if(end != valueStr.size())
					{
						throw std::invalid_argument("Invalid int");
					}
				}
				catch(const std::invalid_argument& ia)
				{
					vtkErrorMacro("Error while parsing line " << iLine << ":\n" <<
						"\"" << curLine << "\"\n" <<
						"Cannot convert \"" << valueStr << "\" to int\n" <<
						"(" << ia.what() << ")")
					this->clearData();
					return;
				}
				vtkDebugMacro("Int value read " << val)
				this->AddIntData(fieldName, val);
			}
		}
		iLine++;
		iGlobal = lineEnd+1;
	}
	this->Modified();
}
