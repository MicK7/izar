#include "vtkIzarXRCutFromCSV.h"

#include "vtkObjectFactory.h"

#include <fstream>
#include <stdexcept>

vtkStandardNewMacro(vtkIzarXRCutFromCSV)

vtkIzarXRCutFromCSV::vtkIzarXRCutFromCSV(): vtkIzarGenericXRCut(),
	FileName("")
{
}

vtkIzarXRCutFromCSV::~vtkIzarXRCutFromCSV()
{
}

void vtkIzarXRCutFromCSV::PrintSelf(std::ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarXRCutFromCSV\n";
}

int vtkIzarXRCutFromCSV::RequestData(vtkInformation* request,
	vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	if(!this->FillXRFromCSVFile())
	{
		return 0;
	}
	return this->Superclass::RequestData(request, inVector, outVector);
}

bool vtkIzarXRCutFromCSV::FillXRFromCSVFile()
{
	// Not efficient but this function is not a bottleneck as we
	// only read small files
	
	this->XList.clear();
	this->RList.clear();
	
	std::ifstream ifs;
	ifs.open(this->FileName);
	if(ifs.fail())
	{
		vtkErrorMacro("Could not open file " << this->FileName);
		return false;
	}
	std::string number_chars("0123456789.");
	ifs.seekg(0);
	int iLine = 0;
	while(!ifs.eof())
	{
		double x, r;
		char dummy;
		ifs >> x;
		if(ifs.eof())
		{
			break;
		}
		if(ifs.fail())
		{
			vtkErrorMacro("Line " << iLine << ": unable to read the x value");
			return false;
		}
		ifs.get(dummy); // Separator
		if(ifs.fail() || (number_chars.find(dummy) != std::string::npos)) // If the separator is not really a separator
		{
			vtkErrorMacro("Line " << iLine << ": could not read after " << x);
			return false;
		}
		ifs >> r;
		if(ifs.fail())
		{
			vtkErrorMacro("Line " << iLine << ": unable to read the r value");
			return false;
		}
		this->XList.push_back(x);
		this->RList.push_back(r);
		
		// Go to the next line (or the end of file)
		dummy = '\0';
		while((!ifs.eof()) && (dummy != '\n'))
		{
			ifs.get(dummy);
		}
		iLine++;
	}
	return true;
}
	/*
	std::string content;
	ifs >> content;
	std::cout << content << std::endl;
	ifs.close();
	
	// Detects the delimiter by reading the first line
	int posDelim = content.find_first_not_of("0123456789.");
	if(posDelim == std::string::npos)
	{
		vtkErrorMacro("In file " << this->FileName << ": not enough columns");
		return false;
	}
	char delim = content[posDelim];
	
	// Detects whether we are using \r\n end of lines or just \n
	bool eolWin = true;
	int posEOL = content.find('\n');
	if(posEOL == std::string::npos)
	{
		vtkErrorMacro("In file " << this->FileName << ": not enough lines");
		return false;
	}
	if((posEOL == 0) || (content[posEOL-1] != '\r'))
	{
		eolWin = false;
	} 
	char eolDetector = '\n';
	if(eolWin)
	{
		eolDetector = '\r';
	}
	
	// Read the data
	int pos = 0;
	int iLine = 0;
	while(pos < content.size())
	{
		// Read a line
		// First : look for the delimiter
		posDelim = content.find(delim, pos);
		if(posDelim == std::string::npos)
		{
			vtkErrorMacro("In file " << this->FileName << ": line " << iLine << ": cannot find the delimiter " << delim);
			return false;
		}
		// We've read the x value
		double x;
		try
		{
			x = std::stod(content.substr(pos, posDelim-pos));
		}
		catch(const std::invalid_argument& e)
		{
			vtkErrorMacro("In file " << this->FileName << ": could not read value '" << content.substr(pos, posDelim-pos) << "'");
			return false;
		}
		// Look for the eol character
		pos = posDelim+1;
		posEOL = content.find(eolDetector, pos);
		if(posEOL == std::string::npos)
		{
			// This is the last element of the file
			posEOL = content.size();
		}
		double r;
		try
		{
			r = std::stod(content.substr(pos, posEOL-pos));
		}
		catch(const std::invalid_argument& e)
		{
			vtkErrorMacro("In file " << this->FileName << ": could not read value '" << content.substr(pos, posEOL-pos) << "'");
			return false;
		}
		if(eolWin)
		{
			pos += 2;
		}
		else
		{
			pos += 1;
		}
		iLine += 1;
		this->XList.push_back(x);
		this->RList.push_back(r);
	}
	* */

