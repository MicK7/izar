#ifndef IZARHELPERS_H
#define IZARHELPERS_H

#include <stdexcept>
#include <string>

#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkDataArray.h"

#include <iostream>

class vtkMultiBlockDataSet;
class vtkStructuredGrid;
class vtkDataObject;
class vtkTransform;
class vtkPointSet;

namespace Izar
{
/**
 * Performs a shallowcopy on each block (the original behaviour of
 * vtkMultiBlockDataSet::ShallowCopy is to copy the structure but link
 * the nodes to the same data sets)
 */
void ShallowCopyBlocks(const vtkMultiBlockDataSet* orig, vtkMultiBlockDataSet* dest);

/**
 * Extracts a block in a multi block data set from its name
 */
vtkDataObject* GetBlockByName(vtkMultiBlockDataSet* in, const std::string& name);

/**
 * The ExtractGrid filter has difficulties when dealing with cell data.
 * This function does the same thing, but does work.
 */
vtkStructuredGrid* NewStructuredGridFromVOI(vtkStructuredGrid* in, int* VOI);

// Workers for the previous function
template <typename ValueType>
void ExtractStructuredSubsetPointDataWorker(ValueType* in, ValueType* out,
	int* oldSize, int* newSize, int* VOI, int nDims)
{
	for(int kk = 0; kk < newSize[2]; kk++)
	{
		for(int jj = 0; jj < newSize[1]; jj++)
		{
			for(int ii = 0; ii < newSize[0]; ii++)
			{
				for(int dim = 0; dim < nDims; dim++)
				{
					out[dim + nDims*(ii + newSize[0]*(jj + newSize[1]*kk))] =
						in[dim + nDims*(ii+VOI[0] + oldSize[0]*(jj+VOI[2] + oldSize[1]*(kk+VOI[4])))];
				}
			}
		}
	}
}
template <typename ValueType>
void ExtractStructuredSubsetCellDataWorker(ValueType* in, ValueType* out,
	int* oldSize, int* newSize, int* VOI, int nDims)
{
	for(int kk = 0; kk < std::max(1, newSize[2]-1); kk++)
	{
		// We take the value at the cell center of the cell at k = kk+VOI[4].
		// But we need to handle the case when this is a plane at k = cst and when k = kmax
		// In this case we need to take the value at cell = VOI[4]-1
		// Is it efficient ? 
		int oldK = kk + VOI[4];
		if((newSize[2] == 1) && (VOI[4] == oldSize[2]-1))
		{
			oldK = VOI[4]-1;
		}
		for(int jj = 0; jj < std::max(1, newSize[1]-1); jj++)
		{
			int oldJ = jj + VOI[2];
			if((newSize[1] == 1) && (VOI[2] == oldSize[1]-1))
			{
				oldJ = VOI[2]-1;
			}
			for(int ii = 0; ii < std::max(1, newSize[0]-1); ii++)
			{
				int oldI = ii + VOI[0];
				if((newSize[0] == 1) && (VOI[0] == oldSize[0]-1))
				{
					oldI = VOI[0]-1;
				}
				for(int dim = 0; dim < nDims; dim++)
				{
					out[dim + nDims*(ii + std::max(1, newSize[0]-1)*(jj + std::max(1, newSize[1]-1)*kk))] =
						in[dim + nDims*(oldI + (oldSize[0]-1)*(oldJ + (oldSize[1]-1)*oldK))];
				}
			}
		}
	}
}

vtkStructuredGrid* GetStructuredBlockFromPath(vtkMultiBlockDataSet* in, const std::string& path);
vtkDataObject* GetBlockFromPath(vtkMultiBlockDataSet* in, const std::string& path);
vtkDataObject* GetBlockFromPathWorker(vtkMultiBlockDataSet* in, const std::string& path, const std::string& fullPath);

/**
 * Fetch a fieldData of a given type
 * Returns true if it is found
 */
template <typename dataType, int dataTypeVTK>
bool GetFieldData(vtkDataSet* data, const char* name, dataType& res)
{
	if(!data) { return false; }
	vtkDataArray* ar = data->GetFieldData()->GetArray(name);
	if(!ar) { return false; }
	if(ar->GetDataType() != dataTypeVTK) { return false; }
	if(ar->GetNumberOfComponents() != 1) { return false; }
	if(ar->GetNumberOfTuples() != 1) { return false; }
	res = reinterpret_cast<dataType*>(ar->GetVoidPointer(0))[0];
	return true;
}

bool GetIntFieldData(vtkDataSet* data, const char* name, int& res);
bool GetDoubleFieldData(vtkDataSet* data, const char* name, double& res);

/**
 * Fetch a fieldData array of a given type
 * Returns true if it is found
 */
template <typename dataType, int dataTypeVTK>
bool GetFieldDataArray(vtkDataSet* data, const char* name, dataType** res, int& size)
{
	size = 0;
	if(!data) { return false; }
	vtkDataArray* ar = data->GetFieldData()->GetArray(name);
	if(!ar) { return false; }
	if(ar->GetDataType() != dataTypeVTK) { return false; }
	if(ar->GetNumberOfComponents() != 1) { return false; }
	size = ar->GetNumberOfTuples();
	*res = reinterpret_cast<dataType*>(ar->GetVoidPointer(0));
	return true;
}

bool GetDoubleFieldDataArray(vtkDataSet* data, const char* name, double** res, int& size);

// To rotate data sets
/**
 * This class abstracts the calculation of a transformation from a
 * given data set
 */
class TransformGenerator
{
public:
	virtual vtkTransform* operator() (vtkDataSet*) = 0;
};
/**
 * This class abstracts the calculation of a transformation from a
 * given data set
 */
class ConstantTransformGenerator : public TransformGenerator
{
public:
	vtkTransform* trans;
	inline virtual vtkTransform* operator() (vtkDataSet*) { return this->trans; }
};

/**
 * Generates a new object, of the same type than orig, rotated with
 * a given angle.
 */
vtkDataObject* NewRotatedDataObject(vtkDataObject* orig, double angle);

/**
 * Generates a new object, of the same type than orig, transformed with
 * transGen
 */
vtkDataObject* NewTransformedDataObject(vtkDataObject* orig, TransformGenerator& transGen);

/**
 * Generates a new multiblock data set, transformed from orig.
 * The transform is given by a functor, which returns a vtkTransform*
 * from a given vtkDataSet.
 */
vtkMultiBlockDataSet* NewTransformedMultiBlockDataSet(vtkMultiBlockDataSet* orig, TransformGenerator& transGen);
vtkPointSet* NewTransformedPointSet(vtkPointSet* orig, TransformGenerator& transGen);
void TransformDataSetAtrributes(vtkDataSetAttributes* data, vtkTransform* transform);


inline int PositiveModulo(int a, int b) { int res = a % b; if(a < 0) { return res + b; } else { return res; } }

// Exceptions

class Exception : public std::runtime_error
{
public:
	Exception(const std::string& arg): std::runtime_error(arg) {}
};

class TreeTraversalException : public Exception
{
public:
	TreeTraversalException(const std::string& arg): Exception(build_arg(arg)) {}
private:
	static const std::string build_arg(const std::string& arg)
	{
		std::string output = "Izar: error while traversing the tree:\n";
		output += arg;
		return output;
	}
};

class BlockNotFound : public Exception
{
public:
	BlockNotFound(const std::string& arg): Exception(arg) {}
};

class BadFormat : public Exception
{
public:
	BadFormat(const std::string& arg): Exception(arg) {}
};

}
#endif
