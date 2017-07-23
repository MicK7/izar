#include "IzarHelpers.h"

#include "vtkMultiBlockDataSet.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkCompositeDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkInformation.h"
#include "vtkSetGet.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkPoints.h"
#include "vtkTransform.h"
#include "vtkPointSet.h"
#include "vtkMath.h"

#include <sstream>

namespace Izar
{
	
void ShallowCopyBlocks(vtkMultiBlockDataSet* orig, vtkMultiBlockDataSet* dest)
{
	dest->Initialize();
	dest->CopyStructure(orig);
	
	vtkDataObjectTreeIterator* iterOrig = orig->NewTreeIterator();
	vtkDataObjectTreeIterator* iterDest = dest->NewTreeIterator();
	iterOrig->SkipEmptyNodesOff();
	iterOrig->TraverseSubTreeOn();
	iterOrig->VisitOnlyLeavesOn();
	iterDest->SkipEmptyNodesOff();
	iterDest->TraverseSubTreeOn();
	iterDest->VisitOnlyLeavesOn();
	iterOrig->InitTraversal();
	iterDest->InitTraversal();
	
	while(!iterOrig->IsDoneWithTraversal())
	{
		vtkDataObject* origObject = orig->GetDataSet(iterOrig);
		if(origObject != NULL)
		{
			vtkDataObject* destObject = origObject->NewInstance();
			destObject->ShallowCopy(origObject);
			dest->SetDataSet(iterDest, destObject);
			destObject->Delete();
		}
		iterOrig->GoToNextItem();
		iterDest->GoToNextItem();
	}
	iterOrig->Delete();
	iterDest->Delete();
}

vtkDataObject* GetBlockByName(vtkMultiBlockDataSet* in, const std::string& name)
{
	int nBks = in->GetNumberOfBlocks();
	bool found = false;
	int ii = 0;
	int iBk = 0;
	while((ii < nBks) && (!found))
	{
		if(in->HasMetaData(ii))
		{
			vtkInformation* meta = in->GetMetaData(ii);
			if(meta->Has(vtkCompositeDataSet::NAME()))
			{
				std::string curName = meta->Get(vtkCompositeDataSet::NAME());
				if(name == curName)
				{
					found = true;
					iBk = ii;
				}
			}
		}
		ii++;
	}
	if(!found)
	{
		std::stringstream err;
		err << "Block " << name << " not found";
		throw Izar::BlockNotFound(err.str());
	}
	return in->GetBlock(iBk);
}

vtkStructuredGrid* NewStructuredGridFromVOI(vtkStructuredGrid* in, int* VOI)
{
	if((in == NULL) || (VOI[1] < VOI[0]) || (VOI[3] < VOI[2]) || (VOI[5] < VOI[4]))
	{
		return NULL;
	}
	int size[3];
	int newSize[3];
	int newNPts = 1;
	int newNCells = 1;
	in->GetDimensions(size);
	
	// Clip the VOI to acceptable values
	int clippedVOI[6];
	for(int ii = 0; ii < 6; ii++)
	{
		if(VOI[ii] < 0)
		{
			clippedVOI[ii] = 0;
		}
		else if(VOI[ii] >= size[ii/2])
		{
			clippedVOI[ii] = size[ii/2]-1;
		}
		else
		{
			clippedVOI[ii] = VOI[ii];
		}
	}
	for(int ii = 0; ii < 3; ii++)
	{
		newSize[ii] = clippedVOI[2*ii+1]-clippedVOI[2*ii]+1;
		newNPts *= newSize[ii];
		if(newSize[ii] != 1)
		{
			newNCells *= newSize[ii]-1;
		}
	}
	vtkStructuredGrid* res = vtkStructuredGrid::New();
	res->SetDimensions(newSize);
	
	// Generate the new points
	vtkPoints* pts = vtkPoints::New();
	vtkPoints* oldPts = in->GetPoints();
	pts->SetDataType(oldPts->GetDataType());
	pts->SetNumberOfPoints(newNPts);
	
	switch(pts->GetDataType())
	{
		vtkTemplateMacro(ExtractStructuredSubsetPointDataWorker<VTK_TT>(
			static_cast<VTK_TT*>(oldPts->GetVoidPointer(0)),
			static_cast<VTK_TT*>(pts->GetVoidPointer(0)),
			size, newSize, clippedVOI, 3));
	}
	
	res->SetPoints(pts);
	pts->Delete();
	
	// Generate the new point data
	vtkPointData* newPtData = res->GetPointData();
	vtkPointData* oldPtData = in->GetPointData();
	for(int iAr = 0; iAr < oldPtData->GetNumberOfArrays(); iAr++)
	{
		vtkDataArray* oldAr = oldPtData->GetArray(iAr);
		vtkDataArray* newAr = oldAr->NewInstance();
		int nDims = oldAr->GetNumberOfComponents();
		newAr->SetNumberOfComponents(nDims);
		newAr->SetNumberOfTuples(newNPts);
		switch(newAr->GetDataType())
		{
			vtkTemplateMacro(ExtractStructuredSubsetPointDataWorker<VTK_TT>(
				static_cast<VTK_TT*>(oldAr->GetVoidPointer(0)),
				static_cast<VTK_TT*>(newAr->GetVoidPointer(0)),
				size, newSize, clippedVOI, nDims));
		}
		newAr->SetName(oldAr->GetName());
		newPtData->AddArray(newAr);
		newAr->Delete();
	}
	
	
	// Generate the new cell data
	vtkCellData* newCData = res->GetCellData();
	vtkCellData* oldCData = in->GetCellData();
	for(int iAr = 0; iAr < oldCData->GetNumberOfArrays(); iAr++)
	{
		vtkDataArray* oldAr = oldCData->GetArray(iAr);
		vtkDataArray* newAr = oldAr->NewInstance();
		int nDims = oldAr->GetNumberOfComponents();
		newAr->SetNumberOfComponents(nDims);
		newAr->SetNumberOfTuples(newNCells);
		switch(newAr->GetDataType())
		{
			vtkTemplateMacro(ExtractStructuredSubsetCellDataWorker<VTK_TT>(
				static_cast<VTK_TT*>(oldAr->GetVoidPointer(0)),
				static_cast<VTK_TT*>(newAr->GetVoidPointer(0)),
				size, newSize, clippedVOI, nDims));
		}
		newAr->SetName(oldAr->GetName());
		newCData->AddArray(newAr);
		newAr->Delete();
	}
	
	// Shallow copy the field data
	vtkFieldData* newFData = res->GetFieldData();
	vtkFieldData* oldFData = in->GetFieldData();
	for(int iAr = 0; iAr < oldFData->GetNumberOfArrays(); iAr++)
	{
		vtkDataArray* oldAr = oldFData->GetArray(iAr);
		newFData->AddArray(oldAr);
	}
	
	return res;
}

vtkStructuredGrid* GetStructuredBlockFromPath(vtkMultiBlockDataSet* in, const std::string& path)
{
	vtkDataObject* obj = Izar::GetBlockFromPath(in, path);
	if(!obj)
	{
		std::stringstream err;
		err << "Block " << path << " is NULL";
		throw Izar::BlockNotFound(err.str());
	}
	vtkStructuredGrid* res = vtkStructuredGrid::SafeDownCast(obj);
	if(!res)
	{
		std::stringstream err;
		err << "Block " << path << " is not a structured block (it is a " << obj->GetClassName();
		throw Izar::BlockNotFound(err.str());
	}
	return res;
}

vtkDataObject* GetBlockFromPath(vtkMultiBlockDataSet* in, const std::string& path)
{
	return Izar::GetBlockFromPathWorker(in, path, path);
}

vtkDataObject* GetBlockFromPathWorker(vtkMultiBlockDataSet* in, const std::string& path, const std::string& fullPath)
{
	// Remove beginning and trailing spaces, and remove trailing /
	int iBegin = path.find_first_not_of(' ');
	int iEnd = path.find_last_not_of(" /"); // Also remove trailing '/'
	if((iBegin == std::string::npos) || (iEnd == std::string::npos))
	{
		std::stringstream err;
		err << "Empty path: '" << path << "' in '" << fullPath << "'";
		throw Izar::BlockNotFound(err.str());
	}
	std::string cleanPath = path.substr(iBegin, iEnd+1-iBegin);
	if(cleanPath[0] != '/')
	{
		std::stringstream err;
		err << "Path must start by '/' (current: '" << fullPath << "')";
		throw Izar::BlockNotFound(err.str());
	}
	int endBkName = cleanPath.find('/', 1);
	if(endBkName == std::string::npos)
	{
		endBkName = cleanPath.size();
	}
	std::string bkName = cleanPath.substr(1, endBkName-1);
	std::string nextBkName = cleanPath.substr(endBkName);
	// If nextBkName is empty, this must be a structured grid. Else, this must be a multiblockdataset
	vtkDataObject* block;
	try
	{
		block = Izar::GetBlockByName(in, bkName);
	}
	catch(const Izar::BlockNotFound& e)
	{
		std::stringstream err; // Enhance the original error message
		err << "Block " << bkName << " not found in " << fullPath;
		throw Izar::BlockNotFound(err.str());
	}
	if(nextBkName.size() == 0)
	{
		return block;
	}
	else if(! block->IsA("vtkMultiBlockDataSet"))
	{
		std::stringstream err;
		err << "Block " << bkName << " in " << fullPath << " is not a multi block data set";
		throw Izar::BlockNotFound(err.str());
	}
	else
	{
		return GetBlockFromPathWorker(vtkMultiBlockDataSet::SafeDownCast(block),
			nextBkName, fullPath);
	}
}

bool GetIntFieldData(vtkDataSet* data, const char* name, int& res) { return Izar::GetFieldData<int, VTK_INT>(data, name, res); }

bool GetDoubleFieldData(vtkDataSet* data, const char* name, double& res) { return Izar::GetFieldData<double, VTK_DOUBLE>(data, name, res); }

bool GetDoubleFieldDataArray(vtkDataSet* data, const char* name, double** res, int& size) { return Izar::GetFieldDataArray<double, VTK_DOUBLE>(data, name, res, size); }

vtkDataObject* NewRotatedDataObject(vtkDataObject* orig, double angle)
{
	// Builds the transform object
    vtkTransform* transform = vtkTransform::New();
    transform->Identity();
    transform->RotateX(angle*180./vtkMath::Pi());
    
    ConstantTransformGenerator transGen;
    transGen.trans = transform;
    
    vtkDataObject* res = NewTransformedDataObject(orig, transGen);
    
    transform->Delete();
    return res;
}

vtkDataObject* NewTransformedDataObject(vtkDataObject* orig, TransformGenerator& transGen)
{  
    vtkDataObject* res = NULL;    
    // For now, this function handles only 2 data types : point sets and multi block data sets
    vtkPointSet* origAsPointSet = vtkPointSet::SafeDownCast(orig);
    if(origAsPointSet)
    {
        // This is a vtkPointSet
        res = NewTransformedPointSet(origAsPointSet, transGen);
    }
    else
    {
        vtkMultiBlockDataSet* origMB = vtkMultiBlockDataSet::SafeDownCast(orig);
        if(origMB)
        {
            res = NewTransformedMultiBlockDataSet(origMB, transGen);
        }
        else
        {
            std::cout << "Error: trying to rotate an unsupported data type" << std::endl;
        }
    }
    return res;
}

vtkPointSet* NewTransformedPointSet(vtkPointSet* orig, TransformGenerator& transGen)
{
    vtkPointSet* res = orig->NewInstance();
    res->ShallowCopy(orig);
    
    vtkTransform* transform = transGen(orig);
    
    // Initialize a new vtkPoints object
    vtkPoints* newPts = vtkPoints::New();
    vtkPoints* origPts = orig->GetPoints();
    if(origPts)
    {
		newPts->SetDataType(origPts->GetDataType());
		
		// Apply the transformation to the points
		transform->TransformPoints(origPts, newPts);
		res->SetPoints(newPts);
		newPts->Delete();
		
		// Shallow copy and transform the point and cell data
		TransformDataSetAtrributes(res->GetPointData(), transform);
		TransformDataSetAtrributes(res->GetCellData(), transform);
	}
    
    return res;
}

vtkMultiBlockDataSet* NewTransformedMultiBlockDataSet(vtkMultiBlockDataSet* orig, TransformGenerator& transGen)
{
    vtkMultiBlockDataSet* res = vtkMultiBlockDataSet::New();
    res->CopyStructure(orig);
    
    vtkDataObjectTreeIterator* it = orig->NewTreeIterator();
    it->VisitOnlyLeavesOn();
    it->TraverseSubTreeOn();
    it->SkipEmptyNodesOff();
    
    while(!it->IsDoneWithTraversal())
    {
        vtkDataObject* curObj = orig->GetDataSet(it);
        if(curObj)
        {
            vtkPointSet* origPtSet = vtkPointSet::SafeDownCast(curObj);
            if(origPtSet)
            {
                vtkPointSet* newPtSet = NewTransformedPointSet(origPtSet, transGen);
                res->SetDataSet(it, newPtSet);
                newPtSet->Delete();
            }
            else
            {
                std::cout << "Error: trying to transform a vtkMultiBlockDataSet containing an unsupported data type" << std::endl;
            }
        }
        it->GoToNextItem();
    }
    return res;
}

void TransformDataSetAtrributes(vtkDataSetAttributes* data, vtkTransform* transform)
{   
    // To handle the normals
    vtkDataArray* normal = data->GetNormals();
    
    // Traverses the arrays
    vtkIdType nArrays = data->GetNumberOfArrays();
    
    std::vector<int> toDelete;
    toDelete.reserve(nArrays);
    std::vector<vtkDataArray*> toAdd;
    toAdd.reserve(nArrays);
    vtkDataArray* normalToAdd = NULL;
    
    for(vtkIdType ii = 0; ii < nArrays; ii++)
    {
        vtkDataArray* origAr = data->GetArray(ii);
        if(origAr->GetNumberOfComponents() == 3)
        {
			// This is a vector: we need to rotate it
            vtkDataArray* newAr = origAr->NewInstance();
            newAr->SetNumberOfComponents(3);
            newAr->SetName(origAr->GetName());
            if(origAr == normal)
            {
                // This is a normal
                transform->TransformNormals(origAr, newAr);
                normalToAdd = newAr;
            }
            else
            {
                // This is a standard vector
                transform->TransformVectors(origAr, newAr);
                toAdd.push_back(newAr);
            }          
            toDelete.push_back(ii);
        } 
    }
    
    // Delete the modified arrays
    for(int ii = toDelete.size()-1; ii >= 0; ii--)
    {
		data->RemoveArray(toDelete[ii]);
	}
	if(normalToAdd)
	{
		data->SetNormals(normalToAdd);
		normalToAdd->Delete();
	}
	for(auto it = std::begin(toAdd) ; it != std::end(toAdd) ; it++)
	{
		data->AddArray(*it);
		(*it)->Delete();
	}
}

}
