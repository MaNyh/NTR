"""
Script for writing files
 - vtk object in .vtk or .vtu
 - numpy array in h5-Format like .h5, .hdf5, .cgns

Author: Lovisa Wonnemann
Thieve and editor: Malte Nyhuis (Thanks, Lovis!)
Date: Sep 2021
"""

import vtk
from vtkmodules.util import numpy_support as ah
#import h5py
import numpy as np
import pyvista as pv

# https://blog.kitware.com/hdf5-reader-and-writer-for-unstructured-grids/

def cgnsReader(file):
    reader = vtk.vtkCGNSReader()
    reader.SetFileName(file)
    reader.UpdateInformation()
    reader.EnableAllCellArrays()
    reader.EnableAllPointArrays()
    reader.Update()

    appendFilter = vtk.vtkAppendFilter()
    # Points with same coordinates are merged
    # with tolerance 0.0000001 GMC GLOBAL Properties
    appendFilter.MergePointsOn()
    appendFilter.SetTolerance(0.0000001)
    # Set Base range

    # loop through all bases
    for nameBaseId in range(reader.GetNumberOfBaseArrays()):
        base = reader.GetOutput().GetBlock(nameBaseId)
        # loop through all zones in base
        for nameZoneId  in range(base.GetNumberOfBlocks()):
            zone = base.GetBlock(nameZoneId)
            # merge zone data and update
            appendFilter.AddInputData(zone)
            appendFilter.Update()

    mergeZones = appendFilter.GetOutput()
    # Add to data container

    return mergeZones

def vtkUnstructuredGridReader(file):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file)
    reader.UpdateInformation()
    reader.Update()
    dataset = reader.GetOutput()
    return dataset

def vtkFLUENTReader(file):
    reader = vtk.vtkFLUENTReader
    reader.SetFileName(file)
    reader.UpdateInformation()
    reader.Update()
    dataset = reader.GetOutput()
    return dataset

def gradient(dataset, array_name, result_array_name):
    """
    #Berechnung des Gradienten mithilfe des Filters "vtkGradientFilter".
    #Die Spalten des Arrays Gradient sind: du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, dw/dx, dw/dy, dw/dz.
    #Der Funktion muss die zusammengefügte Zone (dataset), der Name der Größe,
    #dessen Gradient berechnet werden soll sowie der Name, unter der der
    #berechnete Gradient in Base abgespeichert werden soll, angegeben werden.
    #Die Ausgabe ist die Zone (dataset), die den hinzugefügten Gradienten besitzt.
    """
    gradient = vtk.vtkGradientFilter()
    gradient.SetInputData(dataset)
    gradient.SetInputScalars(FIELD_ASSOCIATION, array_name)
    gradient.SetResultArrayName(result_array_name)
    gradient.Update()
    computed = gradient.GetOutput()
    return computed


def get_array(dataset, array_number):
    """
    #Umwandlung des vtk-Arrays un ein Numpy-Array.
    #Die Arrays sind durch Nummern abgespeichert. Also z.B. Array0 = Density ...
    #Die Nummer des Arrays muss angegeben werden.
    """
    array = ah.vtk_to_numpy(dataset.GetCellData().GetArray(array_number))
    return array


def add_array(dataset, array, array_name):
    """
    #Wandelt einen Numpy-Array in einen vtk-Array um und fügt diesen in die Zone ein.
    #Als Eingabe ist die Zone, der numpy_Array und der Name, unter dem der Array
    #abgespeichert werden soll nötig. Die Ausgabe ist die Zone mit dem hinzugefügten Array.
    """
    vtk_arr = ah.numpy_to_vtk(array)
    vtk_arr.SetName(array_name)
    dataset.GetCellData().AddArray(vtk_arr)
    return dataset


def get_neighbor(dataset, cellId):
    """
    #cellId: #ZellenId, dessen Nachbarn gefunden werden soll
    """
    neighborCellIdListAll = vtk.vtkIdList()  # Id-Liste von den Nachbarn der Zelle
    number_of_neighbors = 0
    for i in range(
        dataset.GetCell(cellId).GetNumberOfFaces()):  # Für jede Fläche der Zelle soll der Nachbar gefunden werden
        face = dataset.GetCell(cellId).GetFace(i)  # Face ist eine Fläche der Zelle
        pointIdList = face.GetPointIds()  # Um die Nachbarn zu finden, müssen die Punkte-IDs der Fläche in pointIdList gespeichert werden
        neighborCellIdList = vtk.vtkIdList()  # Erstellung einer Liste, in der die ID des Nachbarn mit der gleichen Fläche gespeichert wird
        # Bestimmung der Nachbarn, cellId: Zellen-Id, desssen Nachbar gefunden werden soll,
        # pointIdList: Liste der Punkte-Ids der Fläche, neighborCellIdList: Liste in der die Id der Nachbarzelle gespeichert wird
        dataset.GetCellNeighbors(cellId, pointIdList, neighborCellIdList)
        # for j in range(neighborCellIdList.GetNumberOfIds()):
        if neighborCellIdList.GetNumberOfIds() == 1:
            # Der Nachbar jeder Fläche wird hier in die Gesamtliste der Nachbarn hinzugefügt
            # neighborCellIdListAll.InsertNextId(neighborCellIdList.GetId(j))
            neighborCellIdListAll.InsertNextId(neighborCellIdList.GetId(0))
            number_of_neighbors = number_of_neighbors + 1
        else:
            neighborCellIdListAll.InsertNextId(-1)

    return neighborCellIdListAll, number_of_neighbors


def get_points(dataset):
    """
    #Die Punkte eines Datensatzes werden als numpyArray ausgegeben
    """
    vtk_array = dataset.GetPoints().GetData()
    numpy_array = ah.vtk_to_numpy(vtk_array)
    return numpy_array


def dict_neigbour(dataset):
    cellIdListAll = []
    number_of_neighbors_list = []
    for cellId in range(dataset.GetNumberOfCells()):
        (neighborCellIdList, number_of_neighbors) = get_neighbor(dataset, cellId)
        number_of_neighbors_list.append(number_of_neighbors)
        idList = []
        for i in range(neighborCellIdList.GetNumberOfIds()):
            idList.append(neighborCellIdList.GetId(i))
        cellIdListAll.append(idList)
    return cellIdListAll, number_of_neighbors_list


def get_cell_center(dataset):
    center_filter = vtk.vtkCellCenters()
    center_filter.SetInputData(dataset)
    center_filter.Update()
    center = center_filter.GetOutput()
    points_center = get_points(center)
    return points_center


class writerH5:
    def __init__(self, filename=None):
        self.filename = filename
        if filename is not None:
            self.__setWriter__()

    def setFilename(self, filename: str):
        """
        Define output filename and init the writer
        :param filename: output filename like this path/name.format (type string)
        """
        # Set filename and format
        self.filename = filename
        # Init the writer
        self.__setWriter__()

    def __setWriter__(self):
        """
        Initialization the file and create first group
        """
        # Open the new file
        self.fileOut = h5py.File(self.filename, 'w')

    def createGroup(self, groupName: str, groupPath=None, setData=False, dataName='data', data=None, shape=None,
                    dataType=None):
        """
        Create a group optional set data
        :param groupPath: path for the new group e.g. /root node/slices/ (type string)
        :param groupName: name of the new group (type string)
        :param setData: flag for set data, default False (type bool)
        :param dataName: name for the data (type string)
        :param data: data (e.g. np.array or string)
        :param shape: data shape
        :param dataType: data type e.g. np.array -> np.uint8, strings -> 'utf-8'
        """

        if groupPath is None:
            group = self.fileOut.create_group(groupName)
            groupPath = groupName
        else:
            # Get the parent group
            group = self.fileOut.get(groupPath)
            # create a sub group
            group.create_group(groupName)
            groupPath = groupPath + '/' + groupName

        # Write data in sub group
        if setData is True:
            self.addDataToGroup(groupPath=groupPath, data=data, dataName=dataName, dataType=dataType, shape=shape)

    def addDataToGroup(self, groupPath: str, data, dataName='data', dataType=None, shape=None):
        """
        Set the data for specific group
        :param groupPath: path for the new group e.g. /root node/slices/ (type string)
        :param data: data (e.g. np.array or string)
        :param dataName: name for the data (type string)
        :param dataType: data type e.g. np.array -> np.uint8, strings -> 'utf-8'
        :param shape: data shape
        """
        # Get the group
        group = self.fileOut.get(groupPath)
        # Check if data is string
        if type(data) is str:
            dataType = h5py.string_dtype(encoding=dataType)
        # add data to group
        group.create_dataset(dataName, data=data, dtype=dataType, shape=shape)

    def removeAnything(self, path):
        """
        Deleting group or data
        :param path: path of anything, which should be delete (type string)
        """
        del self.fileOut[path]

    def templateBaseDict(self, baseDict: dict, groupName='meta data', groupPath=None):
        """
        This is a template for zone or zone merge arrays from cgnsReader
        :param baseDict: Dict with surface and flow solution keys, data type is np.array
        :param groupName: name of group, where the base will be storage (type string)
        :param groupPath: path for the new group e.g. /root node/slices/ (type string)
        """
        # Create the new groups
        if groupPath is None:
            group = self.fileOut.create_group(groupName)
            groupPath = groupName
        else:
            # Get the parent group
            group = self.fileOut.get(groupPath)
            # create a group
            group.create_group(groupName)
            groupPath = groupPath + '/' + groupName

        # Loop through base dict and read it, key [0] und value [1]
        for base in list(baseDict.items()):
            # Create the new group
            group.create_group(base[0])
            # Loop through zone dict and read it, key [0] und value [1]
            for zone in list(base[1].items()):
                self.templateSliceDict(sliceDict=zone[1], groupName=zone[0], groupPath=groupPath + '/' + base[0])

    def templateSliceDict(self, sliceDict: dict, groupName='slice', groupPath=None):
        """
        This is a template for zone or zone merge arrays from cgnsReader
        :param sliceDict: Dict with surface and flow solution keys, data type is np.array
        :param groupName: name of group, where the slice will be storage (type string)
        :param groupPath: path for the new group e.g. /root node/slices/ (type string)
        """
        # Create the new groups
        if groupPath is None:
            sliceGroup = self.fileOut.create_group(groupName)
        else:
            # Get the parent group
            group = self.fileOut.get(groupPath)
            # create a sub group
            sliceGroup = group.create_group(groupName)

        coordGroup = sliceGroup.create_group('GridCoordinates')
        # Loop through dict and read surface data, key [0] und value [1]
        for coord in list(sliceDict['surface'].items()):
            # Set data
            coordGroup.create_dataset(coord[0], data=coord[1], shape=coord[1].shape, dtype=np.uint8)
        # Create the new group
        solutionGroup = sliceGroup.create_group('FlowSolution')
        # Loop through dict and read floe solution data, key [0] und value [1]
        for solution in list(sliceDict['flow solution'].items()):
            # Set data
            solutionGroup.create_dataset(solution[0], data=solution[1], shape=solution[1].shape, dtype=np.uint8)

    def save(self):
        """
        Close the file to save it
        """
        self.fileOut.close()


class writerVTK:
    def __init__(self, data=None, filename=None):
        # Init file
        self.vtkObj = data
        self.filename = filename
        if data is not None:
            # Transform vtk to pyvista grid
            self.__getGrid__()
        else:
            self.grid = None

    def __getGrid__(self):
        """
        Prepare the input data for plotting
        Select the right grid type dependent on mesh type
        Define: vtkStructuredGrid, vtkUnstructuredGrid and vtkMultiBlockDataSet
        :return: grid
        """
        if str(type(self.vtkObj)).find('vtkStructuredGrid') != -1:
            self.grid = pv.StructuredGrid(self.vtkObj)
        elif str(type(self.vtkObj)).find('vtkUnstructuredGrid') != -1:
            self.grid = pv.UnstructuredGrid(self.vtkObj)
        elif str(type(self.vtkObj)).find('vtkMultiBlockDataSet') != -1:
            self.grid = pv.MultiBlock(self.vtkObj)
        else:
            print('Mesh type is not define.')
            self.grid = None

    def setFilename(self, filename: str):
        """
        Define output filename and specific output format
        :param filename: output filename like this path/name.format (type string)
        """
        # Set filename and format
        self.filename = filename

    def setData(self, data):
        """
        Set the data and transform it to a grid, that can be save
        E.g. vtkStructuredGrid, vtkUnstructuredGrid or vtkMultiBlockDataSet
        :param data: a vtk object with points and cell data
        """
        self.vtkObj = data
        # Transform vtk to pyvista grid
        self.__getGrid__()

    def save(self):
        """
        Save the data
        """
        self.grid.save(self.filename)
