import vtk
from vtkmodules.util import numpy_support as ah


def cgnsReader(file):
    reader = vtk.vtkCGNSReader()
    reader.SetFileName(file)
    reader.UpdateInformation()
    reader.EnableAllCellArrays()
    reader.EnableAllPointArrays()
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
