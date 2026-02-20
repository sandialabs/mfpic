import numpy as np
import xml.etree.ElementTree as ET

import vtk
import vtk.util.numpy_support as vtk_numpy_support

def read_timesteps_and_files():
    xml_tree = ET.parse("MeshOutput/MeshOutput.pvd")
    root = xml_tree.getroot()
    vtk_collection = root[0]

    timesteps = []
    filenames = []
    for data_set in vtk_collection:
        timestep = float(data_set.attrib['timestep'])
        timesteps.append(timestep)

        file = data_set.attrib['file']
        filenames.append(file)

    timesteps = np.array(timesteps)
    return timesteps, filenames

def read_single_timestep(filename):
    grid_reader = vtk.vtkXMLPUnstructuredGridReader()
    grid_reader.SetFileName(f"MeshOutput/{filename}")
    grid_reader.Update()
    output = grid_reader.GetOutput()
    points = vtk_numpy_support.vtk_to_numpy(output.GetPoints().GetData())

    timestep_data = {'points': points}

    data_on_points = output.GetPointData()
    for key in data_on_points.keys():
        timestep_data[key] = vtk_numpy_support.vtk_to_numpy(data_on_points[key])

    return timestep_data

def read_mesh_data():
    timesteps, filenames = read_timesteps_and_files()

    mesh_data = []
    for filename in filenames:
        mesh_data.append(read_single_timestep(filename))

    return timesteps, mesh_data
