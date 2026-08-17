import os
import numpy as np

import pyvista as pv


def vis2vtk(visout, aoa=0.0, save_path=None, file_prefix="jetfoil", untwist=False, paraview=True):

    # ------------- OUTPUT VTK WITH FLOW FIELD --------------------------
    
    points = np.column_stack((visout.VTX, np.zeros(visout.VTX.shape[0])))
    cells = np.column_stack((np.full(visout.TRI.shape[0], 3), visout.TRI - 1))
    cells = cells.astype(np.int64)

    mesh = pv.PolyData(points, cells.ravel())

    for (field_name, field_data) in visout.data.items():
        mesh.point_data[field_name] = field_data

    mesh.point_data["U"] = np.column_stack((visout.data["u"], visout.data["v"], np.zeros(points.shape[0])))
    
    if untwist:
        mesh.rotate_z(aoa, inplace=True)

        # TODO: ROTATE U FIELD AS WELL!!!!!
        print("TODO: ROTATE U FIELD AS WELL!!!!!")

    if save_path != None:
        file_name = os.path.join(save_path, file_prefix + "-appc.vtk")
        
        mesh.save(file_name)
        
        if paraview:
            os.system(f"paraview --data={file_name};")

    else:
        file_name = ""

    plotter = pv.Plotter()
    plotter.add_mesh(mesh, scalars="p", lighting=False)
    plotter.add_arrows(points, mesh.point_data["U"], mag=0.03, color="black")

    plotter.show()

    return mesh, file_name