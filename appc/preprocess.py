import os, math
import numpy as np

import matplotlib.pyplot as plt

def close_TE(surface, tol=1e-8, force_TE_upper=False):
    
    """
    Returns a copy of the surface with a closed trailing edge
    """
    
    distance = math.sqrt( (surface[0, 0]-surface[-1, 0])**2 + (surface[0, 1]-surface[-1, 1])**2 )
    
    # Case that surface already has a closed TE
    if distance == 0:
        return surface
    
    # Case that surface already has a TE almost closed
    elif distance <= tol:
        
        surface = np.ndarray.copy(surface)
        
        if surface[-1, 0] >= surface[0, 0]:
            surface[0, :] = surface[-1, :]
        else:
            surface[-1, :] = surface[0, :]
            
        return surface
    
    # Close the TE
    # elif surface[-1, 0] >= surface[0, 0] and False:  # This causes randonmness to discretization so it's forced to be false
    elif force_TE_upper:
        return np.vstack((surface[-1], surface))
    
    else:
        return np.vstack((surface, surface[0]))



def plot_jetfoil(upper, center, lower, save_path=None, file_prefix="jetfoil"):

    # ------------- VISUALIZE SURFACE GRID ------------------------------
    fig = plt.figure()
    ax = fig.gca()
    
    # Plot contours
    ax.plot(upper[:, 0], upper[:, 1], ".-", label="Upper element")
    ax.plot(center[:, 0], center[:, 1], ".-", label="Centerbody")
    ax.plot(lower[:, 0], lower[:, 1], ".-", label="Lower element")

    nmarker_up = int(0.25 * upper.shape[0])
    nmarker_cb = int(0.25 * center.shape[0])
    nmarker_lo = int(0.25 * lower.shape[0])
    
    # Put markers at begining of contour to see start and orientation of contour
    ax.plot(upper[[0, nmarker_up], 0], upper[[0, nmarker_up], 1], "*k")
    ax.plot(center[[0, nmarker_cb], 0], center[[0, nmarker_cb], 1], "*k")
    ax.plot(lower[[0, nmarker_lo], 0], lower[[0, nmarker_lo], 1], "*k")
    
    # Format plot
    ax.set_xlabel(r"$x$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    
    ax.set_aspect("equal")
    ax.legend(loc="best", frameon=False, fontsize=10)
    [ax.spines[side].set_visible(False) for side in ["top", "right"]];
    
    fig.tight_layout()
    
    if save_path != None:
        for ext in [".png", ".svg"]:
            fig.savefig(os.path.join(save_path, file_prefix + "-geom" + ext), dpi=300, transparent=True)

    return fig, ax