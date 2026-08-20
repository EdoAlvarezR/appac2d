import os
import math
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def plot_wake(foils, wakes, aoa=0.0, delta=0.0, title="",
              elements=["upper", "center", "lower"], save_path=None, file_prefix="jetfoil",
              ground_h=np.inf):
    
    # ------------- PLOT WAKE SHAPE ----------------------------------
    fig = plt.figure(figsize=np.array([10,5]))
    ax = fig.gca()
    
    fig.suptitle(title)
    
    # Plot contours
    for (eli, element) in enumerate(elements):
        
        # Fetch indices of panels in this element
        if element=="upper":
            indices = range(int(foils.m[0][0]), int(foils.m[0][0]+foils.m[0][1]))
            wakeindices = range(int(wakes.m[0][0]), int(wakes.m[0][0]+wakes.m[0][1]))
            
        elif element=="center":
            indices = range(int(foils.m[0][0]+foils.m[0][1]), int(foils.m[0][0]+foils.m[0][1]+foils.m[0][2]))
            
        else: 
            indices = range(0, int(foils.m[0][0]))
            wakeindices = range(0, int(wakes.m[0][0]))
            
        elementi = 1 if element=="upper" else 2 if element=="center" else 0
        
        # Fetch jetfoil contour
        xs = foils.xo[indices]
        ys = foils.yo[indices]
        
        # Fetch wake edge
        xs_wake = wakes.xo[wakeindices]
        ys_wake = wakes.yo[wakeindices]
        
        # Untwist APPC geometry
        if not np.isfinite(ground_h):
            aoa_rad = aoa * math.pi/180
            xs, ys = xs*math.cos(-aoa_rad) + ys*math.sin(-aoa_rad), -xs*math.sin(-aoa_rad) + ys*math.cos(-aoa_rad)
            xs_wake, ys_wake = xs_wake*math.cos(-aoa_rad) + ys_wake*math.sin(-aoa_rad), \
                                                -xs_wake*math.sin(-aoa_rad) + ys_wake*math.cos(-aoa_rad)
        
        # Plot jetfoil contour
        ax.plot(xs, ys, "-", color="black", linewidth=1.0, alpha=1.0, clip_on=False)
        ax.plot(xs_wake, ys_wake, "-", color="green", linewidth=1.0, alpha=1.0, 
                                        label="APPC" if eli==0 else "", clip_on=True)
        
        
    # Format plot
    ax.set_xlabel(r"$x$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    
    ax.set_xlim([-0.5, 3])
    ax.set_ylim([-0.75, 0.2])
    
    ax.set_aspect("equal")
    ax.legend(loc="best", frameon=False, fontsize=10)
    [ax.spines[side].set_visible(False) for side in ["top", "right"]]
    
    fig.tight_layout()
    
    # Save plot
    if save_path != None:
        for ext in [".png"]:
            fig.savefig(os.path.join(save_path, file_prefix + "-wake" + ext), dpi=300, transparent=True)

    return fig, ax

def process_surface_pressure(foils, wakes, xc, yc, Cp, aoa=0.0, 
              elements=["upper", "center", "lower"], save_path=None, file_prefix="jetfoil",
              export_pressure_csv=False, print_pressure_data=False,
              ground_h=np.inf):

    pressure_data = {}
    centerofpressure_data = {}
    geometry_data = {}
    
    # Plot each jetfoil element
    for element in elements:
    
        # Fetch indices of panels in this element
        if element=="upper":
            indices = range(int(foils.m[0][0]), int(foils.m[0][0]+foils.m[0][1]))
        elif element=="center":
            indices = range(int(foils.m[0][0]+foils.m[0][1]), int(foils.m[0][0]+foils.m[0][1]+foils.m[0][2]))
        else: 
            indices = range(0, int(foils.m[0][0]))
            
        elementi = 1 if element=="upper" else 2 if element=="center" else 0
            
        # Fetch jetfoil contour
        xs = foils.xo[indices].flatten()
        ys = foils.yo[indices].flatten()
        
        # Untwist jetfoil contour
        if not np.isfinite(ground_h):
            aoa_rad = aoa * math.pi/180
            xs, ys = xs*math.cos(-aoa_rad) + ys*math.sin(-aoa_rad), \
                                                -xs*math.sin(-aoa_rad) + ys*math.cos(-aoa_rad)
        
        # Fetch APPC2D data
        xs_appc = np.array(xc.tolist()[elementi][0]).flatten()
        ys_appc = np.array(yc.tolist()[elementi][0]).flatten()
        Cps_appc = np.array(Cp.tolist()[elementi][0]).flatten()
        
        # Untwist APPC geometry
        if not np.isfinite(ground_h):
            xs_appc, ys_appc = xs_appc*math.cos(-aoa_rad) + ys_appc*math.sin(-aoa_rad), \
                                            -xs_appc*math.sin(-aoa_rad) + ys_appc*math.cos(-aoa_rad)

        # Store pressure locations per element and optionally export to CSV.
        pressure_df = pd.DataFrame({
            "x (m)": xs_appc,
            "y (m)": ys_appc,
            "Cp": Cps_appc,
        })
        pressure_data[element] = pressure_df

        if save_path is not None and export_pressure_csv:
            csv_name = f"{file_prefix}-Cp-{element}.csv"
            pressure_df.to_csv(os.path.join(save_path, csv_name), index=False)

        # Calculate center of pressure
        xs_midpoint = (xs[1:] + xs[:-1]) / 2
        ys_midpoint = (ys[1:] + ys[:-1]) / 2
        dxs = xs[1:] - xs[:-1]
        dys = ys[1:] - ys[:-1]
        Cps = (Cps_appc[1:] + Cps_appc[:-1]) / 2

        Cfx = sum(Cps*dxs)
        Cfy = sum(Cps*dys)

        x_centerofpressure = sum(xs_midpoint * Cps*dxs) / Cfx
        y_centerofpressure = sum(ys_midpoint * Cps*dys) / Cfy

        centerofpressure_data[element] = (np.array([x_centerofpressure, y_centerofpressure]), np.array([Cfx, Cfy]))
        
        if print_pressure_data:
            print()
            print("#"*73)
            print(f'#\t"{element}" Element Pressure Data')
            print("#"*73)
            print()
            print(f'Center of pressure (x, y) = ({x_centerofpressure}, {y_centerofpressure})')
            print()
            print(pressure_df.to_csv(index=False))

        geometry_data[element] = {'xs': xs, 'ys': ys, 'xs_appc': xs_appc, 'ys_appc': ys_appc, 'Cps_appc': Cps_appc}

    return pressure_data, centerofpressure_data, geometry_data

def plot_surface_pressure(*args, title="",
              elements=["upper", "center", "lower"], save_path=None, file_prefix="jetfoil",
              return_pressure_data=False, export_pressure_csv=True,
              invert_y=True, plot_centerofpressure=False, 
              ground_h=np.inf, **kwargs):


    # ------------- PLOT SURFACE PRESSURE ---------------------------------

    figs = []
    axs = []

    (pressure_data, centerofpressure_data, 
         geometry_data) = process_surface_pressure(*args, elements=elements, save_path=save_path, 
                                                           file_prefix=file_prefix, export_pressure_csv=export_pressure_csv,
                                                           **kwargs)
    
    # Plot each jetfoil element
    for element in elements:

        d = geometry_data[element]
        xs, ys, xs_appc, ys_appc, Cps_appc = d['xs'], d['ys'], d['xs_appc'], d['ys_appc'], d['Cps_appc']

        (centerofpressure, Cf) = centerofpressure_data[element]
        
        # Initialize plot
        fig = plt.figure(figsize=np.array([7*1.3, 5*0.8])*2/3)
        ax = fig.gca()
        ax2 = ax.twinx()
        
        fig.suptitle(title)
        
        # Plot
        ax.plot(xs_appc, Cps_appc, "-", alpha=0.8, label="APPC", clip_on=False)
        ax.plot(xs_appc[[0, 10]], Cps_appc[[0, 10]], "*r", alpha=0.25, clip_on=False)
        ax.plot(xs_appc[[0, 10]], Cps_appc[[0, 10]], "*r", alpha=0.25, clip_on=False)
        
        # Plot jetfoil contour
        ax2.plot(xs, ys - (ground_h if np.isfinite(ground_h) else 0), ":", color="black", linewidth=1.0, alpha=1.0, clip_on=False)
        ax2.plot(xs[[0, 10]], ys[[0, 10]] - (ground_h if np.isfinite(ground_h) else 0), "*r", alpha=0.25, clip_on=False)

        # Plot center of pressure
        if plot_centerofpressure:
            ax2.plot(centerofpressure[0], centerofpressure[1], "x", color="black", linewidth=1.0, alpha=1.0, clip_on=False)
    
        # Format twin axis
        ax2.set_aspect(1.0)
        yaux = 0.65
        if element=="lower":
            ax2.set_ylim(np.array([-yaux, 0]) + 0.25)
        elif element=="center":
            ax2.set_ylim(np.array([0, yaux]) - 0.20)
        else:
            ax2.set_ylim(np.array([0, yaux]) - 0.35)
        ax2.set_yticks([])
    
        # Beautify the plot
        xlims, dx = [-0.2, 1.0], 0.2
        ylims, dy = [-2, 1], 1.0
        
        if element in np.array(["upper", "center", "lower"])[[0]]: 
            #ylims = [ylims[1], ylims[0]]
            dy *= -1
            
        ax.set_xlim(xlims)
        ax.set_xticks(np.arange(xlims[0], xlims[1], dx))
        #ax.set_ylim(ylims)
        #ax.set_yticks(np.arange(ylims[0], 1.01*ylims[1], dy))

        if invert_y and element=="lower":
            ax.invert_yaxis()
    
        if not element=="lower":
    
            for a in [ax, ax2]:
                a.spines["bottom"].set_visible(False)
                a.set_xticks([])
                
        else:
            ax.set_xlabel(r"$x$-position (m)")
            
        ax.set_ylabel(r"Pressure $C_p$")
    
        for a in [ax, ax2]:
            a.spines["right"].set_visible(False)
            a.spines["top"].set_visible(False)
            
        if element=="center":
            
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1],
                           loc="center left", bbox_to_anchor=(0.5, 0.3), frameon=False, fontsize=10)
    
        # ax2.plot(xlims, zeros(2), ":k", alpha=0.25, linewidth=1, clip_on=false)
    
    
        # fig.tight_layout
        faux = 0.0
        fig.tight_layout(rect=[0.0, 0.0, 1.0, 1.0])
    
        # Save plots
        if save_path != None:
            for ext in [".png"]:
                fig.savefig(os.path.join(save_path, file_prefix + "-Cp-" + element + ext), dpi=300, transparent=True)

        figs.append(fig)
        axs.append((ax, ax2))

    if return_pressure_data:
        return figs, axs, pressure_data, centerofpressure_data, geometry_data

    return figs, axs