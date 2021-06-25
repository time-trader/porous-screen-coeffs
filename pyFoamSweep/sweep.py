import os
import numpy as np
import pandas as pd

from PyFoam.RunDictionary.SolutionFile import SolutionFile
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.Execution.ConvergenceRunner import ConvergenceRunner
from PyFoam.LogAnalysis.BoundingLogAnalyzer import BoundingLogAnalyzer


u_mag = 10
rho = 1
porous_dimensions = [0.5, 1, 1]
#  aoa_sweep = range(-85, 86)

aoa_sweep=[-80, -60, -40, -20, 0, 20, 40, 60, 80]
theta=90

case = os.curdir
dire = SolutionDirectory(case, archive="AoA")

dire.addBackup("PyFoamSolve.logfile")
dire.addBackup("PyFoamSolve.analyzed")
dire.addBackup("postProcessing")
dire.addBackup(dire.initialDir())

u_initial = SolutionFile(dire.initialDir(), "U")

sweep_results = dire.makeFile("sweep_results")

sweep_results.writeLine("Psi", "Fx", "Fy", "Fz")

for aoa in aoa_sweep:

    psi = np.deg2rad(aoa)
    theta = np.deg2rad(theta)

    u_0 = [u_mag * np.cos(psi) * np.sin(theta), u_mag * np.sin(psi) * np.sin(theta), u_mag * np.cos(theta)]

    u_initial.replaceInternal("(%f %f %f)" % (u_0[0], u_0[1], u_0[2]))

    # Run the solver
    run = ConvergenceRunner(BoundingLogAnalyzer(), argv=["myporousSimpleFoam", "-case", case], silent=True)
    run.start()

    # Post-process:

    # Fetch some data
    p_data=pd.read_csv(case + '/postProcessing/pressure_inlet/0/surfaceFieldValue.dat',
                sep="\t", header=None, comment='#')
    v_data=pd.read_csv(case + '/postProcessing/velocity_outlet/0/surfaceFieldValue.dat',
                sep="\t", header=None, comment='#')

    # Parse velocity
    u_i = [float(x)/(porous_dimensions[1] * porous_dimensions[2]) for x in v_data.iloc[-1,1].strip('[(|)]').split()]

    # Dump forces csv
    f_x = p_data.iloc[-1, 1]
    f_y = (u_0[1] - u_i[1]) * rho * u_0[0] * (porous_dimensions[1] * porous_dimensions[2])
    f_z = (u_0[2] - u_i[2]) * rho * u_0[0] * (porous_dimensions[1] * porous_dimensions[2])
    sweep_results.writeLine((aoa, f_x, f_y, f_z))

    # Archive
    dire.lastToArchive("AoA_%g_deg" % (aoa))
    dire.clearResults(additional=['postProcessing'], verbose=True)


