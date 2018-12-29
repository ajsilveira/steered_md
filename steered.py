import os
import copy
import time
import pathlib
import logging
import argparse
import parmed as pmd

import mdtraj as md
import numpy as np

from simtk import openmm, unit
from simtk.openmm import app
from mdtraj.reporters import NetCDFReporter

from simtk.openmm import XmlSerializer

logger = logging.getLogger(__name__)
platform = openmm.Platform.getPlatformByName('CUDA')
integrator = openmm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2*unit.femtosecond)
#integrator.setRandomNumberSeed(129)

with open('openmm_system.xml', 'r') as infile:
    openmm_system = XmlSerializer.deserialize(infile.read())
pdbx = app.PDBxFile('mem_prot_md_system.pdbx')
positions = [(pos[0], pos[1], pos[2] + 1*unit.nanometers) for pos in pdbx.positions]
openmm_simulation = app.Simulation(pdbx.topology, openmm_system, integrator, platform)
####### Reading positions from mdtraj trajectory ##########################
#topology = md.Topology.from_openmm(pdbx.topology)
#t = md.load('traj.nc',top=topology)
#positions = t.openmm_positions(1999)
############################################################################
openmm_simulation.context.setPositions(positions)
mdtraj_reporter = NetCDFReporter('steered_forward.nc',500)
openmm_simulation.reporters.append(mdtraj_reporter)
openmm_simulation.reporters.append(app.StateDataReporter('state_forward.log', 500, step=True,
    potentialEnergy=True, temperature=True))

cvforce = []
for force in openmm_system_d.getForces():
    forceType = force.__class__.__name__
    if forceType == 'CustomCVForce':
        cvforce.append(force)

openmm_simulation.context.setParameter('K_parallel', 1250*unit.kilojoules_per_mole/unit.nanometer**2)
structure = pmd.openmm.load_topology(pdbx.topology, system=openmm_system, xyz=positions)
pe = pmd.openmm.energy_decomposition(structure, openmm_simulation.context)
logger.debug(pe)
# This is top ----> bottom
for i in range(0,160):
   lbda = i/159
   logger.debug('lambda ({})= {}'.format(i,lbda))
   openmm_simulation.context.setParameter('lambda_restraints', lbda)
   logger.debug('state, r_parallel, r_orthogonal')
   logger.debug('{}, {}, {}'.format(i,
                cvforce[0].getCollectiveVariableValues(openmm_simulation.context)[0],
                cvforce[1].getCollectiveVariableValues(openmm_simulation.context)[0]))
   for step in range(125):
       openmm_simulation.step(50)
       logger.debug('{}, {}, {}'.format(i, cvforce[0].getCollectiveVariableValues(openmm_simulation.context)[0],
                    cvforce[1].getCollectiveVariableValues(openmm_simulation.context)[0]))
positions = openmm_simulation.context.getState(getPositions=True).getPositions()
app.PDBxFile.writeFile(pdbx.topology,positions,open("final.pdbx", 'w'))
structure = pmd.openmm.load_topology(pdbx.topology, system=openmm_system_d, xyz=positions)
pe = pmd.openmm.energy_decomposition(structure, openmm_simulation.context)
logger.debug(pe)
