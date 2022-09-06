#!/usr/bin/env python
"""
This is modified from Rusty's example script in stompy/examples/dwaq_adjust_bc_fluxes.py 

Command-line tool to infer missing boundary fluxes, and update them
in a flow file.

This applies to non-flow BCs in DFM, i.e. sea boundary conditions and
discharges.

"""
from __future__ import print_function

import argparse
import sys,os
import numpy as np

import stompy_rusty.model.delft.waq_scenario as waq

# Check mass conservation and if not make corrections to the velocity field. 
parser = argparse.ArgumentParser(description='Adjust DFM D-WAQ output to add missing BC fluxes.')

parser.add_argument('hyd_fn', metavar='dfm_out.hyd', type=str,
                    help='path to hyd file')

args = parser.parse_args()

hyd_fn=args.hyd_fn

print("Opening hydro from %s"%hyd_fn)
hydro=waq.HydroFiles(hyd_fn)

print("Adjusting flows, updating %s in place"%hydro.get_path('flows-file'))

hydro.adjust_boundaries_for_conservation()

print("Done")


# To simply check mass conservation without correction
#bc_exchs=np.nonzero( hydro.pointers[:,0]<0 )[0]
#bc_adj_segs=hydro.pointers[bc_exchs,1]-1 # make 0-based
#tidx_select = slice(None) # choose the entire domain. 
#summary=hydro.check_volume_conservation_incr(seg_select=bc_adj_segs,
#                                            tidx_select=tidx_select,
#                                            verbose=True)
                                     
