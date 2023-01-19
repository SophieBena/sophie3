#!/usr/bin/env
""" file:thermoengine/__init__.py
    author: Aaron S. Wolf; Mark S. Ghiorso
    date: Tuesday June 27, 2017

    description: Python package interface to the PhaseObjC library.
"""
# Load all core methods and place in thermoengine namespace
from thermoengine import core
from thermoengine.core import *
#from thermoengine.core import chem
from thermoengine.chemistry import OxideMolComp, OxideWtComp, Oxides, ElemMolComp, Comp
from thermoengine.phases import Phase


from thermoengine import phases
from thermoengine import samples
from thermoengine import model
from thermoengine import calibrate
from thermoengine import equilibrate
from thermoengine import chemistry

__all__ = [s for s in dir() if not s.startswith('_')]
