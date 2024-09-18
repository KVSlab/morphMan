"""Top-level package for morphMan."""

from importlib.metadata import metadata

from .manipulate_area import *
from .manipulate_bend import *
from .manipulate_bifurcation import *
from .manipulate_branch import *
from .manipulate_curvature import *
from .manipulate_surface import *

meta = metadata("morphman")
