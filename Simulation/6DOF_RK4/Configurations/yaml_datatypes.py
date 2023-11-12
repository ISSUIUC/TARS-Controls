import yaml 
import numpy

# This file defines all custom yaml datatypes that we wish to include in the YAML processor.

def numpy_array_constructor(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode) -> numpy.array:
  """Create a numpy array from a yaml !numpy/array tag. (!numpy/array [1, 2, 3])"""
  return numpy.array(loader.construct_sequence(node))

def evaluate_expr(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode):
  """Evaluates a python expression passed to the !eval tag. (!eval 1 + 2)"""
  return eval(loader.construct_scalar(node))

def get_loader():
  """Add the data parsers above into the PyYAML parser."""
  loader = yaml.SafeLoader
  loader.add_constructor("!numpy/array", numpy_array_constructor)
  loader.add_constructor("!array", numpy_array_constructor)
  loader.add_constructor("!eval", evaluate_expr)
  return loader
