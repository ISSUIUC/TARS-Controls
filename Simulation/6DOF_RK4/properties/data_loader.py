import yaml 
import numpy
import os

# Exception to handle errors raised by the config validator.
class ConfigurationException(Exception):
    pass

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

def rec_validate_config_section(traceback, cfg_section):
    """Recursively validates a config output by the yaml parser by checking every field in the dictionary
    for _None or _None:Required data fields, warning and erroring accordingly."""
    if type(cfg_section) is dict:
      # propagate recursion
      for field in cfg_section.keys():
        rec_validate_config_section(traceback+"."+field, cfg_section[field])
    else:
        if type(cfg_section) is str and cfg_section == "_None":
            print("PySim Configuration WARN: " + traceback + " is inherited from a structure but has an unset value.")
        if type(cfg_section) is str and cfg_section == "_None:Required":
            raise ConfigurationException("PySim Configuration ERR: " + traceback + " is inherited from a structure but has an unset:required value.")


def validate_config(parsed_yaml):
   # Top level for config validation
   sections = parsed_yaml.keys()
   for section in sections:
        if(section == "define"):
            continue
        rec_validate_config_section(section, parsed_yaml[section])

def load_config(config_path):
  """Loads and validates a config given by the config path (relative to main.py)
  Returns config data and constants as a python dictionary.
  """

  raw_yaml = ""

  with open(os.path.join(os.path.dirname(__file__), "../properties/typedef.yaml")) as typedef: 
    # Load yaml type definitions
    raw_yaml += typedef.read()


  with open(os.path.join(os.path.dirname(__file__), config_path)) as cfg_raw: 
    # Append the config file to the typedef header
    raw_yaml += "\n" + cfg_raw.read()
 
  yaml_content = yaml.load(raw_yaml, Loader=get_loader()) # PyYaml parser
  validate_config(yaml_content)
  return yaml_content