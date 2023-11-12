Before Avionics ticket 1008, all simulation configuration was done through a single file called `properties.py`, which exported all necessary variables by just being included in the header of every file which needed to access simulation-specific variables.

To expand the abilities of Pysim, we have now pivoted away from using only `properties.py`, and instead the new system uses Pyyaml to parse `.yaml` files and uses constants from these new config files instead. This gives us multiple advantages:

`.yaml` files are easy to read and understand, even for someone who isn't familiar with programming in general. Additionally, they allow us to implement config validation, ensuring that instead of giving a cryptic `KeyError`, we can now get descriptive messages about what simulation variables are missing in particular.

### PySim data loading

We have retained `properties.py` to include constants which are calculated using language features, as well as to specify which configuration file should be used.

```python
...
# Config file
sim_config = "../properties/configs/test.yaml"
...
# Gravitational const
G = 6.6743*10**(-11)
# mass of earth
m_e = 5.9722*10**24
# radius of earth
r_e = 6.3781*10**6
...
```

It is important to note that the `.yaml` parser for Pysim *does* include the ability to also calculate constants on the fly (`!eval`), but that will be covered in a later section.

Because of the dynamic nature of the new configuration files, they have slight differences in access and loading from the old `properties.py` file. Instead of accessing a property like so:

```python
import properties.properties as prop
print(prop.dry_rocket_mass)
```

We access the config file for the current process. To do so, import `data_loader.py`:

```python
import properties.data_loader as dataloader
```

And access config variables through Python's dictionary access operation: `['key_name']`. The specifics of YAML config structuring are defined below.

**NOTE `9/10/2023`:** Current convention for Pysim indicates that there should be one simulation run per execution of the program. Because of this, `data_loader` will print a warning if the dynamic load function `load_config()` is called more than once in the same process. However, it will still function as intended.

##### Accessing properties defined in .yaml:
The spec for **YAML** can be found here: https://yaml.org/

Given a YAML tree given by the following
```yml
object:
	field1: 123
	field2: 456
	field3:
		foo: "bar"
justanode: 3.14
```

`data_loader.config` will allow us to access the data within like a dictionary, where every node of the tree is accessed using Python's built in dictionary access. In other words, to retrieve the value stored in "field1", we would use the following:

```python
dataloader.config["object"]["field1"]      # Returns 123
```

This works very similarly with every other field.

```python
dataloader.config["object"]["field3"]["foo"]      # Returns "bar"
dataloader.config["justanode"]                    # Returns 3.14
dataloader.config["object"]                       # Returns a dictionary with all subnodes
```

All objects are automatically casted to their expected type during parsing, but they can be manually type-casted using one of the Python `!flags` (explored below). It should be noted that it is **not possible** (by design) to access type data from this dictionary.

Adding data to a `.yaml` file is just as easy as adding another node to the tree, and it can immediately be accessed wherever `config` is defined using dictionary notation. As will be discussed later in documentation, however, it is better to specify structures for your data if not having that piece of data will cause the program to throw an exception.
##### Custom datatypes
Sometimes you need to include additional types in a config file's fields to make code more readable or to make it compatible with legacy code. If this is the case, you can add a constructor for those types within `yaml_datatypes.py` and a custom `!flag` to tell the data loader that this field is of a specific type. For user convenience, `yaml_datatypes` for GNC comes with 2 types already pre-implemented:

`!numpy/array` or `!array` will be typecasted into `numpy.array` on load.
`!eval` will execute the expression given and return it into the field.

```yml
object:
	my_arr: !numpy/array [1, 2, 3]   # This will be typecasted to numpy.array on load.
	calculated_value: !eval 1 + 2    # This will return the number 3.
```

##### YAML Flags
By now there have been several mentions of **flags**. YAML flags are written using the convention `!flag_name`. The Pyyaml parser provides multiple flags for typecasting:

```yml
my_object:
	my_string: !python/str 123       # Casted to "123"
	my_tuple: !python/tuple [1,2]    # Casted to a python tuple
	my_long: !python/long 5          # Casted to a long
```

While the Pyyaml parser does its best to properly decode types, you should use flags if you want to be absolutely sure that a key in the config is of a specific type.

## Config validation
To make sure that we give comprehensive errors whenever an expected value is missing in the config, `data_loader` provides the capability to validate configs that are loaded for the simulation.

**Note:** Config validation is completely optional, but it is highly recommended as it gives more descriptive errors than not using validation.

Whenever a config is loaded, it is prepended with the `typedef.yaml` header, which defines the structure for a config that we know is valid. It is built up of 3 main pieces of data:

```yml
*unset ==> Not having this variable set will result in a warning.
*optional ==> Not having this variable set will add a key with _Ignore to the dict.
*unset_required ==> Not having this variable set will result in an error.
```

Whenever a variable has its value set to one of these, it will be handled by `data_loader` as defined above.

To make validation easy, `typedef` is also meant to define the structure of either a config file as a whole, or any part of it.

It is important to note that for config validation to work properly, you must include the structure of the type you are trying to define in your configuration file from `typedef`, which can be done like so:

```yml
<<: *type_name
```

Using this directive will include all fields (and sub-fields!) of the structure. For instance, given a type `my_type` defining fields `hello: 2` and `world: "foo"`, inserting the following into a config:

```yml
my_obj:
	<<: *my_type
```

is equivalent to the config below:

```yml
my_obj:
	hello: 2
	world: "foo"
```

By using a combination of structures, as well as the predefined types `*unset`, `*optional`, and `*unset_required`, we can build up structures that are automatically validated when they are loaded.

Listed below are some of the pre-defined structures in `typedef.yaml`, as of `9/8/2023`.

`*simulated_rocket` - Collection for all rocket data excluding data about the motor.
`*simulated_motor` - Collection of data for a given motor.
`*simulated_body` - Contains center of mass/pressure data for a rocket body.
`*simulated_flaps` - Contains the required data for the flaps actuation simulation.
`*pysim_simulation` - Combines the above types into a full simulation config file.

Ideally, the structure definition for `pysim_simulation` would be such that not having any warnings/errors during data loading will guarantee a successful simulation run. As of `9/8/2023`, this is not guaranteed due to uncertainty regarding which sensors will be used and other simulation data. However, it is always good practice for the first directive in a config file to be

```yml
<<: *pysim_simulation
```

as it will at least help you get most of the way there.

### Outdated configuration files
If, while loading a configuration file, you recieve the following error: 
```
PySim Configuration ERR: The configuration file loaded was marked with the OUTDATED_CONFIG flag, and the simulation will crash if it continues with this config. Please select a different configuration file.
```

This means that the current version of the simulator does not support your config, and someone has marked the config as `OUTDATED`. If you know of a config which will not work with current PySim code, and you do not with to delete the config, simply mark it as `OUTDATED` by prepending the file with this code:

```
<<: *OUTDATED_CONFIG
---

...rest of the config here...
```


### Object templates
If you have an object of a certain type that you know you will be re-using (such as a specific COTS motor), it is better to add it to **templates** instead of re-writing the data for that object in each configuration file.

Object templates are stored in various files in `./templates`, and all `.yaml` files within the directory are loaded into the YAML parser at load-time. This means that any type or template you define within `./templates` will also be available in your config file.

It is important to note that all templates are loaded **after** `typedef.yaml`, so attempting to reference a template from within `typedef` will always throw an error. Also, referencing types from other template files is not recommended, as the order in which templates are loaded is **not guaranteed.**

To create an object template, select a relevant file (`motors.yaml` for motors, for instance), or create a new file for a new category of object. Then, define the object as a YAML type using the `define` directive:

```YAML
define: &cool_srad_motor
```

Then, inherit from any already existing structure or add new fields as is seen fit:

```YAML
define: &cool_srad_motor
	<<: *simulated_motor
	impulse: 30000.0
	motor_mass: 17.0
	delay: 0
	motor_lookup_file: "../lookup/srad_motor.csv"
	cm: !numpy/array [0.3755, 0., 0.]
```

You may then use this template instead of overriding `*simulated_motor` in your configuration file:

```YAML
> in test_config.yaml

...
motor:
	<<: *cool_srad_motor
```
### Advanced: Adding custom data types
To define a custom datatype, you must be able to retrieve the data given in the `.yaml` file (usually in string, sequence, or number form) and construct that data type with it.

Shown below is the implementation for the flag `!numpy/array`:
```python
def numpy_array_constructor(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode):
  """Create a numpy array from a yaml !numpy/array tag. (!numpy/array [1, 2, 3])"""
  return numpy.array(loader.construct_sequence(node))
```

Because the `numpy.array` can be constructed by just passing in a sequence, we can retrieve the sequence using `SafeLoader`'s `construct_sequence` function to retrieve the sequence associated with the node, and return a new numpy array. 

After a constructor is created for a data type, we just need to register it in `data_loader` within the `get_loader()` function by adding the following line:

```python
loader.add_constructor("!your_custom_type", your_constructor)
```

