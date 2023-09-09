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

Because of the dynamic nature of the new configuration files, they must now be **dynamically loaded**. This means that instead of using a `properties.py` variable like so:

```python
import properties.properties as prop
print(prop.dry_rocket_mass)
```

We now must load the config file into a python dictionary using the functions exposed in `properties/data_loader.py`. 

Here is a simple way to dynamically load a config file:

```python
import properties.properties as prop
import properties.data_loader as dataloader

config = dataloader.load_config(prop.sim_config)
```

After executing these lines, the config is now loaded in `config` and we are ready to access the config variables in our code.

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

`data_loader` will allow us to access the data within like a dictionary, where every node of the tree is accessed using Python's built in dictionary access. In other words, to retrieve the value stored in "field1", we would use the following (Assuming `config` has already been dynamically loaded.)

```python
config["object"]["field1"]      # Returns 123
```

This works very similarly with every other field.

```python
config["object"]["field3"]["foo"]      # Returns "bar"
config["justanode"]                    # Returns 3.14
config["object"]                       # Returns a dictionary with all subnodes
```

All objects are automatically casted to their expected type during parsing, but they can be manually type-casted using one of the Python `!flags` (explored below)

Adding data to a `.yaml` file is just as easy as adding another node to the tree, and it can immediately be accessed wherever `config` is defined using dictionary notation. As will be discussed later in documentation, however, it is better to specify structures for your data if not having that piece of data will cause the program to throw an exception.
##### Custom datatypes
Sometimes you need to include additional types in a config file's fields to make code more readable or to make it compatible with legacy code. If this is the case, you can add a constructor for those types within `data_loader` and a custom `!flag` to tell the data loader that this field is of a specific type. For user convenience, `data_loader` for GNC comes with 2 types already pre-implemented:

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

