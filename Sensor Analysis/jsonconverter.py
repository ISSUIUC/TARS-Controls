import json
import csv
import os
from typing import Callable, Any, NamedTuple, cast
import lark
from lark import Lark
from lark.visitors import Interpreter, v_args
from argparse import ArgumentParser
import struct
import pathlib
#######################################
# PLACE JSON IN SENSOR ANALYSIS FOLDERS
#######################################
def construct_header_string(jsonobject: dict[str, Any], path: str = "") -> str:
    """
    Given a json object that is representative of the whole file, recursively create a header string containing every
    key in the json as a period-separated path.
    :param jsonobject: A JSON object to construct a header string for
    :param path: The path this JSON object is at
    :return: A string for the keys in this json object
    """
    items: list[str] = []
    for key, value in jsonobject.items():
        if isinstance(value, dict):
            # If it's a dictionary, continue recursively but including the trace to the datatype
            items.append([construct_header_string(value, f"{path}.{key}")])
        else:
            temp = f"{path}.{key}"
            # idk = ''.join(char for char in temp if char not in ["[", "]", ","])
            items.append(temp)
    return items
def construct_data_string(jsonobject: dict[str, Any]):
    """
    Construct a string turning a single JSON object into csv columns
    :param jsonobject: A JSON object to convert
    :return: a string for the values in this json object as csv rows
    """
    items: list[str] = []
    for key, value in jsonobject.items():
        if isinstance(value, dict):
            new_item = [construct_data_string(value)]
            items.append(new_item)
        elif isinstance(value, list):
            # If it's a list, take only the first value
            # (Lists are only used for FSM, anyway, someone could alter this behavior though if they wish)
            # todo make this less fragile, but it should be good for now
            new_item = [str(value[0])]
            items.append(new_item)
        else:
            new_item = [str(value)]
            items.append(new_item)
    
    # return ",".join(items)
    return items
with open('Sensor Analysis\data122.json') as json_file:
    jsondata = json.load(json_file)
header_string = construct_header_string(jsondata[0])
# new csv file to be created
data_file = open('Sensor Analysis\csvfile.csv', 'w', newline='')
csv_writer = csv.writer(data_file)
csv_writer.writerow(header_string)
for data_object in jsondata:
        for i in data_object:
            i.replace("[","")
        csv_writer.writerow(construct_data_string(data_object))
data_file.close()
# Removing extra characters to get it to a readable state
input_file_path = 'Sensor Analysis\csvfile.csv'
# final csv file that will be used by sensor analysis
output_file_path = 'Sensor Analysis\csvfinal.csv'

# Characters to remove
characters_to_remove = ['[', ']', "'",'"']

# Open the input file in read mode and output file in write mode
with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    # Read the content of the input file
    file_content = input_file.read()

    # Remove characters "a" and "b" from the content
    filtered_content = ''.join(char for char in file_content if char not in characters_to_remove)

    # Write the filtered content to the output file
    output_file.write(filtered_content)
# Deleting intermediate file
if os.path.exists(input_file_path):
    os.remove(input_file_path)
