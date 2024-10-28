import json

class Constants:
    def __init__(self, json_file):
        # Load the constants from the JSON file
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Dynamically set the attributes of the class
        for key, value in data.items():
            setattr(self, key, value)

# Initialize the constants object by loading from constants.json
constants = Constants('./satxplor/utils/constants.json')
