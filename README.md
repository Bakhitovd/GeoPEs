# Geology and Petroleum Engineering
GeoPEs is a Python-based tool designed to perform specialized geophysical calculations based on user-defined configurations. The repository is structured to ensure seamless data processing, from reading configurations to outputting results in an Excel format.

## Key Components:
- calc_row.py:
Responsible for performing calculations on individual data rows based on user-defined configurations and functions.
- calculator.py:
Houses the Calculator class, a core component that orchestrates the calculation process and outputs results to an Excel file.
- config.py:
Introduces the Config class, a robust configuration handler designed to read and provide configuration settings throughout the application.
- functions.py:
A utility module packed with essential functions that assist in retrieving values, formulas, and other functionalities vital for the calculation process.
- main.py:
The application's entry point. It integrates various modules and classes, driving the primary operations of the tool.
- read_config.py:
Dedicated to reading and parsing configuration files, ensuring that user-defined settings are accurately captured and utilized.
