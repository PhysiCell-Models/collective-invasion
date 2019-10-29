# Distributed Automated Parameter Testing (DAPT)

## Install
To install dowload the [repo](https://github.com/BenSDuggan/DAPT) or clone it on your machine `git clone https://github.com/BenSDuggan/DAPT`.  Once downloaded navigate to the root of the project (DAPT) and run `pip install -r requirements.txt` to install all of the dependences.  You can then test to make sure everything installed by starting a python session and then running:
```
import dapt
dapt.__version__
```

## Project structure
```
.
├── dapt                 			# The folder where the DAPT library is housed
├── docs             				# Documentation for project
└── examples          				# Python scripts showing how to use the program
```

## Documentation
Documentation is done using [Sphinx](http://www.sphinx-doc.org/en/master/).  This allows for easy automatic documentation.  The *docs* folder holds all of the resources to document the code.  Here is a good tutorial on Sphinx if your not familiar <https://medium.com/@eikonomega/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365>.  Google docstrings are used for inline commenting inside each file.

## Unit tests
Unit tests are ran using [Pytest](pytest.org).  You can install it by running `pip install -U pytest`.  The tests are located in the `tests` folder inside of the `DAPT` module.  The tests can be run by opening a python session and then running:
```
import dapt
dapt.test()
```

or by running `pytest` in the main project directory.

For more information on the tests go to the [`tests`](dapt/tests) folder.