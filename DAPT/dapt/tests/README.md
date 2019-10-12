# DAPT Testing

Testing for DAP is done using [pytest](pytest.org) and can be installed by entering `pip install -U pytest` in the terminal. 

To run the tests open a python session in the terminal: `python3`.  Then type:
```
	import dapt
	dapt.test()
```

Additionally, you can run `pytest` in the root directory of the project.

You should now see all the tests run successfully.

## Command line arguments
| arg | options | usage |
| ---- |:----:| ----:|
| --color | [yes, no] | Turn terminal coloring on or off |


## ToDo
- `Delimited_file`
    - Add tests that make class fail like trying to write to invalid row