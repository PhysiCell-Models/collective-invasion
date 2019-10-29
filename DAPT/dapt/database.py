"""
Database
========

The ``Database`` class is the basic interface for adding parameter set hosting services.  The idea is that the core methods (``get_table``, ``get_keys``, ``update_row`` and ``update_cell``) stay the same so that the inner workings can use multiple sources to access the parameter sets.  In general, there shouldn't be any more arguments added.  The exception to this is ``__init__``.  It may be necessary or desirable to add additional, optional, arguments.  This should be done by overloading the method or providing a default option to the argument.  

When prepareing for a parameter sweep the collection of parameter sets can be thought of as a database where the name of each paramate is a column name, each row is a paramater set and the value at the i-th column and j-th row is the value.  This is the approach of DAPT.
"""

class Database:
	"""
		An interface for accessing and setting parameter set data.  
	"""
		
	def __init__(self):

		pass

	def get_table(self):
		"""
	        Get the table from the database.

	        Returns:
	            An array with each element being a dictionary of the key-value pairs for the row in the database.
	    """

		pass

	def get_keys(self):
		"""
	        Get the keys of the parameter set
			
	        Returns:
	            Array of strings with each element being a key (order is preserved if possible)
	    """

		pass

	def update_row(self, row_id, values):
		"""
            Get the keys of the parameter set

            Args:
                row_id (int): the row id to replace
                values (OrderedDict): the key-value pairs that should be inserted
            
            Returns:
                A boolean that is True if successfully inserted and False otherwise.
        """

		pass

	def update_cell(self, row_id, key, value):
		"""
	        Get the keys of the parameter set

	        Args:
	            row_id (int): the row id to replace
	            key (str): the key of the value to replace
	            value (str): the value to insert into the cell
			
	        Returns:
	            A boolean that is True if successfully inserted and False otherwise.
	    """

		pass
