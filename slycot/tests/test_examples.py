"""

test_examples.py


"""

from inspect import getmembers, isfunction
import unittest, io
from parameterized import parameterized
from contextlib import redirect_stdout, redirect_stderr

from slycot import examples

examplefunctions = [(fname[:-len("_example")], fun)
                    for (fname, fun) in getmembers(examples)
                    if isfunction(fun) and "_example" in fname]

class TestExamples(unittest.TestCase):

    @parameterized.expand(examplefunctions)
    def test_example(self, name, examplefun):
        #def test_example(examplefun, capsys, recwarn):
        """
        Test the examples.
        
        Test that all the examples work, produce some (unchecked) output but no
        exceptions or warnings.
        """
        out, err = io.StringIO(), io.StringIO()
        with redirect_stderr(err), redirect_stdout(out):
            examplefun()
        self.assertTrue(len(out.getvalue()) > 0)
        self.assertFalse(err.getvalue())
        #assert not recwarn

if __name__ == "__main__":
    unittest.main()
