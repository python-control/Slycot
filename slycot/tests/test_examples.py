"""

test_examples.py


"""

from inspect import getmembers, isfunction
import pytest

from slycot import examples

examplefunctions = [fun for (fname, fun) in getmembers(examples)
                    if isfunction(fun) and "_example" in fname]


@pytest.mark.parametrize('examplefun', examplefunctions)
def test_example(examplefun, capsys, recwarn):
    """
    Test the examples.

    Test that all the examples work, produce some (unchecked) output but no
    exceptions or warnings.
    """
    examplefun()
    captured = capsys.readouterr()
    assert len(captured.out) > 0
    assert not captured.err
    assert not recwarn
