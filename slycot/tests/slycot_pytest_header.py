"""
slycot_pytest_header.py

pytest module to display slycot information in the test report header
"""

from slycot import __version__, __path__, _wrapper


def pytest_report_header(config):
    """Add slycot environment info to pytest report header"""
    return ['Slycot version: {}'.format(__version__),
            'Slycot path: {}'.format(__path__),
            'Slycot wrapper library: {}'.format(_wrapper.__file__)]
