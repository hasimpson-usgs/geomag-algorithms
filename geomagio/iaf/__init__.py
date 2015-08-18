"""IO Module for IAF Format

Based on documentation at:
  http://www.intermagnet.org/data-donnee/formats/iaf-eng.php
"""

from IAFFactory import IAFFactory
from StreamIAFFactory import StreamIAFFactory
from IAFParser import IAFParser
from IAFWriter import IAFWriter


__all__ = [
    'IAFFactory',
    'StreamIAFFactory',
    'IAFParser',
    'IAFWriter'
]
