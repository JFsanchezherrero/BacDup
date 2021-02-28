from pkg_resources import get_distribution
try:
    __version__ = get_distribution('BacDup').version
except:
    __version__ = 'local'

__all__ = [
	'modules',
	'scripts',
	'config'
]

from BacDup import *

