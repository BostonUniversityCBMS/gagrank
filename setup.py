try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

config = {
	'description': 'GAGrank',
	'author': 'J.D. Hogan',
	'url': 'URL to get it at.',
	'download_url': 'Where to download it.',
	'author_email': 'jdhogan@bu.edu',
	'version': '0.1',
	'install_requires': ['nose'],
	'packages': ['gagrank'],
	'scripts': [],
	'name': 'gagrank'
}

setup(**config)
