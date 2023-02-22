from setuptools import setup

setup(name='graph_read_simulator',
      version='0.0.15',
      description='Graph Read Simulator',
      url='http://github.com/ivargr/graph_read_simulator',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=['graph_read_simulator'],
      zip_safe=False,
      install_requires=['numpy', 'simple_read_mutator', 'pyfaidx', 'obgraph>=0.0.33', "biopython"],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
             'console_scripts': ['graph_read_simulator=graph_read_simulator.command_line_interface:main'],
      })

"""
To update package:
#Update version number manually in this file

rm -rf dist
python3 setup.py sdist
python3 setup.py bdist_wheel
twine upload dist/graph_read_simulator-*.tar.gz

"""
