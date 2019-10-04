from setuptools import setup

setup(name='graph_read_simulator',
      version='0.0.1',
      description='Graph Read Simulator',
      url='http://github.com/ivargr/graph_read_simulator',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=['graph_read_simulator'],
      zip_safe=False,
      install_requires=['numpy', 'pyvg==1.1.2', 'offsetbasedgraph', 'simple_read_mutator', 'pyfaidx',
                        'graph_peak_caller==1.2.3'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
             'console_scripts': ['graph_read_simulator=graph_read_simulator.command_line_interface:main'],
      })

"""
To update package:
#Update version number manually in this file

sudo python3 setup.py sdist
sudo python3 setup.py bdist_wheel
twine upload dist/graph_read_simulator_X.tar.gz
"""
