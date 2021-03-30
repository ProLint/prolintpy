from setuptools import setup, find_packages

setup(
    name='prolint',
    version='0.9.0',
    description='Lipid-protein analysis',
    author='Besian I. Sejdiu',
    license='MIT',
    python_requires='>=3.6',
    packages=find_packages(),
    install_requires=[
	'numpy', 
	'pandas==0.24.0', 
	'scipy', 
	'mdtraj', 
	'pyyaml', 
	'colorcet', 
	'bokeh==1.4.0', 
	'networkx', 
	'jupyterlab', 
	'nglview==2.7.7', 
	'matplotlib'],
)

