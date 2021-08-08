from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='prolintpy',
    version='0.9.2',
    description='Automated analyis and visualization of lipid-protein interactions.',
    url="https://github.com/ProLint/prolintpy",
    author='Besian I. Sejdiu',
    author_email="besian.sejdiu@stjude.org",
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",
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
	'nglview', 
	'matplotlib'],
)

