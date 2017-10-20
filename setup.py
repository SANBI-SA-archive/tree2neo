from setuptools import setup, find_packages

setup(
    name='tree2neo',
    version='0.0.3',
    description='Parses FastTree file and load it to the graph database.',
    keywords='neo4j,and fasttree',
    py_modules=['tree2neo'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click',
        'bioservices',
        'py2neo',
        'tqdm'
    ],
    entry_points={
        'console_scripts': ['tree2neo=tree2neo.cli:cli']
    },
)
