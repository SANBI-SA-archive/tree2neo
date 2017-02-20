from setuptools import setup

setup(
    name='tree2neo',
    version='0.0.1',
    description='Parses FastTree file and load it to the graph database.',
    keywords='neo4j,and fasttree',
    py_modules=['tree2neo'],
    install_requires=[
        'click',
        'bioservices',
        'py2neo'
    ],
    entry_points={
        'console_scripts': ['tree2neo=tree2neo.cli:cli']
    },
)
