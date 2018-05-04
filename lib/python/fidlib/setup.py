from setuptools import setup

setup(
    name='fidlib',
    version='0.1',
    py_modules=['fidlib'],
    install_requires=[
        'pyfasta>=0.5.2',
        'pandas>=0.18',
        'bx-python>=0.7.1'
    ],
    author='Ian Fiddes',
    license='Apache 2.0',
)
