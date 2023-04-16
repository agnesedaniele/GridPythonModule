from setuptools import setup

setup(
    name='GridPythonModule',
    version='0.1.0',    
    description='A Sage compatible Python module to manipulate and simplify grid diagrams.',
    url='https://github.com/agnesedaniele/GridPythonModule',
    author='Agnese Barbensi and Daniele Celoria',
    author_email='dceloria.maths@gmail.com',
    license='GNU general public',
    packages=['GridPythonModule'],
    install_requires=['simpy',
                      'random2',
                      'matplotlib'
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
)

