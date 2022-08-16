from setuptools import setup, find_packages, Extension


setup(
    name='myspk',
    version='0.1.0',
    author='Ori and Itamar',
    author_email='orisofer@mail.tau.ac.il',
    description='Spectral Clustering',
    install_requires=['invoke'],
    packages=find_packages(),
    license='GPL-2',
    classifiers=[
        'Development status:: 4 - beta',
        'Liscense:: OSI',
    ],
    ext_modules=[
        Extension(
            'myspk',
            ['spkmeans.c']
        )
    ]
)