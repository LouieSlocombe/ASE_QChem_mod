from setuptools import setup, find_packages

setup(
    name='ASE_QChem_mod',
    version='0.1.0',
    author='Louie Slocombe',
    author_email='louies@hotmail.co.uk',
    description='Modified QChem calculator file.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/LouieSlocombe/ASE_QChem_mod',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy',
        'ase',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
        ],
    },
)
