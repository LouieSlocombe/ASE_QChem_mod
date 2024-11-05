from setuptools import setup, find_packages

setup(
    name='nqetools',
    version='0.0.0',
    author='Louie Slocombe',
    author_email='louies@hotmail.co.uk',
    description='A centralised set of tools for doing nuclear quantum calculations.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/LouieSlocombe/nqetools',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
    install_requires=[
        'numpy',
        'matplotlib',
        'ase',
        'scipy',
        'ipi',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
        ],
    },
)
