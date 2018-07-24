# coding=utf-8
import setuptools

try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except:
    long_description = '''geppy is developed for fast prototyping of gene expression programming in Python, and it is 
    built on top of the DEAP package.'''

setuptools.setup(
    name="geppy",
    version="0.1.0a",
    author="Shuhua Gao",
    author_email="nus.gao.shuhua@gmail.com",
    description="A package for gene expression programming in Python.",
    long_description=long_description,
    url="https://github.com/ShuhuaGao/geppy",
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: LGPL-v3.0 License",
        "Operating System :: OS Independent",
    ),
    keywords='gene expression programming',
    packages=['geppy'],
    install_requires=['deap']
)

