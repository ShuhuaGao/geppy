# coding=utf-8
import setuptools

try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except:
    long_description = '''*geppy* is developed for fast prototyping of **gene expression programming** (GEP) in Python, which is
    built on top of the DEAP package. Please the [GitHub repository](https://github.com/ShuhuaGao/geppy) for details.'''

setuptools.setup(
    name="geppy",
    version="0.1.3",
    license='LGPL-3.0 License',
    author="Shuhua Gao",
    author_email="nus.gao.shuhua@gmail.com",
    description="A package for gene expression programming in Python.",
    long_description=long_description,
    url="https://github.com/ShuhuaGao/geppy",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
	long_description_content_type='text/markdown',
    keywords=['evolutionary computation', 'gene expression programming',
              'computational intelligence', 'genetic programming'],
    packages=setuptools.find_packages(),
    install_requires=['deap'],
    download_url='https://github.com/ShuhuaGao/geppy/archive/v0.1.3.tar.gz'
)

