from setuptools import Extension, setup, find_packages
import os

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: BSD License",
    "Operating System :: Linux, MacOS",
    "Programming Language :: Python :: 3.6+"
]

setup(
    name="metapathpredict", 
    description="Tool for predicting the presence or absence of KEGG modules in bacterial genomes",
    author="D. Geller-McGrath, K.M. Konwar, V.P. Edgcomb, M. Pachiadaki, J.W. Roddy, T.J. Wheeler, J.E. McDermott",
    author_email="dgellermcgrath@gmail.com, kishori82@gmail.com",
    package_dir={"": "src"},
    packages=["metapathpredict"],
    install_requires=[
    ],
    entry_points={
        "console_scripts": [
            # trains and tests the CPR models
            "MetaPathTrain = metapathpredict.cmdline_models:Models.train", 
            "MetaPathPredict = metapathpredict.cmdline_models:Models.predict", 
        ]
    },
    classifiers=CLASSIFIERS,
    include_package_data=True,
    #ext_modules=cythonize("src/metapathpredict/cpp_mods.pyx")
 )
