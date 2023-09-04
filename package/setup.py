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
    description="Tool for predicting and metabolic pathway on metagenomes",
    author="ABCD",
    author_email="davids  and kishoris email",
    package_dir={"": "src"},
    packages=["metapathpredict"],
    install_requires=[
    ],
    entry_points={
        "console_scripts": [
            # trains and tests the CPR models
            "MetaPathPredictTrain = metapathpredict.cmdline_models:Models.train", 
            "MetaPathPredictPredict = metapathpredict.cmdline_models:Models.predict", 
        ]
    },
    classifiers=CLASSIFIERS,
    include_package_data=True,
    #ext_modules=cythonize("src/metapathpredict/cpp_mods.pyx")
 )
