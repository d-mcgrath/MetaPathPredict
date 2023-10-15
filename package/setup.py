from setuptools import Extension, setup, find_packages
import os

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: BSD License",
    "Operating System :: Linux, MacOS, Windows",
    "Programming Language :: Python :: 3.10.6+"
]

setup(
    name="metapathpredict", 
    description="Tool for predicting the presence or absence of KEGG modules in bacterial genomes",
    author="D. Geller-McGrath, K.M. Konwar, V.P. Edgcomb, M. Pachiadaki, J.W. Roddy, T.J. Wheeler, J.E. McDermott",
    author_email="dgellermcgrath@gmail.com, kishori82@gmail.com",
    package_dir={"": "src"},
    packages=["metapathpredict"],
    package_data={"metapathpredict": ["data/*.*"]},
    install_requires=[
      "scikit-learn>=1.1.3",
      "tensorflow>=2.10.0",
      "numpy>=1.23.4",
      "pandas>=1.5.2",
      "keras>=2.10.0",
      "torchvision>=0.15.2",
      "torch>=2.0.1",
    ],
    entry_points={
        "console_scripts": [
            "MetaPathTrain = metapathpredict.MetaPathPredict:Models.train", 
            "MetaPathPredict = metapathpredict.MetaPathPredict:Models.predict", 
            "MetaPathModules = metapathpredict.MetaPathPredict:Models.show_available_modules",
            "DownloadModels = metapathpredict.download_models:Download.download_models",
            "PredictFromTable = metapathpredict.MetaPathPredict:Models.predict_from_feature_table",
            "PredictFromTableFs = metapathpredict.MetaPathPredict:Models.predict_from_feature_table_fs_models"
        ]
    },
    classifiers=CLASSIFIERS,
    include_package_data=True,
    #ext_modules=cythonize("src/metapathpredict/cpp_mods.pyx")
 )
