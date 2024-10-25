#import pyxet
import importlib
import shutil
from importlib import resources
from huggingface_hub import hf_hub_download


class Download:
    """Functions to download MetaPathPredict's machine learning models"""
    
    @classmethod
    def download_models(cls):
      """Downloads MetaPathPredict's models.

      Returns:
        None

      """
      print("Downloading MetaPathPredict models...")
      module_dir = resources.files('metapathpredict')
      data_dir = module_dir.joinpath("data/")
      # model_0_dl_path = "xet://dgellermcgrath/MetaPathPredict/main/package/src/metapathpredict/data/model_0.keras"
      # model_1_dl_path = "xet://dgellermcgrath/MetaPathPredict/main/package/src/metapathpredict/data/model_1.keras"
      model_0_install_path = module_dir.joinpath("data/MetaPathPredict_model_0.keras")
      model_1_install_path = module_dir.joinpath("data/MetaPathPredict_model_1.keras")
      
      model_0_renamed_dir_path = module_dir.joinpath("data/model_0.keras_directory")
      model_1_renamed_dir_path = module_dir.joinpath("data/model_1.keras_directory")
      
      model_0_initial_path = module_dir.joinpath("data/model_0.keras_directory/MetaPathPredict_model_0.keras")
      model_1_initial_path = module_dir.joinpath("data/model_1.keras_directory/MetaPathPredict_model_1.keras")
      
      model_0_final_path = module_dir.joinpath("data/model_0.keras")
      model_1_final_path = module_dir.joinpath("data/model_1.keras")

      hf_hub_download(repo_id="dgellermcgrath/MetaPathPredict", filename="MetaPathPredict_model_0.keras", local_dir=model_0_install_path, force_download=True)
      hf_hub_download(repo_id="dgellermcgrath/MetaPathPredict", filename="MetaPathPredict_model_1.keras", local_dir=model_1_install_path, force_download=True)
      
      # rename the model directories downloaded from HuggingFace
      shutil.move(model_0_install_path, model_0_renamed_dir_path)
      shutil.move(model_1_install_path, model_1_renamed_dir_path)
      
      # move the models out of their directories and rename them
      shutil.move(model_0_initial_path, model_0_final_path)
      shutil.move(model_1_initial_path, model_1_final_path)
      
      # remove the directories downloaded from HuggingFace
      shutil.rmtree(model_0_renamed_dir_path)
      shutil.rmtree(model_1_renamed_dir_path)

      # fs = pyxet.XetFS()  # fsspec filesystem
      # fs.get(model_0_dl_path, str(model_0_install_path))
      # fs.get(model_1_dl_path, str(model_1_install_path))
      print("All done. Use MetaPathPredict -h to see how to make predictions.")
