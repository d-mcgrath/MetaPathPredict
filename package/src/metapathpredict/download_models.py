#import pyxet
import importlib
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
      model_0_install_path = module_dir.joinpath("data/model_0.keras")
      model_1_install_path = module_dir.joinpath("data/model_1.keras")

      hf_hub_download(repo_id="dgellermcgrath/MetaPathPredict", filename="model_0.keras", cache_dir=model_0_install_path)
      hf_hub_download(repo_id="dgellermcgrath/MetaPathPredict", filename="model_1.keras", cache_dir=model_1_install_path)

      # fs = pyxet.XetFS()  # fsspec filesystem
      # fs.get(model_0_dl_path, str(model_0_install_path))
      # fs.get(model_1_dl_path, str(model_1_install_path))
      print("All done. Use MetaPathPredict -h to see how to make predictions.")
