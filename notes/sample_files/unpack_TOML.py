import toml
from collections import defaultdict

metadata_keys = ["is_MC", "is_UL", "is_SMS"]

def _format_fileset(toml_file, dataset_names: list[tuple[str, int]]): # [(dataset_name, num_files), ...]
    """
    A function that returns a dict of a particular format that the Runner class in coffea takes.
    
    Many examples of this fileset format can be found through coffea's documentation, but I prefer TOML's
    
    """

    
    open_toml = toml.load(toml_file)
    fileset = defaultdict(dict)
    
    for dataset_name in dataset_names[0]:
        
        files_dict = {}
        for file in open_toml[dataset_name]["files"]:
            files_dict[file] = open_toml[dataset_name]["branch"]
        fileset[dataset_name]["files"] = files_dict
        
        metadata_dict = {}
        for key in metadata_keys:
            metadata_dict[key] = open_toml[dataset_name][key]
            fileset[dataset_name]["metadata"] = metadata_dict
    
        
    return dict(fileset)  # Convert back to regular dict

def _format_reduced_fileset(toml_file, dataset_names: list[tuple[str, int]])


xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/SMS-TChiWZ_ZToLL_mZMin-0p1_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/80000/0B041116-DADF-CF45-8C16-DF2920CC756D.root .
    
def _format_xroot_query(toml_file):

    