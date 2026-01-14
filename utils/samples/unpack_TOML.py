import toml
from collections import defaultdict

metadata_keys = ["is_MC", "is_UL", "is_SMS"]

def _format_fileset(
    toml_file, 
    dataset: str,
    num_files: int,
    redirector = "CERN"
): 
    """
    A function that returns a dict of a particular format that the Runner class in coffea takes.
    
    Many examples of this fileset format can be found through coffea's documentation, but I prefer TOML's

    Example of call: unpack_TOML._format_fileset(toml_path, ["JetMET2022C", -1])

    where -1 means use all available files
    
    """

    
    
    open_toml = toml.load(toml_file)
    
    fileset = defaultdict(dict)

    red_key = redirector.lower()

    formatted_file_list = []

    for file in open_toml[dataset]["files"]:
        formatted_file_list.append(open_toml["redirectors"][red_key] + file)

    metadata_dict = {}
    for key in metadata_keys:
        metadata_dict[key] = open_toml[dataset][key]
        fileset[dataset]["metadata"] = metadata_dict

    return {
        dataset: {
            "files": formatted_file_list,
            "treename": open_toml[dataset]["treename"],
            "metadata": metadata_dict,
        }
    }
    

def _format_reduced_fileset(toml_file, dataset_names: list[tuple[str, int]]):
    return


"""
xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17NanoAODv9/SMS-TChiWZ_ZToLL_mZMin-0p1_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/80000/0B041116-DADF-CF45-8C16-DF2920CC756D.root .
"""

def _format_xroot_query(toml_file):
    return
    