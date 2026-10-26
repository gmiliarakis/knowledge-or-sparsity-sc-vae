# Which dataset, gene set collection and benchmark a script works on (imported by the scripts, not run on its own).
# SC_DATA selects the dataset (kang, the default, or a file data/processed/<name>.h5ad such as cytokine_ifn_beta), SC_SETS the collection
# (hallmark, the default, or reactome) and SC_TAG a derived data benchmark (drop250). Files that depend on the data carry the
# dataset and SC_TAG suffixes. Results and tables go to results/<data>[_<tag>]/<collection>/ and reports/tables/<data>[_<tag>]/<collection>/.
import os

collection = os.environ.get("SC_SETS", "hallmark")
# the interferon sets used as the known answer: type I (IFN-beta is a type I interferon) and type II
targets = {"hallmark": ("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
           "reactome": ("REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "REACTOME_INTERFERON_GAMMA_SIGNALING")}
assert collection in targets, f"SC_SETS must be one of {list(targets)}"
alpha_name, gamma_name = targets[collection]
if os.environ.get("SC_TARGET") == "il2":  # IL-2 contrast: the known answer is the IL-2 set, not the type I interferon set (only for 08_evaluate.py)
    alpha_name = {"hallmark": "HALLMARK_IL2_STAT5_SIGNALING", "reactome": "REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING"}[collection]

data = os.environ.get("SC_DATA", "kang")
tag = "_" + os.environ["SC_TAG"] if os.environ.get("SC_TAG") else ""
file_tag = ("" if data == "kang" else "_" + data) + tag  # suffix of the mask files
run = f"{data}{tag}/{collection}"  # folder of the results and tables of this benchmark, for example kang/hallmark or cytokine_il_2/reactome
main = run == "kang/hallmark"  # the main experiment
