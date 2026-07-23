import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csc_matrix, csr_matrix

# First 100 genes from ../../src/gene_list/chemistry/xenium_v1_human.rs
GENES = [
    ("ENSG00000262902", "MTCO1P40"),
    ("ENSG00000225493", "LINC01107"),
    ("ENSG00000231934", "FAM242A"),
    ("ENSG00000243104", "MTND4LP14"),
    ("ENSG00000241838", "AP000529.1"),
    ("ENSG00000279048", "AC080080.1"),
    ("ENSG00000196961", "AP2A1"),
    ("ENSG00000276077", "CU633904.1"),
    ("ENSG00000273356", "LINC02019"),
    ("ENSG00000287897", "AP000465.1"),
    ("ENSG00000278905", "AC106818.1"),
    ("ENSG00000223897", "AC069157.1"),
    ("ENSG00000285748", "AC090337.2"),
    ("ENSG00000198081", "ZBTB14"),
    ("ENSG00000119899", "SLC17A5"),
    ("ENSG00000140675", "SLC5A2"),
    ("ENSG00000278601", "AL158163.2"),
    ("ENSG00000228223", "HCG11"),
    ("ENSG00000148158", "SNX30"),
    ("ENSG00000185043", "CIB1"),
    ("ENSG00000134323", "MYCN"),
    ("ENSG00000213757", "AC020898.1"),
    ("ENSG00000214857", "SEM1P1"),
    ("ENSG00000099991", "CABIN1"),
    ("ENSG00000274177", "AC004528.2"),
    ("ENSG00000143536", "CRNN"),
    ("ENSG00000102898", "NUTF2"),
    ("ENSG00000223598", "AL138733.1"),
    ("ENSG00000146360", "GPR6"),
    ("ENSG00000175309", "PHYKPL"),
    ("ENSG00000272543", "MIR4787"),
    ("ENSG00000175573", "C11orf68"),
    ("ENSG00000275207", "MIR6740"),
    ("ENSG00000275457", "AL117332.1"),
    ("ENSG00000226400", "AC010385.1"),
    ("ENSG00000247228", "AC009060.1"),
    ("ENSG00000129965", "INS-IGF2"),
    ("ENSG00000087470", "DNM1L"),
    ("ENSG00000254261", "AL451137.2"),
    ("ENSG00000253887", "AC022784.3"),
    ("ENSG00000201487", "SNORD45B"),
    ("ENSG00000252548", "RNU7-149P"),
    ("ENSG00000242524", "OR2U2P"),
    ("ENSG00000199395", "RNA5SP93"),
    ("ENSG00000105642", "KCNN1"),
    ("ENSG00000263435", "AC024610.1"),
    ("ENSG00000199673", "SNORD16"),
    ("ENSG00000255337", "AP001830.1"),
    ("ENSG00000169016", "E2F6"),
    ("ENSG00000270762", "FXYD6P1"),
    ("ENSG00000226138", "AC004801.1"),
    ("ENSG00000144488", "ESPNL"),
    ("ENSG00000174514", "MFSD4A"),
    ("ENSG00000133142", "TCEAL4"),
    ("ENSG00000243505", "RN7SL240P"),
    ("ENSG00000163798", "SLC4A1AP"),
    ("ENSG00000237664", "LINC00316"),
    ("ENSG00000255088", "AC013828.1"),
    ("ENSG00000252481", "SCARNA13"),
    ("ENSG00000068305", "MEF2A"),
    ("ENSG00000174353", "STAG3L3"),
    ("ENSG00000254037", "AC021733.3"),
    ("ENSG00000226140", "LINC02654"),
    ("ENSG00000235455", "IQCF5-AS1"),
    ("ENSG00000261709", "AC051619.9"),
    ("ENSG00000274308", "AC244093.1"),
    ("ENSG00000180475", "OR10Q1"),
    ("ENSG00000285803", "AL442003.1"),
    ("ENSG00000228340", "MIR646HG"),
    ("ENSG00000276874", "AC011718.1"),
    ("ENSG00000101052", "IFT52"),
    ("ENSG00000215644", "GCGR"),
    ("ENSG00000231152", "MTND2P15"),
    ("ENSG00000253668", "AC103778.1"),
    ("ENSG00000177427", "MIEF2"),
    ("ENSG00000271974", "RDM1P4"),
    ("ENSG00000271269", "AL353778.1"),
    ("ENSG00000223335", "RNU6-603P"),
    ("ENSG00000237186", "AC092418.1"),
    ("ENSG00000234219", "CDCA4P4"),
    ("ENSG00000100387", "RBX1"),
    ("ENSG00000286306", "AL139811.3"),
    ("ENSG00000284221", "AC099654.10"),
    ("ENSG00000135974", "C2orf49"),
    ("ENSG00000277397", "AL606760.4"),
    ("ENSG00000229671", "LINC01150"),
    ("ENSG00000252050", "AL606748.2"),
    ("ENSG00000287610", "AC009403.3"),
    ("ENSG00000206970", "RNU6-474P"),
    ("ENSG00000199165", "MIRLET7A1"),
    ("ENSG00000254633", "AL512590.2"),
    ("ENSG00000223866", "AC002486.1"),
    ("ENSG00000186310", "NAP1L3"),
    ("ENSG00000215912", "TTC34"),
    ("ENSG00000273145", "BX537318.1"),
    ("ENSG00000206672", "Y_RNA"),
    ("ENSG00000249610", "AC010255.2"),
    ("ENSG00000257025", "AC023595.1"),
    ("ENSG00000280433", "FP565260.6"),
    ("ENSG00000238007", "AC024619.1"),
]


def main():
    # Adapted from https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html
    ensembl_ids, gene_names = zip(*GENES)
    n_cells = 10

    rng = np.random.default_rng()
    counts = rng.integers(0, 10, size=(n_cells, len(GENES)))
    # Set the first element = 10 so we can assert against it in tests
    counts[0] = 10

    csr_adata = ad.AnnData(csr_matrix(counts, dtype=np.float32))
    csc_adata = ad.AnnData(csc_matrix(counts, dtype=np.float32))

    for name, adata in [("csr_adata", csr_adata), ("csc_adata", csc_adata)]:
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.obs["annotation"] = (["group1"] * 5) + (["group2"] * 5)

        adata.var_names = list(gene_names)
        adata.var["ensembl_id"] = list(ensembl_ids)
        adata.var["gene_name"] = list(gene_names)

        adata.write_h5ad(f"../{name}.h5ad")

    filename = "../WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5"
    tenx_adata = sc.read_10x_h5(filename)
    sc.write(f"{filename}ad", sc.pp.subsample(tenx_adata, 0.01, copy=True))


if __name__ == "__main__":
    main()
