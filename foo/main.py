import h5py
import scanpy as sc


def main():
    adata = sc.read_10x_h5(
        "SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5"
    )
    adata = sc.pp.subsample(adata, 0.01, copy=True)
    sc.write(
        "SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix_small.h5ad",
        adata,
    )


if __name__ == "__main__":
    main()
