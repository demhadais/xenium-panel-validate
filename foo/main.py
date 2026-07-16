import h5py
import scanpy as sc


def main():
    adata = sc.read_h5ad("../crates/core/test-data/sample.h5ad")
    print(adata.)

    # adata.write(
    #     "../crates/core/test-data/SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5ad"
    # )

    # f = h5py.File(
    #     "../crates/core/test-data/SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix.h5ad"
    # )
    # print(f)
    # print(*(f.attrs.items()))
    # print(*(f.items()))


if __name__ == "__main__":
    main()
