use hdf5_metno::{Dataset, File, Group, types::VarLenUnicode};
use ndarray::Array1;
use serde::Serialize;
use sprs::CsMatBase;

pub type UmiCounts = CsMatBase<u32, usize, Vec<usize>, Vec<usize>, Vec<u32>>;

pub fn read_umi_counts_from_h5ad(file: &File) -> Result<UmiCounts, Error> {
    const X: &str = "X";
    // The actual counts (stored by scanpy as the highly-descriptive name "X") are usually stored as in compressed-sparse matrix format in a group, but they might also be in a dataset
    match file.group(X) {
        Ok(g) => read_x_group(&g),
        Err(_) => read_x_dataset(&file.dataset(X)?),
    }
}

fn read_x_group(x: &Group) -> Result<UmiCounts, Error> {
    // First, read the data and fallibly cast every element to a u32 (which indicates untransformed counts). If it fails, we know we have some kind of transformed data
    let data: Array1<f32> = x.dataset("data").map(|ds| ds.read_1d()).flatten()?;
    let data: Vec<_> = data.into_iter().map(f32_to_u32).collect::<Result<_, _>>()?;

    let nondata_parts: Vec<Vec<usize>> = ["indptr", "indices"]
        .into_iter()
        .flat_map(|key| x.dataset(key).map(|ds| ds.read_raw()))
        .collect::<Result<_, _>>()?;
    let [indptr, indices] = nondata_parts.try_into().unwrap();

    // It's very nice that scanpy decides to store the shape as an attribute rather than a dataset, the way 10x Genomics does it :) I don't understand this decision
    let shape = x.attr("shape").map(|sh| sh.read_1d()).flatten()?;
    let shape = (shape[0], shape[1]);

    let encoding_type: VarLenUnicode =
        x.attr("encoding-type").map(|a| a.read_scalar()).flatten()?;

    let counts = match encoding_type.as_str() {
        "csr_matrix" => UmiCounts::new(shape, indptr, indices, data),
        "csc_matrix" => UmiCounts::new_csc(shape, indptr, indices, data),
        _ => {
            return Err(Error::UnknownEncodingType {
                encoding_type: encoding_type.to_string(),
            });
        }
    };

    Ok(counts)
}

fn read_x_dataset(x: &Dataset) -> Result<UmiCounts, Error> {
    todo!()
}

fn f32_to_u32(f: f32) -> Result<u32, Error> {
    let is_integer = f.round() == f;
    let is_nonnegative = f >= 0.0;

    if is_integer && is_nonnegative {
        Ok(f as u32)
    } else {
        Err(Error::TransformedCounts)
    }
}

#[derive(Debug, thiserror::Error, Serialize)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum Error {
    #[error("non-raw counts found")]
    TransformedCounts,
    #[error("{reason}")]
    Hdf5 { reason: String },
    #[error("unknown encoding-type")]
    UnknownEncodingType { encoding_type: String },
}

impl From<hdf5_metno::Error> for Error {
    fn from(err: hdf5_metno::Error) -> Self {
        Self::Hdf5 {
            reason: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use hdf5_metno::File;

    use crate::reference_dataset::umi_counts::read_umi_counts_from_h5ad;

    #[test]
    fn read_h5ad_files() {
        let files = ["csr_adata", "csc_adata", "WT_mouse_spinal_cord_P112_specimen_1_WT_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix"]
            .map(|fname| format!("test-data/{fname}.h5ad"))
            .map(|path| File::open(path).unwrap());

        for f in files {
            let counts = read_umi_counts_from_h5ad(&f).unwrap();
            assert_eq!(counts.data()[0], 10);
        }
    }
}
