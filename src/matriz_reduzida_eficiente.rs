use nalgebra::DMatrix;
use rayon::prelude::*;
#[derive(Debug)]
pub struct CsrMatrix {
    pub n: usize,
    pub row_ptr: Vec<usize>,
    pub col_ind: Vec<usize>,
    pub values: Vec<f64>,
}

pub fn dense_to_csr(mat: &DMatrix<f64>, tol: f64) -> CsrMatrix {
    let nrows = mat.nrows();
    let ncols = mat.ncols();

    // Cada thread acumula contribuições por linha
    let per_thread: Vec<Vec<(usize, f64)>> = (0..ncols)
        .into_par_iter()
        .map(|j| {
            let mut local: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nrows];

            // Acesso contíguo à coluna
            let col = mat.column(j);

            for (i, &v) in col.iter().enumerate() {
                if v.abs() > tol {
                    local[i].push((j, v));
                }
            }

            local
        })
        .reduce(
            || vec![Vec::new(); nrows],
            |mut acc, mut local| {
                for i in 0..nrows {
                    acc[i].extend(local[i].drain(..));
                }
                acc
            },
        );

    // Montagem CSR
    let mut row_ptr = Vec::with_capacity(nrows + 1);
    let mut col_ind = Vec::new();
    let mut values = Vec::new();

    row_ptr.push(0);

    for row in per_thread {
        for (j, v) in row {
            col_ind.push(j);
            values.push(v);
        }
        row_ptr.push(col_ind.len());
    }

    CsrMatrix {
        n: nrows,
        row_ptr,
        col_ind,
        values,
    }
}

pub fn reduzir_csr(k: &CsrMatrix, removidos: &[usize]) -> CsrMatrix {
    let n = k.n;

    // 1. máscara
    let mut livre = vec![true; n];
    for &r in removidos {
        livre[r] = false;
    }

    // 2. mapeamento global -> reduzido
    let mut map = vec![usize::MAX; n];
    let mut count = 0;
    for i in 0..n {
        if livre[i] {
            map[i] = count;
            count += 1;
        }
    }
    let m = count;

    // 3. Processamento paralelo por linhas
    let rows: Vec<(Vec<usize>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .filter(|&i| livre[i])
        .map(|i| {
            let mut cols = Vec::new();
            let mut vals = Vec::new();

            let start = k.row_ptr[i];
            let end = k.row_ptr[i + 1];

            for idx in start..end {
                let j = k.col_ind[idx];
                if livre[j] {
                    cols.push(map[j]);
                    vals.push(k.values[idx]);
                }
            }

            (cols, vals)
        })
        .collect();

    // 4. Montagem CSR final
    let mut row_ptr = Vec::with_capacity(m + 1);
    let mut col_ind = Vec::new();
    let mut values = Vec::new();

    row_ptr.push(0);

    for (cols, vals) in rows {
        col_ind.extend(cols);
        values.extend(vals);
        row_ptr.push(col_ind.len());
    }

    CsrMatrix {
        n: m,
        row_ptr,
        col_ind,
        values,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    fn approx_eq(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| (x - y).abs() < tol)
    }

    #[test]
    fn test_small_matrix() {
        let mat = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0]);

        let csr = dense_to_csr(&mat, 1e-12);

        assert_eq!(csr.n, 3);
        assert_eq!(csr.row_ptr, vec![0, 2, 3, 4]);
        assert_eq!(csr.col_ind, vec![0, 2, 2, 0]);
        assert!(approx_eq(&csr.values, &[1.0, 2.0, 3.0, 4.0], 1e-12));
    }

    #[test]
    fn test_zero_matrix() {
        let mat = DMatrix::<f64>::zeros(4, 4);

        let csr = dense_to_csr(&mat, 1e-12);

        assert_eq!(csr.n, 4);
        assert_eq!(csr.row_ptr, vec![0, 0, 0, 0, 0]);
        assert!(csr.col_ind.is_empty());
        assert!(csr.values.is_empty());
    }

    #[test]
    fn test_dense_matrix() {
        let mat = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);

        let csr = dense_to_csr(&mat, 1e-12);

        assert_eq!(csr.row_ptr, vec![0, 2, 4]);
        assert_eq!(csr.col_ind, vec![0, 1, 0, 1]);
        assert!(approx_eq(&csr.values, &[1.0, 2.0, 3.0, 4.0], 1e-12));
    }

    #[test]
    fn test_tolerance_filtering() {
        let mat = DMatrix::from_row_slice(2, 3, &[1e-10, 1.0, 0.0, 0.0, 1e-11, 2.0]);

        let csr = dense_to_csr(&mat, 1e-9);

        assert_eq!(csr.row_ptr, vec![0, 1, 2]);
        assert_eq!(csr.col_ind, vec![1, 2]);
        assert!(approx_eq(&csr.values, &[1.0, 2.0], 1e-12));
    }

    #[test]
    fn test_rectangular_matrix() {
        let mat = DMatrix::from_row_slice(2, 4, &[0.0, 5.0, 0.0, 6.0, 7.0, 0.0, 0.0, 8.0]);

        let csr = dense_to_csr(&mat, 1e-12);

        assert_eq!(csr.row_ptr, vec![0, 2, 4]);
        assert_eq!(csr.col_ind, vec![1, 3, 0, 3]);
        assert!(approx_eq(&csr.values, &[5.0, 6.0, 7.0, 8.0], 1e-12));
    }

    #[test]
    fn test_row_ptr_monotonicity() {
        let mat = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]);

        println!("matriz = {}",mat);

        let csr = dense_to_csr(&mat, 1e-12);

        // row_ptr deve ser não-decrescente
        for i in 0..csr.row_ptr.len() - 1 {
            assert!(csr.row_ptr[i] <= csr.row_ptr[i + 1]);
        }

        // último valor = número total de nnz
        assert_eq!(csr.row_ptr.last().copied().unwrap(), csr.values.len());
    }
}
