use rayon::prelude::*;

pub struct CsrMatrix {
    pub n: usize,
    pub row_ptr: Vec<usize>,
    pub col_ind: Vec<usize>,
    pub values: Vec<f64>,
}


use rayon::prelude::*;

pub fn dense_to_csr(mat: &[Vec<f64>], tol: f64) -> CsrMatrix {
    let n = mat.len();

    // Processa cada linha em paralelo
    let rows: Vec<(Vec<usize>, Vec<f64>)> = mat
        .par_iter()
        .map(|row| {
            let mut cols = Vec::new();
            let mut vals = Vec::new();

            for (j, &v) in row.iter().enumerate() {
                if v.abs() > tol {
                    cols.push(j);
                    vals.push(v);
                }
            }

            (cols, vals)
        })
        .collect();

    // Montagem CSR final (sequencial para manter ordem)
    let mut row_ptr = Vec::with_capacity(n + 1);
    let mut col_ind = Vec::new();
    let mut values = Vec::new();

    row_ptr.push(0);

    for (cols, vals) in rows {
        col_ind.extend(cols);
        values.extend(vals);
        row_ptr.push(col_ind.len());
    }

    CsrMatrix {
        n,
        row_ptr,
        col_ind,
        values,
    }
}


pub fn reduzir_csr(
    k: &CsrMatrix,
    removidos: &[usize],
) -> CsrMatrix {
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

    fn exemplo_matriz() -> CsrMatrix {
        //[10  0  2  0  0]
        //[ 0 20  0  0  3]
        //[ 2  0 30  4  0]
        //[ 0  0  4 40  5]
        //[ 0  3  0  5 50]
        // Matriz 5x5 em CSR
        CsrMatrix {
            n: 5,
            row_ptr: vec![0, 2, 4, 7, 10, 13],
            col_ind: vec![
                0, 2,        // linha 0
                1, 4,        // linha 1
                0, 2, 3,     // linha 2
                2, 3, 4,     // linha 3
                1, 3, 4,     // linha 4
            ],
            values: vec![
                10.0, 2.0,
                20.0, 3.0,
                2.0, 30.0, 4.0,
                4.0, 40.0, 5.0,
                3.0, 5.0, 50.0,
            ],
        }
    }
    #[test]
    fn testa_reducao_csr() {
        let k = exemplo_matriz();

        // Remover DOFs 1 e 3
        let removidos = vec![1, 3];

        let kr = reduzir_csr(&k, &removidos);

        // Esperado: índices livres {0, 2, 4}
        assert_eq!(kr.n, 3);

        // CSR esperado:
        // [10 2 0]
        // [2 30 0]
        // [0 0 50]

        let row_ptr_esperado = vec![0, 2, 4, 5];
        let col_ind_esperado = vec![0, 1, 0, 1, 2];
        let values_esperado = vec![10.0, 2.0, 2.0, 30.0, 50.0];

        assert_eq!(kr.row_ptr, row_ptr_esperado);
        assert_eq!(kr.col_ind, col_ind_esperado);
        assert_eq!(kr.values, values_esperado);
    }
    #[test]
fn testa_dense_to_csr() {
    let mat = vec![
        vec![10.0, 0.0, 2.0],
        vec![0.0, 20.0, 0.0],
        vec![2.0, 0.0, 30.0],
    ];

    let csr = dense_to_csr(&mat, 1e-12);

    assert_eq!(csr.row_ptr, vec![0, 2, 3, 5]);
    assert_eq!(csr.col_ind, vec![0, 2, 1, 0, 2]);
    assert_eq!(csr.values, vec![10.0, 2.0, 20.0, 2.0, 30.0]);
}
}