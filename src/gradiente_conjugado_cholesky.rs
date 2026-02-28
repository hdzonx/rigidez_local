// ================================================================
// Solver MEF esparso em Rust
// Implementa:
// 1) Estrutura CSR
// 2) Conversão densa -> CSR (paralela)
// 3) Redução de matriz CSR (remoção de DOFs)
// 4) Operações básicas paralelas (SpMV, dot, axpy)
// 5) Fatoração Incomplete Cholesky ICC(0)
// 6) Aplicação do pré-condicionador ICC
// 7) Gradiente Conjugado pré-condicionado (PCG + ICC)
// 8) Testes unitários
//
// Dependências:
// rayon = "1"
//
// Observação:
// - A paralelização é aplicada nas partes mais custosas:
//   * multiplicação matriz-vetor (SpMV)
//   * produtos internos
//   * operações vetoriais
// - A fatoração ICC(0) é majoritariamente sequencial por natureza.
// ================================================================

use rayon::prelude::*;

// ================================================================
// Estrutura CSR (Compressed Sparse Row)
// Armazena matriz esparsa de dimensão n x n
// ================================================================
#[derive(Clone, Debug)]
pub struct CsrMatrix {
    pub n: usize,          // dimensão da matriz
    pub row_ptr: Vec<usize>, // ponteiros de linha (tamanho n+1)
    pub col_ind: Vec<usize>, // índices de coluna
    pub values: Vec<f64>,    // valores não nulos
}

// ================================================================
// Converte matriz densa (Vec<Vec<f64>>) para CSR paralelamente
// tol: tolerância para eliminar valores muito pequenos (≈ zero)
// ================================================================
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

    // Montagem final CSR (sequencial para manter ordem)
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

// ================================================================
// Reduz matriz CSR removendo DOFs especificados
// removidos: índices globais a eliminar (condições de contorno)
// ================================================================
pub fn reduzir_csr(k: &CsrMatrix, removidos: &[usize]) -> CsrMatrix {
    let n = k.n;

    // máscara de DOFs livres
    let mut livre = vec![true; n];
    for &r in removidos {
        livre[r] = false;
    }

    // mapeamento global -> reduzido
    let mut map = vec![usize::MAX; n];
    let mut count = 0;
    for i in 0..n {
        if livre[i] {
            map[i] = count;
            count += 1;
        }
    }
    let m = count;

    // processamento paralelo por linhas
    let rows: Vec<(Vec<usize>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .filter(|&i| livre[i])
        .map(|i| {
            let mut cols = Vec::new();
            let mut vals = Vec::new();

            for idx in k.row_ptr[i]..k.row_ptr[i + 1] {
                let j = k.col_ind[idx];
                if livre[j] {
                    cols.push(map[j]);
                    vals.push(k.values[idx]);
                }
            }
            (cols, vals)
        })
        .collect();

    // montagem CSR final
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

// ================================================================
// Multiplicação matriz-vetor esparsa paralela: y = A * x
// Operação dominante do método CG
// ================================================================
pub fn spmv_csr(a: &CsrMatrix, x: &[f64]) -> Vec<f64> {
    (0..a.n)
        .into_par_iter()
        .map(|i| {
            let mut sum = 0.0;
            for idx in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_ind[idx];
                sum += a.values[idx] * x[j];
            }
            sum
        })
        .collect()
}

// ================================================================
// Produto interno paralelo: x·y
// ================================================================
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.par_iter().zip(y.par_iter()).map(|(a, b)| a * b).sum()
}

// ================================================================
// Operação vetorial paralela: y = y + a*x
// ================================================================
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    y.par_iter_mut()
        .zip(x.par_iter())
        .for_each(|(yi, xi)| *yi += a * xi);
}

// ================================================================
// Fatoração Incomplete Cholesky ICC(0)
// Mantém mesma estrutura de esparsidade da matriz A
// Assume matriz SPD
// ================================================================
pub fn icc0(a: &CsrMatrix) -> CsrMatrix {
    let n = a.n;
    let mut l = a.clone();

    for i in 0..n {
        for idx in l.row_ptr[i]..l.row_ptr[i + 1] {
            let j = l.col_ind[idx];
            if j >= i { break; }

            let mut sum = l.values[idx];
            let mut p = l.row_ptr[i];
            let mut q = l.row_ptr[j];

            while p < l.row_ptr[i + 1] && q < l.row_ptr[j + 1] {
                let cp = l.col_ind[p];
                let cq = l.col_ind[q];
                if cp == cq && cp < j {
                    sum -= l.values[p] * l.values[q];
                    p += 1;
                    q += 1;
                } else if cp < cq {
                    p += 1;
                } else {
                    q += 1;
                }
            }

            l.values[idx] = sum / l_diag(&l, j);
        }

        // diagonal
        let mut diag_sum = 0.0;
        for idx in l.row_ptr[i]..l.row_ptr[i + 1] {
            let j = l.col_ind[idx];
            if j < i {
                diag_sum += l.values[idx] * l.values[idx];
            } else if j == i {
                l.values[idx] = (l.values[idx] - diag_sum).sqrt();
                break;
            }
        }
    }

    l
}

// ================================================================
// Recupera elemento diagonal L(i,i)
// ================================================================
fn l_diag(l: &CsrMatrix, i: usize) -> f64 {
    for idx in l.row_ptr[i]..l.row_ptr[i + 1] {
        if l.col_ind[idx] == i {
            return l.values[idx];
        }
    }
    panic!("Diagonal não encontrada");
}

// ================================================================
// Resolve M z = r, onde M = L Lᵀ (pré-condicionador ICC)
// ================================================================
pub fn icc_solve(l: &CsrMatrix, r: &[f64]) -> Vec<f64> {
    let n = l.n;
    let mut y = vec![0.0; n];

    // Forward: L y = r
    for i in 0..n {
        let mut sum = r[i];
        let mut diag = 1.0;
        for idx in l.row_ptr[i]..l.row_ptr[i + 1] {
            let j = l.col_ind[idx];
            if j < i {
                sum -= l.values[idx] * y[j];
            } else if j == i {
                diag = l.values[idx];
                break;
            }
        }
        y[i] = sum / diag;
    }

    // construir Lᵀ
    let lt = csr_transpose(l);

    // Backward: Lᵀ z = y
    let mut z = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = y[i];
        let mut diag = 1.0;

        for idx in lt.row_ptr[i]..lt.row_ptr[i + 1] {
            let j = lt.col_ind[idx];
            if j > i {
                sum -= lt.values[idx] * z[j];
            } else if j == i {
                diag = lt.values[idx];
            }
        }
        z[i] = sum / diag;
    }

    z
}

// ================================================================
// Gradiente Conjugado Pré-condicionado (PCG + ICC)
// Resolve A x = b para matriz SPD esparsa
// ================================================================
pub fn pcg_icc(a: &CsrMatrix, b: &[f64], tol: f64, max_iter: usize) -> Vec<f64> {
    let n = a.n;
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();

    let l = icc0(a);
    let mut z = icc_solve(&l, &r);
    let mut p = z.clone();

    let mut rz_old = dot(&r, &z);

    for _ in 0..max_iter {
        let ap = spmv_csr(a, &p);
        let alpha = rz_old / dot(&p, &ap);

        axpy(alpha, &p, &mut x);
        axpy(-alpha, &ap, &mut r);

        if dot(&r, &r).sqrt() < tol {
            break;
        }

        z = icc_solve(&l, &r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;

        p.par_iter_mut()
            .zip(z.par_iter())
            .for_each(|(pi, zi)| *pi = zi + beta * *pi);

        rz_old = rz_new;
    }

    x
}


pub fn csr_transpose(a: &CsrMatrix) -> CsrMatrix {
    let n = a.n;
    let mut nnz_per_col = vec![0usize; n];

    for &j in &a.col_ind {
        nnz_per_col[j] += 1;
    }

    let mut row_ptr = vec![0usize; n + 1];
    for i in 0..n {
        row_ptr[i + 1] = row_ptr[i] + nnz_per_col[i];
    }

    let mut col_ind = vec![0usize; a.col_ind.len()];
    let mut values = vec![0.0; a.values.len()];
    let mut counter = row_ptr.clone();

    for i in 0..n {
        for idx in a.row_ptr[i]..a.row_ptr[i + 1] {
            let j = a.col_ind[idx];
            let dest = counter[j];
            col_ind[dest] = i;
            values[dest] = a.values[idx];
            counter[j] += 1;
        }
    }

    CsrMatrix { n, row_ptr, col_ind, values }
}

// ================================================================
// TESTES UNITÁRIOS
// Verifica consistência do solver PCG + ICC
// ================================================================
#[cfg(test)]
mod tests {
    use super::*;

    fn matriz_teste() -> CsrMatrix {
        // Matriz SPD:
        // [4 1 0]
        // [1 3 1]
        // [0 1 2]
        CsrMatrix {
            n: 3,
            row_ptr: vec![0, 2, 5, 7],
            col_ind: vec![0,1, 0,1,2, 1,2],
            values: vec![4.0,1.0, 1.0,3.0,1.0, 1.0,2.0],
        }
    }

    #[test]
    fn testa_pcg_icc() {
        let a = matriz_teste();
        let b = vec![1.0, 2.0, 0.0];

        let x = pcg_icc(&a, &b, 1e-8, 100);

        // solução exata aproximada
        let x_ref = vec![0.09090909, 0.63636363, -0.31818181];

        for i in 0..3 {
            assert!((x[i] - x_ref[i]).abs() < 1e-4);
        }
    }
}