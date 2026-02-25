use rayon::prelude::*;

/// Estrutura CSR (Compressed Sparse Row)
#[derive(Debug)]
struct CsrMatrix {
    values: Vec<f64>,
    col_indices: Vec<usize>,
    row_ptr: Vec<usize>,
    n: usize,
}

impl CsrMatrix {
    /// Produto matriz-vetor paralelo
    fn matvec(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n)
            .into_par_iter()
            .map(|i| {
                let start = self.row_ptr[i];
                let end = self.row_ptr[i + 1];

                let mut sum = 0.0;
                for idx in start..end {
                    sum += self.values[idx] * x[self.col_indices[idx]];
                }
                sum
            })
            .collect()
    }
}

/// Produto interno paralelo
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.par_iter()
        .zip(b.par_iter())
        .map(|(x, y)| x * y)
        .sum()
}

/// Gradiente Conjugado
fn conjugate_gradient(
    a: &CsrMatrix,
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];
    let mut r = b.to_vec(); // r = b - Ax (x=0 → r=b)
    let mut p = r.clone();
    let mut rs_old = dot(&r, &r);

    for _ in 0..max_iter {
        let ap = a.matvec(&p);
        let alpha = rs_old / dot(&p, &ap);

        // x = x + alpha * p
        x.par_iter_mut()
            .zip(p.par_iter())
            .for_each(|(xi, pi)| *xi += alpha * pi);

        // r = r - alpha * Ap
        r.par_iter_mut()
            .zip(ap.par_iter())
            .for_each(|(ri, api)| *ri -= alpha * api);

        let rs_new = dot(&r, &r);

        if rs_new.sqrt() < tol {
            break;
        }

        let beta = rs_new / rs_old;

        // p = r + beta * p
        p.par_iter_mut()
            .zip(r.par_iter())
            .for_each(|(pi, ri)| *pi = ri + beta * (*pi));

        rs_old = rs_new;
    }

    x
}

#[cfg(test)]
mod test{
    use crate::gradiente_conjugado::CsrMatrix;
    use crate::gradiente_conjugado::conjugate_gradient;
    #[test]
    fn conjugate_test(){
    // Matriz SPD simples (4x4)
    // Exemplo clássico
    let values = vec![
        4.0, 1.0,
        1.0, 3.0, 1.0,
        1.0, 2.0,
        2.0
    ];

    let col_indices = vec![
        0, 1,
        0, 1, 2,
        1, 2,
        3
    ];

    let row_ptr = vec![0, 2, 5, 7, 8];

    let a = CsrMatrix {
        values,
        col_indices,
        row_ptr,
        n: 4,
    };

    let b = vec![1.0, 2.0, 3.0, 4.0];

    let solution = conjugate_gradient(&a, &b, 100, 1e-10);

    println!("Solução aproximada:");
    for (i, val) in solution.iter().enumerate() {
        println!("x{} = {:.6}", i, val);
    }
    }

}