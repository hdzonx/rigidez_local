use rayon::prelude::*;
///Gradiente conjugado com pré-condicionado Jacobi
#[derive(Debug, Clone)]
pub struct CsrMatrix {
    values: Vec<f64>,
    col_indices: Vec<usize>,
    row_ptr: Vec<usize>,
    n: usize,
}

impl CsrMatrix {
    pub fn from_dense(dense: &[Vec<f64>]) -> Self {
        let n = dense.len();
        assert!(dense.iter().all(|r| r.len() == n));

        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptr = Vec::with_capacity(n + 1);

        row_ptr.push(0);

        for row in dense {
            for (j, &val) in row.iter().enumerate() {
                if val != 0.0 {
                    values.push(val);
                    col_indices.push(j);
                }
            }
            row_ptr.push(values.len());
        }

        Self {
            values,
            col_indices,
            row_ptr,
            n,
        }
    }

    /// Extrai diagonal (para Jacobi)
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n)
            .into_par_iter()
            .map(|i| {
                let start = self.row_ptr[i];
                let end = self.row_ptr[i + 1];

                for idx in start..end {
                    if self.col_indices[idx] == i {
                        return self.values[idx];
                    }
                }
                panic!("Matriz não possui elemento diagonal na linha {}", i);
            })
            .collect()
    }

    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
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

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.par_iter().zip(b.par_iter()).map(|(x, y)| x * y).sum()
}

/// Gradiente Conjugado Precondicionado (Jacobi)
pub fn conjugate_gradient_jacobi(a: &CsrMatrix, b: &[f64], max_iter: usize, tol: f64) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];

    let mut r = b.to_vec();
    let diag = a.diagonal();

    // z = M⁻¹ r
    let mut z: Vec<f64> = r
        .par_iter()
        .zip(diag.par_iter())
        .map(|(ri, di)| ri / di)
        .collect();

    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);

    for _ in 0..max_iter {
        let ap = a.matvec(&p);

        let alpha = rz_old / dot(&p, &ap);

        x.par_iter_mut()
            .zip(p.par_iter())
            .for_each(|(xi, pi)| *xi += alpha * pi);

        r.par_iter_mut()
            .zip(ap.par_iter())
            .for_each(|(ri, api)| *ri -= alpha * api);

        if dot(&r, &r).sqrt() < tol {
            break;
        }

        // z = M⁻¹ r
        z.par_iter_mut()
            .zip(r.par_iter())
            .zip(diag.par_iter())
            .for_each(|((zi, ri), di)| *zi = ri / di);

        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;

        p.par_iter_mut()
            .zip(z.par_iter())
            .for_each(|(pi, zi)| *pi = zi + beta * (*pi));

        rz_old = rz_new;
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pcg_jacobi() {
        let dense = vec![
            vec![4.0, 1.0, 0.0, 0.0],
            vec![1.0, 3.0, 1.0, 0.0],
            vec![0.0, 1.0, 2.0, 1.0],
            vec![0.0, 0.0, 1.0, 2.0],
        ];

        println!("Threads Rayon used: {}", rayon::current_num_threads());
        let a = CsrMatrix::from_dense(&dense);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        println!("a = {:?}", a);
        println!("b = {:?}", b);

        let x = conjugate_gradient_jacobi(&a, &b, 1000, 1e-12);

        let ax = a.matvec(&x);

        println!("ax = {:?}", ax);

        let error = ax
            .iter()
            .zip(b.iter())
            .map(|(ax_i, b_i)| (ax_i - b_i).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro = {}", error);
        assert!(error < 1e-8);
    }

    #[test]
    fn test_pcg_jacobi_02() {
        let dense = vec![
            vec![1764100.0, -807300.0, 358800.0],
            vec![-807300.0, 1764100.0, 0.0],
            vec![358800.0, 0.0, 2511600.0],
        ];

        println!("Threads Rayon used: {}", rayon::current_num_threads());
        //Sistema [a] [x] =[b]
        let a = CsrMatrix::from_dense(&dense);
        let b = vec![0.0, 0.0, -4450.0];
        println!("a = {:?}", a);
        println!("b = {:?}", b);

        let x = conjugate_gradient_jacobi(&a, &b, 1000, 1e-12);
        println!("x = {:?}", x);

        let ax = a.matvec(&x);
        println!("ax = {:?}", ax);

        let error = ax
            .iter()
            .zip(b.iter())
            .map(|(ax_i, b_i)| (ax_i - b_i).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro = {}", error);
        assert!(error < 1e-8);
    }

    #[test]
    fn test_pcg_jacobi_03() {
        let scale = 3.297e6;

        let dense = vec![
            vec![1.3 * scale, -0.488 * scale, -0.55 * scale, 0.313 * scale],
            vec![-0.488 * scale, 0.894 * scale, 0.338 * scale, -0.631 * scale],
            vec![-0.55 * scale, 0.338 * scale, 1.3 * scale, -0.163 * scale],
            vec![0.313 * scale, -0.631 * scale, -0.163 * scale, 0.894 * scale],
        ];

        println!("Threads Rayon used: {}", rayon::current_num_threads());
        //Sistema [a] [x] = [b]
        let a = CsrMatrix::from_dense(&dense);
        let b = vec![0.0, -50000.0, 50000.0, 0.0];
        println!("a = {:?}", a);
        println!("b = {:?}", b);

        let x = conjugate_gradient_jacobi(&a, &b, 1000, 1e-12);
        println!("x = {:?}", x);

        // Multiplicação [a][x]
        let ax = a.matvec(&x);
        println!("ax = {:?}", ax);

        let error = ax
            .iter()
            .zip(b.iter())
            .map(|(ax_i, b_i)| (ax_i - b_i).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro = {}", error);
        assert!(error < 1e-8);
    }

    #[test]
    fn test_pcg_jacobi_result_x() {
        let scale = 3.297e6;

        let dense = vec![
            vec![1.3 * scale, -0.488 * scale, -0.55 * scale, 0.313 * scale],
            vec![-0.488 * scale, 0.894 * scale, 0.338 * scale, -0.631 * scale],
            vec![-0.55 * scale, 0.338 * scale, 1.3 * scale, -0.163 * scale],
            vec![0.313 * scale, -0.631 * scale, -0.163 * scale, 0.894 * scale],
        ];

        println!("Threads Rayon used: {}", rayon::current_num_threads());
        //Sistema [a] [x] = [b]
        let a = CsrMatrix::from_dense(&dense);
        let b = vec![0.0, -50000.0, 50000.0, 0.0];
        println!("a = {:?}", a);
        println!("b = {:?}", b);

        let x = conjugate_gradient_jacobi(&a, &b, 1000, 1e-12);
        println!("x = {:?}", x);

        let x_expected: Vec<f64> = vec![-2.14863e-3, -4.4481e-2, 1.89118e-2, -2.7195e-2];
        println!(" expected x = {:?}", x_expected);

        let error = x
            .iter()
            .zip(x_expected.iter())
            .map(|(x_i, x_exp_i)| (x_i - x_exp_i).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro = {}", error);
        assert!(error < 1e-6);
    }
}
