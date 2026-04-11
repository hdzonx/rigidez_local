// use nalgebra::DMatrix;
// use rayon::prelude::*;
// #[derive(Debug)]
// pub struct CsrMatrix {
//     pub values: Vec<f64>,
//     col_indices: Vec<usize>,
//     row_ptr: Vec<usize>,
//     n: usize,
// }

// pub fn dense_to_csr(mat: &DMatrix<f64>, tol: f64) -> CsrMatrix {
//     let nrows = mat.nrows();
//     let ncols = mat.ncols();

//     // Cada thread acumula contribuições por linha
//     let per_thread: Vec<Vec<(usize, f64)>> = (0..ncols)
//         .into_par_iter()
//         .map(|j| {
//             let mut local: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nrows];

//             // Acesso contíguo à coluna
//             let col = mat.column(j);

//             for (i, &v) in col.iter().enumerate() {
//                 if v.abs() > tol {
//                     local[i].push((j, v));
//                 }
//             }

//             local
//         })
//         .reduce(
//             || vec![Vec::new(); nrows],
//             |mut acc, mut local| {
//                 for i in 0..nrows {
//                     acc[i].extend(local[i].drain(..));
//                 }
//                 acc
//             },
//         );

//     // Montagem CSR
//     let mut row_ptr = Vec::with_capacity(nrows + 1);
//     let mut col_indices = Vec::new();
//     let mut values = Vec::new();

//     row_ptr.push(0);

//     for row in per_thread {
//         for (j, v) in row {
//             col_indices.push(j);
//             values.push(v);
//         }
//         row_ptr.push(col_indices.len());
//     }

//     CsrMatrix {
//         values,
//         col_indices,
//         row_ptr,
//         n: nrows,
//     }
// }

// pub fn reduzir_csr(k: &CsrMatrix, removidos: &[usize]) -> CsrMatrix {
//     let n = k.n;

//     // 1. máscara
//     let mut livre = vec![true; n];
//     for &r in removidos {
//         livre[r] = false;
//     }

//     // 2. mapeamento global -> reduzido
//     let mut map = vec![usize::MAX; n];
//     let mut count = 0;
//     for i in 0..n {
//         if livre[i] {
//             map[i] = count;
//             count += 1;
//         }
//     }
//     let m = count;

//     // 3. Processamento paralelo por linhas
//     let rows: Vec<(Vec<usize>, Vec<f64>)> = (0..n)
//         .into_par_iter()
//         .filter(|&i| livre[i])
//         .map(|i| {
//             let mut cols = Vec::new();
//             let mut vals = Vec::new();

//             let start = k.row_ptr[i];
//             let end = k.row_ptr[i + 1];

//             for idx in start..end {
//                 let j = k.col_indices[idx];
//                 if livre[j] {
//                     cols.push(map[j]);
//                     vals.push(k.values[idx]);
//                 }
//             }

//             (cols, vals)
//         })
//         .collect();

//     // 4. Montagem CSR final
//     let mut row_ptr = Vec::with_capacity(m + 1);
//     let mut col_indices = Vec::new();
//     let mut values = Vec::new();

//     row_ptr.push(0);

//     for (cols, vals) in rows {
//         col_indices.extend(cols);
//         values.extend(vals);
//         row_ptr.push(col_indices.len());
//     }

//     CsrMatrix {
//         values,
//         col_indices,
//         row_ptr,
//         n: m,
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use nalgebra::DMatrix;

//     fn approx_eq(a: &[f64], b: &[f64], tol: f64) -> bool {
//         a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| (x - y).abs() < tol)
//     }

//     #[test]
//     fn test_small_matrix() {
//         let mat = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0]);

//         let csr = dense_to_csr(&mat, 1e-12);

//         assert_eq!(csr.n, 3);
//         assert_eq!(csr.row_ptr, vec![0, 2, 3, 4]);
//         assert_eq!(csr.col_indices, vec![0, 2, 2, 0]);
//         assert!(approx_eq(&csr.values, &[1.0, 2.0, 3.0, 4.0], 1e-12));
//     }

//     #[test]
//     fn test_zero_matrix() {
//         let mat = DMatrix::<f64>::zeros(4, 4);

//         let csr = dense_to_csr(&mat, 1e-12);

//         assert_eq!(csr.n, 4);
//         assert_eq!(csr.row_ptr, vec![0, 0, 0, 0, 0]);
//         assert!(csr.col_indices.is_empty());
//         assert!(csr.values.is_empty());
//     }

//     #[test]
//     fn test_dense_matrix() {
//         let mat = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);

//         let csr = dense_to_csr(&mat, 1e-12);

//         assert_eq!(csr.row_ptr, vec![0, 2, 4]);
//         assert_eq!(csr.col_indices, vec![0, 1, 0, 1]);
//         assert!(approx_eq(&csr.values, &[1.0, 2.0, 3.0, 4.0], 1e-12));
//     }

//     #[test]
//     fn test_tolerance_filtering() {
//         let mat = DMatrix::from_row_slice(2, 3, &[1e-10, 1.0, 0.0, 0.0, 1e-11, 2.0]);

//         let csr = dense_to_csr(&mat, 1e-9);

//         assert_eq!(csr.row_ptr, vec![0, 1, 2]);
//         assert_eq!(csr.col_indices, vec![1, 2]);
//         assert!(approx_eq(&csr.values, &[1.0, 2.0], 1e-12));
//     }

//     #[test]
//     fn test_rectangular_matrix() {
//         let mat = DMatrix::from_row_slice(2, 4, &[0.0, 5.0, 0.0, 6.0, 7.0, 0.0, 0.0, 8.0]);

//         let csr = dense_to_csr(&mat, 1e-12);

//         assert_eq!(csr.row_ptr, vec![0, 2, 4]);
//         assert_eq!(csr.col_indices, vec![1, 3, 0, 3]);
//         assert!(approx_eq(&csr.values, &[5.0, 6.0, 7.0, 8.0], 1e-12));
//     }

//     #[test]
//     fn test_row_ptr_monotonicity() {
//         let mat = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]);

//         println!("matriz = {}", mat);

//         let csr = dense_to_csr(&mat, 1e-12);

//         // row_ptr deve ser não-decrescente
//         for i in 0..csr.row_ptr.len() - 1 {
//             assert!(csr.row_ptr[i] <= csr.row_ptr[i + 1]);
//         }

//         // último valor = número total de nnz
//         assert_eq!(csr.row_ptr.last().copied().unwrap(), csr.values.len());
//     }
//     #[test]
//     fn reducao_matriz_rig_global() {
//         let matriz_reduzida_esperada = DMatrix::<f64>::from_row_slice(
//             4,
//             4,
//             &[
//                 2.11e3, 2.81e3, -1.88e3, 3.75e3, 2.81e3, 15.0e3, -3.75e3, 5.0e3, -1.88e3, -3.75e3,
//                 1.88e3, -3.75e3, 3.75e3, 5.0e3, -3.75e3, 10.0e3,
//             ],
//         );

//         let rigidez_global_densa = DMatrix::<f64>::from_row_slice(
//             6,
//             6,
//             &[
//                 0.23e3, 0.94e3, -0.23e3, 0.94e3, 0.0, 0.0, 0.94e3, 5.0e3, -0.94e3, 2.50e3, 0.0,
//                 0.0, -0.23e3, -0.94e3, 2.11e3, 2.81e3, -1.88e3, 3.75e3, 0.94e3, 2.5e3, 2.81e3,
//                 15.0e3, -3.75e3, 5.0e3, 0.0, 0.0, -1.88e3, -3.75e3, 1.88e3, -3.75e3, 0.0, 0.0,
//                 3.75e3, 5.0e3, -3.75e3, 10.0e3,
//             ],
//         );

//         //Transforma a matriz densa (esparsa) numa matriz csr
//         let rigidez_global_csr = dense_to_csr(&rigidez_global_densa, 1e-12);

//         println!("rigidez global csr = {:?}", rigidez_global_csr);

//         // Remover DOFs 1 e 3
//         let removidos = vec![0, 1];

//         let k_reduzida = reduzir_csr(&rigidez_global_csr, &removidos);

//         println!("rigidez global reduzida = {:?}", k_reduzida);

//         //Avalie o erro com assertion para cada valor da matriz
//         let tol = 0.01;
//         let mut n = 0;
//         for i in 0..matriz_reduzida_esperada.nrows() {
//             for j in 0..matriz_reduzida_esperada.ncols() {
//                 println!("i = {}, j ={}", i + 1, j + 1);
//                 let relat_err = matriz_reduzida_esperada[(i, j)].abs() - k_reduzida.values[n].abs();
//                 println!("valor esperado = {}", matriz_reduzida_esperada[(i, j)]);
//                 println!("valor calculado = {}", k_reduzida.values[n]);

//                 n += 1;

//                 assert!(
//                     relat_err < tol,
//                     "Erro relativo alto demais: {} (esperado < {})",
//                     relat_err,
//                     tol
//                 );
//             }
//         }
//     }

//     #[test]
//     fn reducao_matriz_rig_global_2() {
//         let matriz_reduzida_esperada =
//             DMatrix::<f64>::from_row_slice(2, 2, &[0.23e3, 0.94e3, 0.94e3, 5.0e3]);

//         let rigidez_global_densa = DMatrix::<f64>::from_row_slice(
//             6,
//             6,
//             &[
//                 0.23e3, 0.94e3, -0.23e3, 0.94e3, 0.0, 0.0, 0.94e3, 5.0e3, -0.94e3, 2.50e3, 0.0,
//                 0.0, -0.23e3, -0.94e3, 2.11e3, 2.81e3, -1.88e3, 3.75e3, 0.94e3, 2.5e3, 2.81e3,
//                 15.0e3, -3.75e3, 5.0e3, 0.0, 0.0, -1.88e3, -3.75e3, 1.88e3, -3.75e3, 0.0, 0.0,
//                 3.75e3, 5.0e3, -3.75e3, 10.0e3,
//             ],
//         );

//         //Transforma a matriz densa (esparsa) numa matriz csr
//         let rigidez_global_csr = dense_to_csr(&rigidez_global_densa, 1e-12);

//         println!("rigidez global csr = {:?}", rigidez_global_csr);

//         // Remover DOFs 1 e 3
//         let removidos = vec![2, 3, 4, 5];

//         let k_reduzida = reduzir_csr(&rigidez_global_csr, &removidos);

//         println!("rigidez global reduzida = {:?}", k_reduzida);

//         //Avalie o erro com assertion para cada valor da matriz
//         let tol = 0.01;
//         let mut n = 0;
//         for i in 0..matriz_reduzida_esperada.nrows() {
//             for j in 0..matriz_reduzida_esperada.ncols() {
//                 println!("i = {}, j ={}", i + 1, j + 1);
//                 let relat_err = matriz_reduzida_esperada[(i, j)].abs() - k_reduzida.values[n].abs();
//                 println!("valor esperado = {}", matriz_reduzida_esperada[(i, j)]);
//                 println!("valor calculado = {}", k_reduzida.values[n]);

//                 n += 1;

//                 assert!(
//                     relat_err < tol,
//                     "Erro relativo alto demais: {} (esperado < {})",
//                     relat_err,
//                     tol
//                 );
//             }
//         }
//     }
// }
