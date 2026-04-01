// Elemento CST (Triângulo de deformação constante)
//
//      |\ (4)-[0,20]
//      |  \
//      |    \
//      |      \
//      |   2    \ (3)-[10,15]
//      |       / |
//      |      /  |
//      |     /   |
//      |    /  1 | (2)-[10,5]
//      |   /    /
//      |  /   /
//      | /  /
//      |/ /
//      |/ (1)-[0,0]
use sprs::{CsMat, TriMat};
#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    //use approx::relative_eq;
    use nalgebra::DMatrix;
    use nalgebra::Matrix3;
    use nalgebra::Matrix6;
    use nalgebra::SMatrix;

    use crate::area_triangulo;
    use crate::constitutive_matrix;
    use crate::matriz_reduzida_eficiente;
    use crate::matriz_rigidez_local;
    use crate::rigidez_global;

    use sprs::{CsMat, TriMat};

    #[test]
    fn calcula_rigidez_local_el_1() {
        let x_coords = vec![0.0, 10.0, 10.0];
        let y_coords = vec![0.0, 5.0, 15.0];

        let area = area_triangulo(&x_coords, &y_coords);
        let espessura = 0.1;

        let escalar = 3.297e6;

        let matriz_nua: Matrix6<f64> = Matrix6::new(
            0.5, 0.0, -0.75, 0.15, 0.25, -0.15, 0., 0.175, 0.175, -0.2627, -0.175, 0.088, -0.75,
            0.1750, 1.3, -0.488, -0.55, 0.313, 0.15, -0.2627, -0.488, 0.894, 0.338, -0.631, 0.25,
            -0.175, -0.55, 0.338, 0.3, -0.163, -0.15, 0.088, 0.313, -0.631, -0.163, 0.544,
        );

        let rigidez_local_esperada = escalar * matriz_nua;

        let poisson = 0.3;
        let elasticidade = 30.0e6;
        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let rigidez_local_calc =
            matriz_rigidez_local(x_coords, y_coords, espessura, constitutive_matrix);

        //Avalie o erro com assertion para cada valor da matriz
        for i in 0..rigidez_local_esperada.nrows() {
            for j in 0..rigidez_local_esperada.ncols() {
                let tol = 1e4;
                println!("i = {}, j ={}", i + 1, j + 1);
                let relat_err = (rigidez_local_esperada[(i, j)] - rigidez_local_calc[(i, j)]).abs();

                println!("calculado = {}", rigidez_local_calc[(i, j)]);
                println!("esperado = {}", rigidez_local_esperada[(i, j)]);

                assert!(
                    relat_err < tol,
                    "Erro relativo alto demais: {} (esperado < {})",
                    relat_err,
                    tol
                );
            }
        }

        println!("{:?}", rigidez_local_calc);
    }

    #[test]
    fn calcula_rigidez_local_el_2() {
        let x_coords = vec![0.0, 10.0, 0.0];
        let y_coords = vec![0.0, 15.0, 20.0];

        let area = area_triangulo(&x_coords, &y_coords);
        let espessura = 0.1;

        let escalar = 3.297e6;

        let matriz_nua: Matrix6<f64> = Matrix6::new(
            0.15, 0.081, -0.25, -0.175, 0.1, 0.094, 0.081, 0.272, -0.15, -0.088, 0.069, -0.184,
            -0.25, -0.15, 1.0, 0.0, -0.75, 0.15, -0.175, -0.088, 0.0, 0.35, 0.175, -0.263, 0.1,
            0.069, -0.75, 0.175, 0.65, -0.244, 0.094, -0.184, 0.15, -0.263, -0.244, 0.447,
        );

        let rigidez_local_esperada = escalar * matriz_nua;

        let poisson = 0.3;
        let elasticidade = 30.0e6;
        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let rigidez_local_calc =
            matriz_rigidez_local(x_coords, y_coords, espessura, constitutive_matrix);

        //Avalie o erro com assertion para cada valor da matriz
        for i in 0..rigidez_local_esperada.nrows() {
            for j in 0..rigidez_local_esperada.ncols() {
                let tol = 1e4;
                println!("i = {}, j ={}", i + 1, j + 1);
                let relat_err = (rigidez_local_esperada[(i, j)] - rigidez_local_calc[(i, j)]).abs();

                println!("calculado = {}", rigidez_local_calc[(i, j)]);
                println!("esperado = {}", rigidez_local_esperada[(i, j)]);

                assert!(
                    relat_err < tol,
                    "Erro relativo alto demais: {} (esperado < {})",
                    relat_err,
                    tol
                );
            }
        }

        println!("{:?}", rigidez_local_calc);
    }

    #[test]
    fn test_rigidez_global_esparsa() {
        //Matriz de rigidez global esperada
        let escalar = 3.297e6;
        let matriz_nua = SMatrix::<f64, 8, 8>::from_row_slice(&[
            0.65, 0.081, -0.75, 0.15, 0.0, -0.325, 0.1, 0.094, 0.081, 0.447, 0.175, -0.263, -0.325,
            0.0, 0.069, -0.184, -0.75, 0.175, 1.3, -0.488, -0.55, 0.313, 0.0, 0.0, 0.15, -0.263,
            -0.488, 0.894, 0.338, -0.631, 0.0, 0.0, 0.0, -0.325, -0.55, 0.338, 1.3, -0.163, -0.75,
            0.15, -0.325, 0.0, 0.313, -0.631, -0.163, 0.894, 0.175, -0.263, 0.1, 0.069, 0.0, 0.0,
            -0.75, 0.175, 0.65, -0.244, 0.094, -0.184, 0.0, 0.0, 0.15, -0.263, -0.244, 0.447,
        ]);

        let rigidez_global_esperada = escalar * matriz_nua;

        //Matriz constitutiva
        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;
        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        //Dados do elemento 1
        let x_coords_ele_01 = vec![0.0, 10.0, 10.0];
        let y_coords_ele_01 = vec![0.0, 5.0, 15.0];
        let rigidez_local_element_01 = matriz_rigidez_local(
            x_coords_ele_01,
            y_coords_ele_01,
            espessura,
            constitutive_matrix,
        );

        //Dados do elemento 2
        let x_coords_el_02 = vec![0.0, 10.0, 0.0];
        let y_coords_el_02 = vec![0.0, 15.0, 20.0];
        let rigidez_local_element_02 = matriz_rigidez_local(
            x_coords_el_02,
            y_coords_el_02,
            espessura,
            constitutive_matrix,
        );

        //Mapeamento de cada nó em cada elemento.
        //O índice i = nó-1
        let elements = vec![
            ([0, 1, 2], rigidez_local_element_01),
            ([0, 2, 3], rigidez_local_element_02),
        ];

        //Cálculo da matriz de rigidez global
        let k_global = rigidez_global::assemble_sparse(4, &elements);

        //Imprime os resultados obtidos
        println!("Matriz global (CSR):");
        println!("{:?}", k_global);
        println!("\nConvertendo para matriz densa para visualização:\n");
        let dense_global_k = k_global.to_dense();
        println!("{}", dense_global_k);

        //Avalie o erro com assertion para cada valor da matriz
        let tol = 1e4;
        for i in 0..rigidez_global_esperada.nrows() {
            for j in 0..rigidez_global_esperada.ncols() {
                println!("i = {}, j ={}", i + 1, j + 1);
                let relat_err = (rigidez_global_esperada[(i, j)] - dense_global_k[(i, j)]).abs();

                println!("calculado = {}", dense_global_k[(i, j)]);
                println!("esperado = {}", rigidez_global_esperada[(i, j)]);

                assert!(
                    relat_err < tol,
                    "Erro relativo alto demais: {} (esperado < {})",
                    relat_err,
                    tol
                );
            }
        }
    }

    #[test]
    fn reducao_matriz_rig_global() {
        let escalar = 3.297e6;
        let matriz_reduzida_nua = DMatrix::<f64>::from_row_slice(
            4,
            4,
            &[
                1.3, -0.488, -0.55, 0.313, -0.488, 0.894, 0.338, -0.631, -0.55, 0.338, 1.3, -0.163,
                0.313, -0.631, -0.163, 0.894,
            ],
        );

        let matriz_reduzida_esperada = escalar * matriz_reduzida_nua;

        let matriz_nua = DMatrix::<f64>::from_row_slice(
            8,
            8,
            &[
                0.65, 0.081, -0.75, 0.15, 0.0, -0.325, 0.1, 0.094, 0.081, 0.447, 0.175, -0.263,
                -0.325, 0.0, 0.069, -0.184, -0.75, 0.175, 1.3, -0.488, -0.55, 0.313, 0.0, 0.0,
                0.15, -0.263, -0.488, 0.894, 0.338, -0.631, 0.0, 0.0, 0.0, -0.325, -0.55, 0.338,
                1.3, -0.163, -0.75, 0.15, -0.325, 0.0, 0.313, -0.631, -0.163, 0.894, 0.175, -0.263,
                0.1, 0.069, 0.0, 0.0, -0.75, 0.175, 0.65, -0.244, 0.094, -0.184, 0.0, 0.0, 0.15,
                -0.263, -0.244, 0.447,
            ],
        );
        let rigidez_global_densa = escalar * matriz_nua;

        //Transforma a matriz densa (esparsa) numa matriz csr
        let rigidez_global_csr =
            matriz_reduzida_eficiente::dense_to_csr(&rigidez_global_densa, 1e-12);

        println!("rigidez global csr = {:?}", rigidez_global_csr);

        // Remover DOFs 1 e 3
        let removidos = vec![0, 1, 6, 7];

        let k_reduzida = matriz_reduzida_eficiente::reduzir_csr(&rigidez_global_csr, &removidos);

        println!("rigidez global reduzida = {:?}", k_reduzida);

        //Avalie o erro com assertion para cada valor da matriz
        let tol = 10.0;
        let mut valor_esperado = 0.0;
        let mut valor_calculado = 0.0;
        let mut n = 0;
        for i in 0..matriz_reduzida_esperada.nrows() {
            for j in 0..matriz_reduzida_esperada.ncols() {
                println!("i = {}, j ={}", i + 1, j + 1);
                let relat_err = matriz_reduzida_esperada[(i, j)].abs() - k_reduzida.values[n].abs();
                println!("valor esperado = {}", matriz_reduzida_esperada[(i, j)]);
                println!("valor calculado = {}", k_reduzida.values[n]);

                n += 1;

                assert!(
                    relat_err < tol,
                    "Erro relativo alto demais: {} (esperado < {})",
                    relat_err,
                    tol
                );
            }
        }
    }
}
