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

mod test {
    use nalgebra::DMatrix;
    use nalgebra::Matrix6;
    use nalgebra::SMatrix;

    use crate::area_triangulo;
    use crate::constitutive_matrix;
    use crate::gradiente_conjug_jacobi;
    use crate::matriz_rigidez_local;
    use crate::rigidez_global;

    use sprs::{CsMat, TriMat};

    #[test]
    fn estrutura_completa_teste_1() {
        let x_coords_elem_1 = vec![0.0, 10.0, 10.0];
        let y_coords_elem_1 = vec![0.0, 5.0, 15.0];

        let x_coords_elem_2 = vec![0.0, 10.0, 0.0];
        let y_coords_elem_2 = vec![0.0, 15.0, 20.0];

        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let rigidez_local_element_01 = matriz_rigidez_local(
            x_coords_elem_1,
            y_coords_elem_1,
            espessura,
            constitutive_matrix,
        );

        let rigidez_local_element_02 = matriz_rigidez_local(
            x_coords_elem_2,
            y_coords_elem_2,
            espessura,
            constitutive_matrix,
        );

        let elements = vec![
            ([0, 1, 2], rigidez_local_element_01),
            ([0, 2, 3], rigidez_local_element_02),
        ];

        let k_global = rigidez_global::assemble_sparse(4, &elements);

        // Remover DOFs 1 e 3
        let removidos = vec![0, 1, 6, 7];

        let k_reduzida =
            gradiente_conjug_jacobi::reduzir_matriz_k_global_csr(&k_global, &removidos);

        let b = vec![0.0, -50000.0, 50000.0, 0.0];
        let x = gradiente_conjug_jacobi::conjugate_gradient_jacobi(&k_reduzida, &b, 1000, 1e-12);
        println!("x = {:?}", x);
        //Avaliando o erro do vetor X
        //a.x = b
        let x_expected = [-2.147e-3, -4.455e-2, 1.891e-2, -2.727e-2];
        let error_x = x
            .iter()
            .zip(x_expected.iter())
            .map(|(x_i, x_i_exp)| (x_i - x_i_exp).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro de X = {}", error_x);
        assert!(error_x < 1e-3);
    }
}
