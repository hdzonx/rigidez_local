// Elemento CST (Triângulo de deformação constante)
// A contagem dos nós começa em zero
//
//    X-|\ (3)-[0,20]
//      |  \
//      |    \
//      |      \
//      |   2    \ (2)-[10,15]  -----> 50000
//      |       / |
//      |      /  |
//      |     /   |
//      |    /  1 | (1)-[10,5]
//      |   /    /:
//      |  /   /  :
//      | /  /    V 50000
//      |/ /
//    X-|/ (0)-[0,0]

use sprs::{CsMat, TriMat};
#[cfg(test)]
#[allow(unused)]
mod test {
    use nalgebra::DMatrix;
    use nalgebra::Matrix6;
    use nalgebra::SMatrix;

    use crate::constitutive_matrix;
    use crate::gradiente_conjug_jacobi;
    use crate::rigidez_global;
    use crate::rigidez_local_struct;

    use sprs::{CsMat, TriMat};

    #[test]
    fn estrutura_completa_teste_1() {
        let x_coords_elem_1 = vec![0.0, 10.0, 10.0];
        let y_coords_elem_1 = vec![0.0, 5.0, 15.0];
        let global_node_ele_1: Vec<usize> = vec![0, 1, 2];

        let x_coords_elem_2 = vec![0.0, 10.0, 0.0];
        let y_coords_elem_2 = vec![0.0, 15.0, 20.0];
        let global_node_ele_2: Vec<usize> = vec![0, 2, 3];

        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let mut rigidez_local_element_01 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_1,
            y_coords_elem_1,
            global_node_ele_1,
            espessura,
            constitutive_matrix,
        );

        let mut rigidez_local_element_02 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_2,
            y_coords_elem_2,
            global_node_ele_2,
            espessura,
            constitutive_matrix,
        );

        rigidez_local_element_01.assemble_local();
        rigidez_local_element_02.assemble_local();

        let elements = vec![
            ([0, 1, 2], rigidez_local_element_01),
            ([0, 2, 3], rigidez_local_element_02),
        ];

        let k_global = rigidez_global::assemble_sparse(4, &elements);

        //F = [Rx0, Ry0, 0.0, -50.000, 50.000, 0.0, Rx3, Ry3]
        // Remover DOFs 1 e 3 : Rx0, Ry0, Rx3, Ry3
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

    #[test]
    fn teste_mudando_numeracao_01() {
        // Elemento CST (Triângulo de deformação constante)
        // A contagem dos nós começa em zero
        //
        //    X-|\ (0)-[0,20]
        //      |  \
        //      |    \
        //      |      \
        //      |   1    \ (3)-[10,15]  -----> 50000
        //      |       / |
        //      |      /  |
        //      |     /   |
        //      |    /  2 | (2)-[10,5]
        //      |   /    /:
        //      |  /   /  :
        //      | /  /    V 50000
        //      |/ /
        //    X-|/ (1)-[0,0]

        let x_coords_elem_1 = vec![0.0, 0.0, 10.0];
        let y_coords_elem_1 = vec![20.0, 0.0, 15.0];
        let global_node_ele_1: Vec<usize> = vec![0, 1, 3];

        let x_coords_elem_2 = vec![0.0, 10.0, 10.0];
        let y_coords_elem_2 = vec![0.0, 5.0, 15.0];
        let global_node_ele_2: Vec<usize> = vec![1, 2, 3];

        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let mut rigidez_local_element_01 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_1,
            y_coords_elem_1,
            global_node_ele_1,
            espessura,
            constitutive_matrix,
        );

        let mut rigidez_local_element_02 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_2,
            y_coords_elem_2,
            global_node_ele_2,
            espessura,
            constitutive_matrix,
        );

        rigidez_local_element_01.assemble_local();
        rigidez_local_element_02.assemble_local();

        let elements = vec![
            ([0, 1, 3], rigidez_local_element_01),
            ([1, 2, 3], rigidez_local_element_02),
        ];
        let k_global = rigidez_global::assemble_sparse(4, &elements);

        //F = [Rx0, Ry0, Rx1, Ry1, 0.0, -50.000, 50.000, 0.0]
        // Remover DOFs 0 e 1: Rx0, Ry0, Rx1, Ry1
        let removidos = vec![0, 1, 2, 3];

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

    #[test]
    fn teste_triangulo_horario() {
        // Elemento CST (Triângulo de deformação constante)
        // A contagem dos nós começa em zero
        // A montegem das funções de forma é horária. Isso dará determinante negativo.
        // O programa deve ser capaz de resolver isso automaticamente
        //
        //    X-|\ (0)-[0,20]
        //      |  \
        //      |    \
        //      |      \
        //      |   1    \ (1)-[10,15] ------> 50000
        //      |       / |
        //      |      /  |
        //      |     /   |
        //      |    /  2 | (2)-[10,5]
        //      |   /    /:
        //      |  /   /  :
        //      | /  /    :
        //      |/ /      V 50000
        //    X-|/ (3)-[0,0]

        let x_coords_elem_1 = vec![0.0, 10.0, 0.0];
        let y_coords_elem_1 = vec![20.0, 15.0, 0.0];
        let global_node_ele_1: Vec<usize> = vec![0, 1, 3];

        let x_coords_elem_2 = vec![10.0, 10.0, 0.0];
        let y_coords_elem_2 = vec![15.0, 5.0, 0.0];
        let global_node_ele_2: Vec<usize> = vec![1, 2, 3];

        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let mut rigidez_local_element_01 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_1,
            y_coords_elem_1,
            global_node_ele_1,
            espessura,
            constitutive_matrix,
        );

        let mut rigidez_local_element_02 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_2,
            y_coords_elem_2,
            global_node_ele_2,
            espessura,
            constitutive_matrix,
        );

        rigidez_local_element_01.assemble_local();
        rigidez_local_element_02.assemble_local();

        println!(
            "Elemento 1 - nós = {:?}",
            rigidez_local_element_01.global_nodes()
        );

        println!(
            "Elemento 2 - nós = {:?}",
            rigidez_local_element_02.global_nodes()
        );

        let elements = vec![
            ([0, 1, 3], rigidez_local_element_01),
            ([1, 2, 3], rigidez_local_element_02),
        ];
        let k_global = rigidez_global::assemble_sparse(4, &elements);

        //F = [Rx0, Ry0, 50000, 0.0, 0.0, -50.000, Rx3, Ry3]
        // Remover DOFs 0 e 3: Rx0, Ry0, Rx3, Ry3
        let removidos = vec![0, 1, 6, 7];

        let k_reduzida =
            gradiente_conjug_jacobi::reduzir_matriz_k_global_csr(&k_global, &removidos);

        let b = vec![50000.0, 0.0, 0.0, -50000.0];
        let x = gradiente_conjug_jacobi::conjugate_gradient_jacobi(&k_reduzida, &b, 1000, 1e-12);
        println!("x = {:?}", x);
        //Avaliando o erro do vetor X
        //a.x = b
        let x_expected = [1.891e-2, -2.727e-2, -2.147e-3, -4.455e-2];

        let error_x = x
            .iter()
            .zip(x_expected.iter())
            .map(|(x_i, x_i_exp)| (x_i - x_i_exp).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro de X = {}", error_x);
        assert!(error_x < 1e-3);
    }

    #[test]
    fn test_triangulo_unico() {
        //
        //                      /| (2) (1500.0, 1500.0)
        //                    /  |
        //                  /    |
        //                /      |
        //              /        |
        //            /          |
        //          /            |
        //        /              |
        //      /                |
        //    /__________________|
        //     |(0)-[0.0,0.0]    | (1) (1500.0, 0.0)
        //     x                 x
        // Dimensões em mm
        //Força distribuída perpendicular à face 1-3 com F = 10 N/mm²
        //Apoios simples em (0) e em (1)
        //
        //Coordenadas do elemtno
        let x_coords_elem_1 = vec![0.0, 1500.0, 1500.0];
        let y_coords_elem_1 = vec![0.0, 0.0, 1500.0];
        let global_node_ele_1: Vec<usize> = vec![0, 1, 2];

        //Propriedades
        let espessura = 2.0; //mm
        let poisson = 0.3;
        let elasticidade = 70000.0; //N/mm²

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let mut rigidez_local_element_01 = rigidez_local_struct::RigidezLocal::new(
            x_coords_elem_1,
            y_coords_elem_1,
            global_node_ele_1,
            espessura,
            constitutive_matrix,
        );
        rigidez_local_element_01.assemble_local();

        let elements = vec![([0, 1, 2], rigidez_local_element_01)];
        let k_global = rigidez_global::assemble_sparse(3, &elements);

        //F = [-15000+Rx0, 15000+Ry0, 0.0, Ry1, -15000, 15000]
        // Remover DOFs 0 e 1 (apenas em y): Rx0, Ry0, Ry1
        let removidos = vec![0, 1, 3];

        let k_reduzida =
            gradiente_conjug_jacobi::reduzir_matriz_k_global_csr(&k_global, &removidos);

        println!("K global:");
        k_global.print_dense();

        println!("K reduzida:");
        k_reduzida.print_dense();

        // neste caso deve ser b reduzida
        let b = vec![0.0, -15000.0, 15000.0];
        let x = gradiente_conjug_jacobi::conjugate_gradient_jacobi(&k_reduzida, &b, 1000, 1e-12);
        println!("x = {:?}", x);
        //Avaliando o erro do vetor X
        //a.x = b
        let x_expected = [-0.27857142857, -0.83571428571, 0.27857142857];

        let error_x = x
            .iter()
            .zip(x_expected.iter())
            .map(|(x_i, x_i_exp)| (x_i - x_i_exp).powi(2))
            .sum::<f64>()
            .sqrt();

        println!("Erro de X = {}", error_x);
        assert!(error_x < 1e-10);
    }
}
