use crate::gradiente_conjug_jacobi::CsrMatrix;
use crate::rigidez_local_struct::RigidezLocal;
use nalgebra::{DMatrix, SMatrix};
use sprs::{CsMat, TriMat};

type Matrix6 = SMatrix<f64, 6, 6>;

//Não eficiente para matriz esparsa
pub fn assemble_global_triangle(
    global_size: usize,
    locals: &Vec<(Matrix6, Vec<usize>)>,
) -> DMatrix<f64> {
    let mut k_global = DMatrix::<f64>::zeros(global_size, global_size);

    for (k_local, mapping) in locals {
        for i in 0..6 {
            for j in 0..6 {
                let I = mapping[i];
                let J = mapping[j];

                k_global[(I, J)] += k_local[(i, j)];
            }
        }
    }

    k_global
}
//Eficiente para matriz esparsa
pub fn assemble_sparse(num_nodes: usize, elements: &Vec<([usize; 3], RigidezLocal)>) -> CsrMatrix {
    let total_dofs = num_nodes * 2;

    let mut triplet = TriMat::<f64>::new((total_dofs, total_dofs));

    for (_, rigidez_local) in elements {
        let nodes = rigidez_local.global_nodes();
        let k_local = &rigidez_local.matrix;

        println!("Montando elemento com nós = {:?}", nodes);

        let dofs = [
            2 * nodes[0],
            2 * nodes[0] + 1,
            2 * nodes[1],
            2 * nodes[1] + 1,
            2 * nodes[2],
            2 * nodes[2] + 1,
        ];

        for i in 0..6 {
            for j in 0..6 {
                triplet.add_triplet(dofs[i], dofs[j], k_local[(i, j)]);
            }
        }
    }

    let csr = triplet.to_csr();

    CsrMatrix {
        values: csr.data().to_vec(),
        col_indices: csr.indices().to_vec(),
        row_ptr: csr.indptr().raw_storage().to_vec(),
        n: total_dofs,
    }
}

#[cfg(test)]
mod tests {
    use crate::constitutive_matrix;
    use crate::rigidez_global;
    use crate::rigidez_local_struct;
    use nalgebra::Matrix6;
    #[test]
    fn test_rigidez_global() {
        let k1 = Matrix6::identity();
        let map1 = vec![0, 1, 2, 3, 4, 5];

        let k2 = Matrix6::from_element(2.0);
        let map2 = vec![3, 4, 5, 6, 7, 8];

        let locals = vec![(k1, map1), (k2, map2)];

        let k_global = rigidez_global::assemble_global_triangle(9, &locals);

        println!("{}", k_global);
    }

    use crate::rigidez_global::assemble_sparse;
    #[test]
    fn test_rigidez_global_esparsa() {
        // matriz local simples para teste

        //Matriz constitutiva
        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;
        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        //Dados do elemento 1
        let x_coords_ele_01 = vec![0.0, 10.0, 10.0];
        let y_coords_ele_01 = vec![0.0, 5.0, 15.0];
        let global_node_ele_1: Vec<usize> = vec![0, 1, 2];
        let k_local_1 = rigidez_local_struct::RigidezLocal::new(
            x_coords_ele_01,
            y_coords_ele_01,
            global_node_ele_1,
            espessura,
            constitutive_matrix,
        );

        //Dados do elemento 2
        let x_coords_el_02 = vec![0.0, 10.0, 0.0];
        let y_coords_el_02 = vec![0.0, 15.0, 20.0];
        let global_node_ele_2: Vec<usize> = vec![0, 2, 3];

        let k_local_2 = rigidez_local_struct::RigidezLocal::new(
            x_coords_el_02,
            y_coords_el_02,
            global_node_ele_2,
            espessura,
            constitutive_matrix,
        );
        // elementos e nós
        // (2)-----(3)
        //  |  \    |
        //  |   \   |
        //  |     \ |
        // (0)-----(1)
        // e1 = [0, 1, 2] (elemento 1)
        // e2 = [1, 3, 2] (elemento 2)

        let elements = vec![([0, 1, 2], k_local_1), ([0, 2, 3], k_local_2)];

        let k_global = assemble_sparse(4, &elements);

        println!("Matriz global (CSR):");
        println!("{:?}", k_global);

        println!("\nConvertendo para matriz densa para visualização:\n");

        // let dense = k_global.to_dense();
        //  println!("{}", dense);
    }

    // #[test]
    // fn test_rigidez_global_esparsa_01() {
    //     let k_local: Matrix6<f64> = Matrix6::new(
    //         1764100., -897000., -807300., 358800., -956800., 538200., -897000., 2511600., 538200.,
    //         -2152800., 358800., -358800., -807300., 538200., 807300., 0., 0., -538200., 358800.,
    //         -2152800., 0., 2152800., -358800., 0., -956800., 358800., 0., -358800., 956800., 0.,
    //         538200., -358800., -538200., 0., 0., 358800.,
    //     );

    //     // elementos e nós
    //     // (2)-----(1)
    //     //  |     / |
    //     //  |   /   |
    //     //  | /     |
    //     // (3)-----(0)
    //     // e1 = [0, 1, 3] (elemento 1)
    //     // e2 = [2, 3, 1] (elemento 2)

    //     let elements = vec![([0, 1, 3], k_local), ([2, 3, 1], k_local)];
    //     let k_global = assemble_sparse(4, &elements);

    //     println!("Matriz global (CSR):");
    //     println!("{:?}", k_global);

    //     println!("\nConvertendo para matriz densa para visualização:\n");

    //    // let dense = k_global.to_dense();
    //    // println!("{}", dense);
    // }

    // #[test]
    // fn nova_definicao_assemble_sparse(){

    // }
}
