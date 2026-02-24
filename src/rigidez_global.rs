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
fn assemble_sparse(num_nodes: usize, elements: &Vec<([usize; 3], Matrix6)>) -> CsMat<f64> {
    let total_dofs = num_nodes * 2;

    // matriz em formato triplet (COO)
    let mut triplet = TriMat::<f64>::new((total_dofs, total_dofs));

    for (nodes, k_local) in elements {
        // mapear DOFs globais
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

    // converte para CSR (compactado e eficiente)
    triplet.to_csr()
}

#[cfg(test)]
mod tests {
    use crate::rigidez_global;
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
        let k_local = Matrix6::identity();
        // elementos e nós
        // (2)-----(3)
        //  |  \    |
        //  |   \   |
        //  |     \ |
        // (0)-----(1)
        // e1 = [0, 1, 2] (elemento 1)
        // e2 = [1, 3, 2] (elemento 2)
        let elements = vec![([0, 1, 2], k_local), ([1, 3, 2], k_local)];

        let k_global = assemble_sparse(4, &elements);

        println!("Matriz global (CSR):");
        println!("{:?}", k_global);

        println!("\nConvertendo para matriz densa para visualização:\n");

        let dense = k_global.to_dense();
        println!("{}", dense);
    }

    #[test]
    fn test_rigidez_global_esparsa_01() {
        let k_local: Matrix6<f64> = Matrix6::new(
            1764100., -897000., -807300., 358800., -956800., 538200., -897000., 2511600., 538200.,
            -2152800., 358800., -358800., -807300., 538200., 807300., 0., 0., -538200., 358800.,
            -2152800., 0., 2152800., -358800., 0., -956800., 358800., 0., -358800., 956800., 0.,
            538200., -358800., -538200., 0., 0., 358800.,
        );

        // elementos e nós
        // (2)-----(3)
        //  |  \    |
        //  |   \   |
        //  |     \ |
        // (0)-----(1)
        // e1 = [0, 1, 2] (elemento 1)
        // e2 = [1, 3, 2] (elemento 2)

        let elements = vec![([0, 1, 2], k_local), ([1, 3, 2], k_local)];
        let k_global = assemble_sparse(4, &elements);

        println!("Matriz global (CSR):");
        println!("{:?}", k_global);

        println!("\nConvertendo para matriz densa para visualização:\n");

        let dense = k_global.to_dense();
        println!("{}", dense);
    }
}
