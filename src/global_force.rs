fn assemble_global_force(global_size: usize, f_locals: &Vec<(Vec<f64>, Vec<usize>)>) -> Vec<f64> {
    let mut f_global = vec![0.0; global_size];
    for (forces, mapping) in f_locals {
        for i in 0..4 {
            let i_global = mapping[i];
            f_global[i_global] += forces[i];
        }
    }
    f_global
}

#[cfg(test)]
mod tests {
    use crate::element::element::{Element, Triangle};
    use crate::global_force::assemble_global_force;
    #[test]
    fn global_force_test_01() {
        //nós globais onde agem as forças (mapeamento)
        let global_node_1 = vec![0, 1, 2, 3];
        let global_node_2 = vec![2, 3, 4, 5];
        //vetores de forças
        let force_1 = vec![2.0, 3.0, 4.0, 5.0];
        let force_2 = vec![5.0, 6.0, 7.0, 8.0];

        let f_locals: Vec<(Vec<f64>, Vec<usize>)> =
            vec![(force_1, global_node_1), (force_2, global_node_2)];

        let expected = vec![2.0, 3.0, 9.0, 11.0, 7.0, 8.0];
        let result = assemble_global_force(6, &f_locals);

        for i in 0..expected.len() {
            approx::assert_abs_diff_eq!(result[i], expected[i], epsilon = 0.001);
        }
        println!("result = {:?}", result);
    }
    #[test]
    fn global_force_test_02() {
        //nós globais onde agem as forças (mapeamento)
        let global_node_1 = vec![0, 1, 2, 3];
        let global_node_2 = vec![0, 1, 4, 5];
        //vetores de forças
        let force_1 = vec![2.0, 3.0, 4.0, 5.0];
        let force_2 = vec![5.0, 6.0, 7.0, 8.0];

        let f_locals: Vec<(Vec<f64>, Vec<usize>)> =
            vec![(force_1, global_node_1), (force_2, global_node_2)];

        let expected = vec![7.0, 9.0, 4.0, 5.0, 7.0, 8.0];
        let result = assemble_global_force(6, &f_locals);

        for i in 0..expected.len() {
            approx::assert_abs_diff_eq!(result[i], expected[i], epsilon = 0.001);
        }
        println!("result = {:?}", result);
    }
}
