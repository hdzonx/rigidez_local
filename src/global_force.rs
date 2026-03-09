


fn assemble_global_forces(){
    
}

#[cfg(test)]
mod tests {
    use crate::element::element::{Element, Triangle};
    #[test]
    fn global_force_test_01() {
        //número de elementos
        let ne = 2;
        //nós globais onde agem as forças
        let global_node_1 = vec![8, 9];
        let global_node_2 = vec![9, 10];
        //vetores de forças
        let force_1 = vec![2.0, 3.0, 4.0, 5.0];
        let force_2 = vec![5.0, 6.0, 7.0, 8.0];

        let locals = vec![(force_1, global_node_1), (force_2, global_node_2)];
    }
}
