use crate::element::element::{Element, Triangle};

fn surface_force(elem: &Triangle, global_node_force_act: Vec<u32>, thickness: f64) {
    let coord_x = elem.get_xcoords();
    let coord_y = elem.get_ycoords();

    let global_node_element = elem.get_global_nodes();

    let node_1 = &global_node_force_act[0];
    let node_2 = &global_node_force_act[1];

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut y1: f64 = 0.0;
    let mut y2: f64 = 0.0;

    if let Some(pos) = global_node_element.iter().position(|x| x == node_1) {
        println!("Encontrado no índice: {}", pos);
        x1 = coord_x[pos];
        y1 = coord_y[pos];
    }

    if let Some(pos) = global_node_element.iter().position(|x| x == node_2) {
        println!("Encontrado no índice: {}", pos);
        x2 = coord_x[pos];
        y2 = coord_y[pos];
    }

    let lenght_size = ((x1 - x2).powf(2.0) + (y1 - y2).powf(2.0)).sqrt();
    if lenght_size == 0.0 {
        panic!(
            "value of lenght zero encontred in element = {}",
            elem.get_id()
        );
    }
    println!("lenght size of load surface = {}", lenght_size);
}

#[cfg(test)]
mod tests {
    use crate::{
        element::element::{Element, Triangle},
        forces::surface_force,
    };

    #[test]
    fn nodal_force_test_01() {
        let elem_01 = Triangle::new(
            1,
            vec![3, 4, 7],
            vec![75.0, 0.0, 75.0],
            vec![0.0, 0.0, 50.0],
        );
        let global_node_force_act = vec![3, 4];
        let esp = 5.0;

        surface_force(&elem_01, global_node_force_act, esp);
    }
}
