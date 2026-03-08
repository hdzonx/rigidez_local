use crate::element::element::{Element, Triangle};
#[allow(unused)]
fn surface_force(
    elem: &Triangle,
    global_node_force_act: Vec<u32>,
    thickness: f64,
    act_force_val: Vec<f64>,
) -> Vec<f64> {
    // O vetor de forças atua considerando apenas dois nós,a e b
    let node_a = &global_node_force_act[0];
    let node_b = &global_node_force_act[1];

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut y1: f64 = 0.0;
    let mut y2: f64 = 0.0;

    //Mapeia as coordenadas nos nós locais node_a e node_b
    if let Some(pos) = elem.get_global_nodes().iter().position(|x| x == node_a) {
        x1 = elem.get_xcoords()[pos];
        y1 = elem.get_ycoords()[pos];
    }

    if let Some(pos) = elem.get_global_nodes().iter().position(|x| x == node_b) {
        x2 = elem.get_xcoords()[pos];
        y2 = elem.get_ycoords()[pos];
    }

    let lenght_size = ((x1 - x2).powf(2.0) + (y1 - y2).powf(2.0)).sqrt();
    if lenght_size == 0.0 {
        panic!(
            "value of lenght zero encontred in element = {}",
            elem.get_id()
        );
    }

    let force_xa = -(y2 - y1) / lenght_size * act_force_val[0];
    let force_ya = -(x1 - x2) / lenght_size * act_force_val[0];

    let force_xb = -(y2 - y1) / lenght_size * act_force_val[1];
    let force_yb = -(x1 - x2) / lenght_size * act_force_val[1];

    let scalar = thickness * lenght_size / 6.0;

    let pre_force = vec![
        2.0 * force_xa + force_xb,
        2.0 * force_ya + force_yb,
        force_xa + 2.0 * force_xb,
        force_ya + 2.0 * force_yb,
    ];

    let force_vec: Vec<f64> = pre_force.iter().map(|x| x * scalar).collect();
    force_vec
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
            vec![8, 9, 10],
            vec![85.0, 70.0, 50.0],
            vec![40.0, 60.0, 40.0],
        );
        let global_node_force_act = vec![8, 9];
        let esp = 10.0;
        let force_val = vec![2.0, 3.0];

        let calc_force_vec = surface_force(&elem_01, global_node_force_act, esp, force_val);
        println!("force vecto = {:?}", calc_force_vec);

        let expected_forces = vec![-233.33, -175.0, -266.67, -200.0];

        for i in 0..calc_force_vec.len() {
            approx::assert_abs_diff_eq!(calc_force_vec[i], expected_forces[i], epsilon = 0.1);
        }
    }
    #[test]
    fn nodal_force_test_02() {
        let elem_01 = Triangle::new(
            1,
            vec![7, 8, 9],
            vec![100.0, 85.0, 50.0],
            vec![20.0, 40.0, 20.0],
        );
        let global_node_force_act = vec![7, 8];
        let esp = 10.0;
        let force_val = vec![1.0, 2.0];

        let calc_force_vec = surface_force(&elem_01, global_node_force_act, esp, force_val);
        println!("force vecto = {:?}", calc_force_vec);

        let expected_forces = vec![-133.33, -100.0, -166.67, -125.0];

        for i in 0..calc_force_vec.len() {
            approx::assert_abs_diff_eq!(calc_force_vec[i], expected_forces[i], epsilon = 0.1);
        }


    }
}
