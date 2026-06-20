use nalgebra::{Const, Matrix3, Matrix3x6, Matrix6};
use std::collections::HashMap;
pub struct RigidezLocal {
    x_coords: Vec<f64>,
    y_coords: Vec<f64>,
    global_node_num: Vec<usize>,
    espessura: f64,
    constitutive_matrix: Matrix3<f64>,
}

impl RigidezLocal {
    pub fn new(
        x_coords: Vec<f64>,
        y_coords: Vec<f64>,
        global_node_num: Vec<usize>,
        espessura: f64,
        constitutive_matrix: Matrix3<f64>,
    ) -> RigidezLocal {
        RigidezLocal {
            x_coords,
            y_coords,
            global_node_num,
            espessura,
            constitutive_matrix,
        }
    }

    fn assemble_local(&mut self) {
        let x1 = self.x_coords[0];
        let mut x2 = self.x_coords[1];
        let mut x3 = self.x_coords[2];

        let y1 = self.y_coords[0];
        let mut y2 = self.y_coords[1];
        let mut y3 = self.y_coords[2];

        let mut y23 = y2 - y3;
        let mut y13 = y1 - y3;

        let mut x13 = x1 - x3;
        let mut x23 = x2 - x3;

        let mut det_jacobian = x13 * y23 - x23 * y13;

        let mut global_node_num_updated: Vec<usize> = self.global_node_num.clone();

        if det_jacobian < 0.0 {
            // se o determinante do jacobiano for menor que zero,
            //significa que o giro no triângulo será horário. Neste caso,
            //devemos inverter a ordem das coordenadas locais 2 e 3.
            //As coordenadas globais não são alteradas.
            x2 = self.x_coords[2];
            x3 = self.x_coords[1];

            y2 = self.y_coords[2];
            y3 = self.y_coords[1];

            y23 = y2 - y3;
            y13 = y1 - y3;
            x13 = x1 - x3;
            x23 = x2 - x3;
            global_node_num_updated = vec![
                self.global_node_num[0],
                self.global_node_num[2],
                self.global_node_num[1],
            ];

            det_jacobian = x13 * y23 - x23 * y13;
        }

        let x32 = x3 - x2;
        let x21 = x2 - x1;
        let y31 = y3 - y1;
        let y12 = y1 - y2;

        //let global_num = self.global_nodes_num(&global_node_num_updated);
        self.global_node_num = global_node_num_updated;

        println!("global nodes update = {:?}", self.global_node_num);
        println!("determinante de J = {}", det_jacobian);
    }

    fn global_nodes_num<'a>(&self, global_nodes_num: &'a Vec<usize>) -> &'a Vec<usize> {
        &global_nodes_num
    }
}

#[cfg(test)]
mod tests {
    use crate::rigidez_local_struct;

    use super::*;
    use crate::constitutive_matrix;

    #[test]
    fn local_test_01() {
        let x_coords = vec![75.0, 0.0, 75.0];
        let y_coords = vec![0.0, 0.0, 50.0];
        let global_node_num: Vec<usize> = vec![0, 1, 3];
        let espessura = 1.0;
        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let mut rigidez_1 = rigidez_local_struct::RigidezLocal::new(
            x_coords,
            y_coords,
            global_node_num,
            espessura,
            constitutive_matrix,
        );
        rigidez_1.assemble_local();
        let calc_global_num = rigidez_1.global_node_num;
        let real_global_num: Vec<usize> = vec![0, 3, 1];
        for i in 0..3 {
            assert_eq!(calc_global_num[i], real_global_num[i]);
        }
    }
    #[test]
    #[should_panic]
    fn local_test_02() {
        let x_coords = vec![75.0, 0.0, 75.0];
        let y_coords = vec![0.0, 0.0, 50.0];
        let global_node_num: Vec<usize> = vec![0, 1, 3];
        let espessura = 1.0;
        let espessura = 0.1;
        let poisson = 0.3;
        let elasticidade = 30.0e6;

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let mut rigidez_1 = rigidez_local_struct::RigidezLocal::new(
            x_coords,
            y_coords,
            global_node_num,
            espessura,
            constitutive_matrix,
        );
        rigidez_1.assemble_local();
        let calc_global_num = rigidez_1.global_node_num;
        let real_global_num: Vec<usize> = vec![0, 1, 3];
        for i in 0..3 {
            if calc_global_num[i] != real_global_num[i] {
                panic!("Erro esperado.")
            }
        }
    }
}
