use nalgebra::{Const, Matrix3, Matrix3x6, Matrix6};
use std::collections::HashMap;
pub struct RigidezLocal {
    x_coords: Vec<f64>,
    y_coords: Vec<f64>,
    global_node_num: Vec<usize>,
    espessura: f64,
    constitutive_matrix: Matrix3<f64>,
    pub matrix: Matrix6<f64>,
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
            matrix: Matrix6::zeros(),
        }
    }

    pub fn assemble_local(&mut self) -> Matrix6<f64> {
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

        let mut global_node_num_update: Vec<usize> = self.global_node_num.clone();

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
            global_node_num_update = vec![
                self.global_node_num[0],
                self.global_node_num[2],
                self.global_node_num[1],
            ];

            det_jacobian = x13 * y23 - x23 * y13;
        }
        self.global_node_num = global_node_num_update;

        let x32 = x3 - x2;
        let x21 = x2 - x1;
        let y31 = y3 - y1;
        let y12 = y1 - y2;

        println!("global nodes update = {:?}", self.global_node_num);
        println!("determinante de J = {}", det_jacobian);

        let m: nalgebra::Matrix<f64, Const<3>, Const<6>, nalgebra::ArrayStorage<f64, 3, 6>> =
            Matrix3x6::new(
                y23, 0.0, y31, 0.0, y12, 0.0, 0.0, x32, 0.0, x13, 0.0, x21, x32, y23, x13, y31,
                x21, y12,
            );

        let b: nalgebra::Matrix<f64, Const<3>, Const<6>, nalgebra::ArrayStorage<f64, 3, 6>> =
            (1.0 / (det_jacobian)) * m;
        println!("b = {}", b);
        let b_transposta = &b.transpose();

        //Matriz constitutiva para estado plano de tensões

        let db = self.constitutive_matrix * b;
        println!("dxb = {}", db);

        let area = area_triangulo(&self.x_coords, &self.y_coords);

        // Matriz de rigidez local
        let k = area * self.espessura * b_transposta * db;

        println!("k = {}", k);

        self.matrix = k;

        self.matrix
    }

    pub fn global_nodes(&self) -> &[usize] {
        &self.global_node_num
    }
}

pub fn area_triangulo(x_coords: &Vec<f64>, y_coords: &Vec<f64>) -> f64 {
    if x_coords.len() != 3 || y_coords.len() != 3 {
        panic!("dimension of coordinates vector must be 3");
    }
    let m = Matrix3::new(
        x_coords[0],
        y_coords[0],
        1.0,
        x_coords[1],
        y_coords[1],
        1.0,
        x_coords[2],
        y_coords[2],
        1.0,
    );

    let det = m.determinant();
    let area = 0.5 * det;

    //A área não poderá ser negativa. A sentença
    //abaixo garante isso.
    if area < 0.0 {
        return -area;
    }

    area
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
