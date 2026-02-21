use std::vec;

use nalgebra::{Const, Matrix3, Matrix3x6};

fn main() {
    let x_coords = vec![75.0, 0.0, 75.0];
    let y_coords = vec![0.0, 0.0, 50.0];
    let espessura = 13.0;
    matriz_rigidez_local(x_coords,y_coords, espessura);
}

fn matriz_rigidez_local(x_coords: Vec<f64>, y_coords: Vec<f64>, espessura:f64) {
    let x1 = x_coords[0];
    let mut x2 = x_coords[1];
    let mut x3 = x_coords[2];

    let y1 = y_coords[0];
    let mut y2 = y_coords[1];
    let mut y3 = y_coords[2];

    let mut y23 = y2 - y3;
    let mut y13 = y1 - y3;

    let mut x13 = x1 - x3;
    let mut x23 = x2 - x3;

    let mut det_j = x13 * y23 - x23 * y13;

    if det_j < 0.0 {
        // se o determinante do jacobiano for menor que zero,
        //significa que o giro no triângulo será horário. Neste caso,
        //devemos inverter a ordem das coordenadas locais 2 e 3.
        //As coordenadas globais não são alteradas.
        x2 = x_coords[2];
        x3 = x_coords[1];

        y2 = y_coords[2];
        y3 = y_coords[1];

        y23 = y2 - y3;
        y13 = y1 - y3;
        x13 = x1 - x3;
        x23 = x2 - x3;

        det_j = x13 * y23 - x23 * y13;
    }

    let x32 = x3 - x2;
    let x21 = x2 - x1;
    let y31 = y3 - y1;
    let y12 = y1 - y2;

    println!("determinante de J = {}", det_j);

    let m: nalgebra::Matrix<f64, Const<3>, Const<6>, nalgebra::ArrayStorage<f64, 3, 6>> =
        Matrix3x6::new(
            y23, 0.0, y31, 0.0, y12, 0.0, 0.0, x32, 0.0, x13, 0.0, x21, x32, y23, x13, y31, x21,
            y12,
        );

    let b: nalgebra::Matrix<f64, Const<3>, Const<6>, nalgebra::ArrayStorage<f64, 3, 6>> =
        (1.0 / (det_j)) * m;
    println!("b = {}", b);
    let b_transposta = &b.transpose();

    //Matriz constitutiva para estado plano de tensões
    let d = Matrix3::new(220800., 55200., 0., 55200., 220800., 0., 0., 0., 82800.);

    let db = d * b;
    println!("dxb = {}", db);

    let area = area_triangulo(x_coords, y_coords);

    //Matriz de rigidez local
    let k =area*espessura* b_transposta * db;
    println!("k = {}", k);
}

fn area_triangulo(x_coords: Vec<f64>, y_coords: Vec<f64>) -> f64 {
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
    use approx::assert_relative_eq;
    use nalgebra::Matrix3;

    use crate::area_triangulo;

    #[test]
    fn calcula_triangulo() {
        let x_coords = vec![75.0, 0.0, 75.0];
        let y_coords = vec![0.0, 0.0, 50.0];

        let calculado = area_triangulo(x_coords, y_coords);

        // àrea esperada
        let esperado = 1875.0;

        assert_relative_eq!(calculado, esperado, epsilon = 1e-10);
    }
}
