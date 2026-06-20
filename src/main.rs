use std::vec;
mod constitutive_matrix;
mod element;
//mod rigidez_global_obsoleta;
mod global_force;
mod gradiente_conjug_jacobi;
mod gradiente_conjugado;
mod gradiente_conjugado_cholesky;
mod local_forces;
mod matriz_reduzida;
mod matriz_reduzida_eficiente;
mod rigidez_global;
mod rigidez_local;
mod rigidez_local_struct;
mod test1;
mod test2;
mod tests;

use element::element::{Element, Triangle};

use nalgebra::{Const, Matrix3, Matrix3x6, Matrix6};
use nalgebra::{DMatrix, SMatrix};

fn main() {
    let espessura = 13.0;
    let constitutive_matrix =
        Matrix3::new(220800., 55200., 0., 55200., 220800., 0., 0., 0., 82800.);

    let triangulo_1 = Triangle::new(
        1,
        vec![1, 2, 4],
        vec![75.0, 0.0, 75.0],
        vec![0.0, 0.0, 50.0],
    );

    let triangulo_2 = Triangle::new(
        2,
        vec![3, 4, 2],
        vec![0.0, 0.0, 75.0],
        vec![50.0, 0.0, 50.0],
    );
}
