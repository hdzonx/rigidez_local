//Elemento CST
//  [/]1 ________2______3   --------> 100
//      |\   2  |\   4  |
//      |  \    |  \    |
//      | 1  \  | 3  \  |
//  [/]4|_____\5|_____\6|  ----------> 200
//      |\   5  |\   8  |
//      |  \    |  \    |
//      | 6  \  | 7  \  |
//  [/]7|_____\8|_____\9| ----------> 100
//

use sprs::{CsMat, TriMat};
#[cfg(test)]
mod test {

    use approx::assert_relative_eq;
    //use approx::relative_eq;
    use nalgebra::DMatrix;
    use nalgebra::Matrix3;
    use nalgebra::Matrix6;
    use nalgebra::SMatrix;

    use crate::area_triangulo;
    use crate::constitutive_matrix;
    use crate::gradiente_conjug_jacobi;
    use crate::matriz_reduzida_eficiente;
    use crate::matriz_rigidez_local;
    use crate::rigidez_global;

    use sprs::{CsMat, TriMat};

    #[test]
    fn chapa() {
        //Propriedades de material
        let poisson = 0.3;
        let elasticity = 30.0e6;
        let espessura = 1.0;
        //Elemento 1
        let x_coords_el_1 = vec![0.0, 0.0, 5.0];
        let y_coords_el_1 = vec![10.0, 5.0, 5.0];

        //Elemento 2
        let x_coords_el_2 = vec![0.0, 5.0, 5.0];
        let y_coords_el_2 = vec![10.0, 5.0, 10.0];

        //Elemento 3
        let x_coords_el_3 = vec![5.0, 5.0, 10.0];
        let y_coords_el_3 = vec![10.0, 5.0, 5.0];

        //Elemento 4
        let x_coords_el_4 = vec![5.0, 10.0, 10.0];
        let y_coords_el_4 = vec![10.0, 5.0, 10.0];

        //Elemento 5
        let x_coords_el_5 = vec![0.0, 5.0, 5.0];
        let y_coords_el_5 = vec![5.0, 0.0, 5.0];

        //Elemento 6
        let x_coords_el_6 = vec![0.0, 0.0, 5.0];
        let y_coords_el_6 = vec![5.0, 0.0, 0.0];

        //Elemento 7
        let x_coords_el_7 = vec![5.0, 5.0, 10.0];
        let y_coords_el_7 = vec![5.0, 0.0, 0.0];

        //Elemento 8
        let x_coords_el_8 = vec![10.0, 10.0, 5.0];
        let y_coords_el_8 = vec![0.0, 5.0, 5.0];

        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticity);

        let rigidez_local_el_1 =
            matriz_rigidez_local(x_coords_el_1, y_coords_el_1, espessura, constitutive_matrix);

        println!("rigidez do elemento 1: {:?}", rigidez_local_el_1);
    }
}
