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
        let elasticity = 1.0;
        let espessura = 1.0;
        //Elemento 1: Nós [4, 5, 1]
        let x_coords_el_1 = vec![0.0, 5.0, 0.0];
        let y_coords_el_1 = vec![5.0, 5.0, 10.0];

        //Elemento 2 Nós [5, 2, 1]
        let x_coords_el_2 = vec![5.0, 5.0, 0.0];
        let y_coords_el_2 = vec![5.0, 10.0, 10.0];

        //Elemento 3: Nós [5,6,2]
        let x_coords_el_3 = vec![5.0, 10.0, 5.0];
        let y_coords_el_3 = vec![5.0, 5.0, 10.0];

        //Elemento 4: Nós [2,6,3]
        let x_coords_el_4 = vec![5.0, 10.0, 10.0];
        let y_coords_el_4 = vec![10.0, 5.0, 10.0];

        //Elemento 5: Nós [4,8,5]
        let x_coords_el_5 = vec![0.0, 5.0, 5.0];
        let y_coords_el_5 = vec![5.0, 0.0, 5.0];

        //Elemento 6: Nós [4,7,8]
        let x_coords_el_6 = vec![0.0, 0.0, 5.0];
        let y_coords_el_6 = vec![5.0, 0.0, 0.0];

        //Elemento 7: Nós [5,8,9]
        let x_coords_el_7 = vec![5.0, 5.0, 10.0];
        let y_coords_el_7 = vec![5.0, 0.0, 0.0];

        //Elemento 8: Nós [9,6,5]
        let x_coords_el_8 = vec![10.0, 10.0, 5.0];
        let y_coords_el_8 = vec![0.0, 5.0, 5.0];

        //Matriz constitutiva
        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticity);

        //Rigidez do elemento 1
        let rigidez_local_el_1 =
            matriz_rigidez_local(x_coords_el_1, y_coords_el_1, espessura, constitutive_matrix);
        //Rigidez do elemento 2
        let rigidez_local_el_2 =
            matriz_rigidez_local(x_coords_el_2, y_coords_el_2, espessura, constitutive_matrix);
        //Rigidez do elemento 3
        let rigidez_local_el_3 =
            matriz_rigidez_local(x_coords_el_3, y_coords_el_3, espessura, constitutive_matrix);
        //Rigidez do elemento 4
        let rigidez_local_el_4 =
            matriz_rigidez_local(x_coords_el_4, y_coords_el_4, espessura, constitutive_matrix);
        //Rigidez do elemento 5
        let rigidez_local_el_5 =
            matriz_rigidez_local(x_coords_el_5, y_coords_el_5, espessura, constitutive_matrix);
        //Rigidez do elemento 6
        let rigidez_local_el_6 =
            matriz_rigidez_local(x_coords_el_6, y_coords_el_6, espessura, constitutive_matrix);
        //Rigidez do elemento 7
        let rigidez_local_el_7 =
            matriz_rigidez_local(x_coords_el_7, y_coords_el_7, espessura, constitutive_matrix);
        //Rigidez do elemento 8
        let rigidez_local_el_8 =
            matriz_rigidez_local(x_coords_el_8, y_coords_el_8, espessura, constitutive_matrix);

        //Mapeamento de cada nó em cada elemento.
        //O índice i = nó-1
        let elements = vec![
            ([3, 4, 0], rigidez_local_el_1),
            ([4, 1, 0], rigidez_local_el_2),
            ([4, 5, 1], rigidez_local_el_3),
            ([1, 5, 2], rigidez_local_el_4),
            ([3, 7, 4], rigidez_local_el_5),
            ([3, 6, 7], rigidez_local_el_6),
            ([4, 7, 8], rigidez_local_el_7),
            ([8, 5, 4], rigidez_local_el_8),
        ];
    }
}
