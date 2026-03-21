// Elemento CST (Triângulo de deformação constante)
//
//      |\ (4)-[0,20]
//      |  \
//      |    \
//      |      \
//      |        \ (3)-[10,15]
//      |       / |
//      |      /  |
//      |     /   |
//      |    /    | (2)-[10,5]
//      |   /    /
//      |  /   /
//      | /  /
//      |/ /
//      |/ (1)-[0,0]

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    //use approx::relative_eq;
    use nalgebra::Matrix3;
    use nalgebra::Matrix6;

    use crate::area_triangulo;
    use crate::constitutive_matrix;
    use crate::matriz_rigidez_local;

    #[test]
    fn calcula_rigidez_local_el_1() {
        let x_coords = vec![0.0, 10.0, 10.0];
        let y_coords = vec![0.0, 5.0, 15.0];

        let area = area_triangulo(&x_coords, &y_coords);
        let espessura = 0.1;

        let escalar = 3.297e6;

        let matriz_inicial: Matrix6<f64> = Matrix6::new(
            0.5, 0.0, -0.75, 0.15, 0.25, -0.15, 0., 0.175, 0.175, -0.263, -0.175, 0.088, -0.75,
            0.1750, 1.3, -0.488, -0.55, 0.313, 0.15, -0.263, -0.488, 0.894, 0.338, -0.631, 0.25,
            -0.175, -0.55, 0.338, 0.3, -0.163, -1.5, 0.088, 0.313, -0.631, -0.163, 0.544,
        );

        let rigidez_local_esperada = escalar * matriz_inicial;
        //println!("{:?}", rigidez_local_esperada);

        let poisson = 0.3;
        let elasticidade = 30.0e6;
        let constitutive_matrix =
            constitutive_matrix::constitutive_matrix("Plane stress", poisson, elasticidade);

        let rigidez_local_calc =
            matriz_rigidez_local(x_coords, y_coords, espessura, constitutive_matrix);

        println!("{:?}", rigidez_local_calc);
    }
}
