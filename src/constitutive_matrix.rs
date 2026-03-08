use nalgebra::Matrix3;
#[allow(unused)]

fn constitutive_matrix(plane_state: &str, poisson: f64, elasticity: f64)->Matrix3<f64> {
    let mut constitutive_matrix = Matrix3::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    if plane_state.eq_ignore_ascii_case("Plane Strain") {
        println!("Plane Strain");
        constitutive_matrix = Matrix3::new(
            elasticity * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
            elasticity * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
            0.,
            elasticity * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
            elasticity * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
            0.,
            0.,
            0.,
            elasticity * (1.0 - 2.0 * poisson) / (2.0 * (1.0 + poisson) * (1.0 - 2.0 * poisson)),
        );
    } else if plane_state.eq_ignore_ascii_case("Plane Stress") {
        println!("Plane Stress");
        constitutive_matrix = Matrix3::new(
            elasticity / (1.0 - poisson * poisson),
            poisson * elasticity / (1.0 - poisson * poisson),
            0.,
            poisson * elasticity / (1.0 - poisson * poisson),
            elasticity / (1.0 - poisson * poisson),
            0.,
            0.,
            0.,
            elasticity / (2.0 * (1.0 + poisson)),
        );
    } else {
        panic!(" Plane states must be 'Plane Strain' or 'Plane Stress'.");
    }
    constitutive_matrix
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Matrix3;

    use crate::constitutive_matrix::constitutive_matrix;

    #[test]
    fn plane_stress_test_1() {
        let cm = constitutive_matrix("Plane stress", 0.3, 30.0E6);

        let expected = Matrix3::new(
            3.297E7, 0.9891E7, 0.0,
            0.9891E7, 3.297E7, 0.0,
            0.0, 0.0, 1.15395E7);

            for i in 0..expected.len(){
                println!("{}",cm[i]);

                assert_relative_eq!(cm[i], expected[i], epsilon=0.001E7);
            }

        println!("const. matrx = {:?}", cm);
    }
}
