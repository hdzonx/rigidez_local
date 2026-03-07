use nalgebra::Matrix3;

fn constitutive_matrix(plane_state: &str, poisson: f64, elasticity: f64) {
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
}

#[cfg(test)]
mod tests {
    use crate::constitutive_matrix::constitutive_matrix;

    #[test]
    fn const_matr_test_01() {
        let cm = constitutive_matrix("Plane train", 0.3, 1000000.0);
    }
}
