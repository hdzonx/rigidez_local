pub struct Triangle {
    id: String,
    global_nodes: Vec<String>,
    x_coords: Vec<f64>,
    y_coords: Vec<f64>,
}

impl Triangle {
    pub fn new(
        id: String,
        global_nodes: Vec<String>,
        x_coords: Vec<f64>,
        y_coords: Vec<f64>,
    ) -> Self {
        Self {
            id,
            global_nodes,
            x_coords,
            y_coords,
        }
    }

    pub fn get_xcoords(&self) -> Vec<f64> {
        self.x_coords.clone()
    }
    pub fn get_ycoords(&self) -> Vec<f64> {
        self.y_coords.clone()
    }
}
