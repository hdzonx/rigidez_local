pub mod element{
    pub trait Element{
        type elementType;
        fn new(id: u32,
            global_nodes: Vec<u32>,
            x_coords: Vec<f64>,
            y_coords: Vec<f64>)->Self::elementType;

        fn get_xcoords(&self)->Vec<f64> ;
        fn get_ycoords(&self) -> Vec<f64>;
        fn get_id(&self)->u32;
        fn get_global_nodes(&self)->Vec<u32>;
    }
pub struct Triangle {
        id: u32,
        global_nodes: Vec<u32>,
        x_coords: Vec<f64>,
        y_coords: Vec<f64>,
    }

    impl Element for Triangle {
        type elementType = Triangle;
         fn new(
            id: u32,
            global_nodes: Vec<u32>,
            x_coords: Vec<f64>,
            y_coords: Vec<f64>,
        ) -> Self::elementType {
            Self {
                id,
                global_nodes,
                x_coords,
                y_coords,
            }
        }

        fn get_xcoords(&self)->Vec<f64>  {
          self.x_coords.clone()
           
        }
         fn get_ycoords(&self) -> Vec<f64> {
            self.y_coords.clone()
        }
         fn get_id(&self)->u32{
            self.id.clone()
        }
        fn get_global_nodes(&self)->Vec<u32>{
            self.global_nodes.clone()
        }
    }

}
