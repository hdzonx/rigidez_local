use crate:: element::element::{Triangle, Element};



fn surface_force(elem:&Triangle, global_node_force_act:Vec<u32>, thickness:f64){

    let coord_x = elem.get_xcoords();
    let corrd_y = elem.get_ycoords();

    let global_node_element= elem.get_global_nodes();

    let node_1 = &global_node_force_act[0];
    let node_2 = &global_node_force_act[1];

    let x1:f64;
    let x2:f64;
    let y1:f64;
    let y2:f64;

    let mut counter =0;

    for i in &global_node_element{
        println!("{}", i);
        if node_1==i{
            println!("nó = {} ", node_1);

        }
        counter+=1;

    }

 //   let lenght_side = 

}

#[cfg(test)]
mod tests {
use crate::{ element::element::{Element, Triangle}, forces::surface_force};

    
    #[test]
    fn nodal_force_test_01(){   
            let elem_01 = Triangle::new(1, vec![3,4,7], vec![75.0, 0.0, 75.0], vec![0.0, 0.0, 50.0]);
            let global_node_force_act = vec![3,4];
            let esp = 5.0;

            surface_force(&elem_01, global_node_force_act, esp);


    }
}