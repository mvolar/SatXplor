use bio;
use plotters::prelude::*;
use bio::io::fasta::Record;
use ndarray::Axis;


//kmer in array plot
pub fn plot_kmer_in_array_raw(record: &Record, kmer_in_array_pos_dist: Vec<i32> ) {

    let output_folder = "./output_folder/";
    let filename = output_folder.to_owned() + &record.id().to_string().to_owned() + "_bool.png";


    let root = BitMapBackend::new(&filename, (800, 600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    // Create a chart context
    let mut chart = ChartBuilder::on(&root)
        .caption("Boolean Vector Plot", ("sans-serif", 40))
        .caption("Line Graph Example", ("Arial", 20))
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0..kmer_in_array_pos_dist.len() as i32, 0..30)
        .unwrap();

    chart
        .configure_mesh()
        .draw()
        .unwrap();
    // Plot the boolean vector
    chart.draw_series(LineSeries::new(
        (0..kmer_in_array_pos_dist.len() as i32).zip(kmer_in_array_pos_dist.clone()),
        &BLACK,
    ))
    .unwrap();
}



        //row sum plot
pub fn plot_row_sum_record(record: &Record, dist_mat:   ndarray::prelude::ArrayBase<ndarray::OwnedRepr<u32>, ndarray::prelude::Dim<[usize; 2]>> ) {
    let row_sums: Vec<f64> = dist_mat.axis_iter(Axis(1))
    .map(|row| row.iter().filter(|&&element| element < 5).count() as f64)
    .collect();
    //
    let output_folder = "./output_folder/";


    let filename = output_folder.to_owned() + &record.id().to_string().to_owned() + "_rowsum_plot_.png";
    println!("{:?}",filename);
    let root = BitMapBackend::new(&filename, (740, 780)).into_drawing_area();
    root.fill(&WHITE)
    .unwrap();

    let mut chart = ChartBuilder::on(&root)
    .caption("Line Graph Example", ("Arial", 20))
    .x_label_area_size(30)
    .y_label_area_size(60)
    .build_cartesian_2d(0.0..row_sums.len() as f64, 0.0..30.5)
    .unwrap();

    chart
    .configure_mesh()
    .draw()
    .unwrap();

    chart.draw_series(
    LineSeries::new(row_sums.clone().into_iter().enumerate().map(|(x, y)| (x as f64, y as f64)), &RED))
    .unwrap();
    
}
    

    //roll sum plot
pub fn plot_roll_mean(record: &Record, roll_mean: Vec<f64>,outpath: &str) {
        
        let filename = outpath.to_owned()+ "/pictures/" + &record.id().to_string().to_owned() + "_rol_sum.png";
    
        
        let root = BitMapBackend::new(&filename, (800, 600)).into_drawing_area();
        root.fill(&WHITE)
        .unwrap();
    
        // Create a chart context
        let mut chart = ChartBuilder::on(&root)
            .caption("Boolean Vector Plot", ("sans-serif", 40))
            .caption("Line Graph Example", ("Arial", 20))
            .x_label_area_size(30)
            .y_label_area_size(60)
            .build_cartesian_2d(0..roll_mean.len() as i32, 0.0..30.0)
            .unwrap();
    
        chart
            .configure_mesh()
            .draw()
            .unwrap();
        // Plot the boolean vector
        chart.draw_series(LineSeries::new(
            (0..roll_mean.len() as i32).zip(roll_mean.clone()),
            &BLACK,
        ))
        .unwrap();
        }
    
    
    
    