use bio;

use std::{error::Error, str::FromStr};
use std::fs::File;
use bio::io::fasta;
use std::fs;
use debruijn::*;
use debruijn::kmer::*;
use ndarray::{ Array2,ArrayBase,OwnedRepr,Dim};
mod utils;
mod ploting;
use ploting::plot_roll_mean;
use std::io::Write;
use rayon::prelude::*;
use log::info;
use env_logger;



fn process_fasta(sequence_path: &str, monomer_path: &str,outpath: &str) -> Result<(), Box<dyn Error>> {

    let file = File::open(sequence_path)?;
    let reader = fasta::Reader::new(file);

    info!("Processing {} -> {}",sequence_path,monomer_path);
    let kmers_in_monomer = utils::create_kmers_from_sat(monomer_path).unwrap();
    

    reader.records().par_bridge().for_each(|result| {
        
        let record = result.expect("Error during fasta record parsing");
        
        let kmers = Kmer32::kmers_from_ascii(record.seq());

        let mut dist_mat: ArrayBase<OwnedRepr<u32>,
         Dim<[usize; 2]>> = Array2::zeros((kmers.len(), kmers.len()));
        let mut kmer_in_array_pos_dist: Vec<i32> = Vec::new();

        for i in 0..kmers.len() {

         let mut dist_vec: Vec<u32> = Vec::new();

         for z in 0..kmers_in_monomer.len()-1 {
            
            let dist = kmers[i].hamming_dist(kmers_in_monomer[z]);
            dist_vec.push(dist)

         }
        
         let min = dist_vec.iter().cloned().min().unwrap();
         kmer_in_array_pos_dist.push(min as i32)
         
        }

        let roll_mean: Vec<f64> = utils::calculate_means_around_index(&kmer_in_array_pos_dist);

        let _ = plot_roll_mean(&record,roll_mean.clone(),outpath);
    
     
        //writing the kmer tables, both the pos in array and the roll mean
        let outf = outpath;
        let outfn = outf.to_owned()+ "/data/" + &record.id().to_string().to_owned() + "_kmers_in_mono.txt";
        let mut file = File::create(outfn).unwrap();
        writeln!(file, "index\tactual\troll_mean").unwrap();
        for (i, (elem1, elem2)) in kmer_in_array_pos_dist.iter().zip(roll_mean.clone().iter()).enumerate() {
            writeln!(file, "{}\t{}\t{}", i, elem1, elem2).unwrap();
        }
        

    }
);
    Ok(())

}


fn init_logger() {
    // Read the RUST_LOG environment variable or use a default log level
    let log_level = std::env::var("RUST_LOG").unwrap_or_else(|_| String::from("info"));

    // Initialize the logger with the specified log level
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::from_str(&log_level).unwrap())
        .init();
}


fn remove_all_files_and_folders_in_folder(folder_path: &str) -> std::io::Result<()> {
    // Read directory entries and remove each file or subfolder
    for entry in fs::read_dir(folder_path)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            fs::remove_file(&path)?;
            println!("Deleted file: {:?}", path);
        } else if path.is_dir() {
            fs::remove_dir_all(&path)?;
            println!("Deleted folder: {:?}", path);
        }
    }

    Ok(())
}

fn main() {
    init_logger();

    let folder_path: &str = "./results/kmer_analysis";
    let data_path: &str = "./results/kmer_analysis/data/";
    let pic_path: &str = "./results/kmer_analysis/pictures/";
    

    // Create the folder if it doesn't exist
    if let Err(err) = fs::create_dir(folder_path) {
        if err.kind() != std::io::ErrorKind::AlreadyExists {
            eprintln!("Error creating folder: {:?}", err);
            return;
        }
    }

    // Check if the folder is not empty
    let is_empty = fs::read_dir(folder_path)
        .map(|entries| entries.count() == 0)
        .unwrap_or(true);

    if !is_empty {
        // Delete all files in the folder
        if let Err(err) = remove_all_files_and_folders_in_folder(&folder_path) {
            eprintln!("Error deleting files: {:?}", err);
            return;
        }
        info!("All files in the folder have been deleted.");
    } else {
        println!("The folder is empty.");
    }


    // create new data and folder paths
    if let Err(err) = fs::create_dir(data_path) {
        if err.kind() != std::io::ErrorKind::AlreadyExists {
            eprintln!("Error creating folder: {:?}", err);
            return;
        }
    }    
    
    if let Err(err) = fs::create_dir(pic_path) {
        if err.kind() != std::io::ErrorKind::AlreadyExists {
            eprintln!("Error creating folder: {:?}", err);
            return;
        }
    }    


    let json_file = "pairs.json";

    // Read the JSON file
    let pairs = utils::read_json_file(json_file);

    // Print the Monomer-Array Pairs
    info!(
        "Doing edge finding for the following RU:Array pairs: \n {}",
        pairs.iter()
            .map(|(monomer, array)| format!("{} -> {}", monomer, array))
            .collect::<Vec<String>>()
            .join("\n")
    );
    pairs.iter().for_each(|(monomer, array)| {
        let _ = process_fasta(array,monomer, folder_path);
    });

}
