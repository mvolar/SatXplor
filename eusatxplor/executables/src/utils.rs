use bio;
use serde_json::Error;

use std::fs::File;
use bio::io::fasta;
use debruijn::*;
use debruijn::kmer::*;
use std::io::{self, Write};
use std::collections::HashMap;
use std::io::BufReader;
use serde::{Deserialize, Serialize};
use std::io::{Read};


//function creates a hash table of all kmers in a sattellite, will be expanded into iterator over 
// multiple fasta files
pub fn create_kmers_from_sat(mononomer_path: &str) -> Result<Vec<IntKmer<u64>>, std::io::Error>{


    let fasta_file = File::open( mononomer_path)?;

    let mut kmer_total: Vec<IntKmer<u64>> = Vec::new();

    let reader = fasta::Reader::new(fasta_file);
    for result in reader.records() {
        
        let record = result.expect("Error during fasta record parsing");

        let newseq = Vec::from_iter(record.seq().iter().cloned().chain(record.seq().iter().cloned()));
        let k32_tmp = Kmer32::kmers_from_ascii(&newseq);

        'outer: for i in k32_tmp.iter().to_owned() {
            for j in kmer_total.iter().to_owned(){

                if i.hamming_dist(*j)==0 {
                    
                    break 'outer;
                }
            }
            kmer_total.push(*i)
            

        }
        }

        
        println!("\nNumber of unique kmers: {} in {}\n",kmer_total.len(),mononomer_path);
        return Ok(kmer_total)
    }
    

    
#[derive(Debug, Deserialize, Serialize)]
struct MonomerArrayPair {
        monomer_path: String,
        array_path: String,
}
    
    
pub fn get_monomer_array_pairs() -> Result<HashMap<String,String>,std::io::Error> {
        let json_file = "pairs.json";
    
        // Read the JSON file
        let pairs = read_json_file(json_file);
    
        // Print the Monomer-Array Pairs
        println!("Monomer-Array Pairs:");
        for (monomer_path, array_path) in &pairs {
            println!("{} -> {}", monomer_path, array_path);
        }
        return Ok(pairs)
    }


pub fn calculate_means_around_index(data: &[i32]) -> Vec<f64> {
        let mut means = Vec::new();
    
        for i in 10..(data.len() - 10) {
            let sum: i32 = data[i - 5..=i + 10].iter().cloned().sum();
            let count = 21.0; // Count of elements in the range [i - 5, i + 5]
            let mean = f64::from(sum) / count;
    
            means.push(mean);
        }
    
        means
    }



pub fn read_json_file(json_file: &str) -> HashMap<String, String> {
        // Open the JSON file
        let file = File::open(json_file).expect("Failed to open JSON file");
        let reader = BufReader::new(file);
    
        // Deserialize the JSON content into a HashMap<String, String>
        let pairs: HashMap<String, String> = serde_json::from_reader(reader).expect("Failed to deserialize JSON");
    
        pairs
    }

